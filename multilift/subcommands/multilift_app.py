import argparse
from pathlib import Path
import re
import tarfile
from time import sleep

import streamlit as st

from multilift import __prog_string__
from multilift.utils.state import MultiliftState


# Globals #####################################################################


parser = argparse.ArgumentParser()
parser.add_argument('--cache', type=str)


# Functions ###################################################################


def _session_create() -> None:
    args = parser.parse_args()
    if not 'dir_cache' in st.session_state:
        st.session_state.dir_cache = args.cache
    dir_working = Path(args.cache, st.session_state.session_id)
    (dir_working / 'uploads').mkdir(parents=True, exist_ok=True)


def _session_init() -> None:
    st.session_state.dir_working = \
        Path(st.session_state.dir_cache, st.session_state.session_id)
    st.session_state.multiliftstate = \
        MultiliftState(
            st.session_state.dir_working / 'config.yml', dump_on_modify=True)


def _callback_setter_session_id() -> None:
    if bool(re.match(
            r'^\w{8}-\w{4}-\w{4}-\w{4}-\w{12}$',
            new_token := st.session_state.setter_session_id)):
        if Path(st.session_state.dir_cache, new_token).exists():
            st.report_thread.get_report_ctx().session_id = new_token
            st.session_state.session_id = new_token
            _session_init()
            _callback_update_required()
        else:
            with container_message:
                st.error(
                    f'Session token not found within the cache: {new_token}')
                sleep(3)
            container_message.empty()
    else:
        with container_message:
            st.error(f'Not a valid session token: {new_token}')
            sleep(3)
        container_message.empty()
    st.experimental_rerun()


def _callback_update_required() -> None:
    st.session_state.display_viewer = False
    st.session_state.display_update = True


def _callback_add_genome() -> None:
    if not multiliftstate.genomes:
        multiliftstate.add_genome(
            st.session_state.add_genome, as_reference=True)
    else:
        multiliftstate.add_genome(st.session_state.add_genome)
    st.session_state.add_genome = ''
    _callback_update_required()


def _callback_genome_picker() -> None:
    removed = set(multiliftstate.genomes) - set(st.session_state.genome_picker)
    multiliftstate.del_genome(removed.pop())
    _callback_update_required()


def _callback_reference_picker() -> None:
    multiliftstate.add_genome(
        st.session_state.reference_picker, as_reference=True)
    _callback_update_required()


def _callback_file_uploader(application: str) -> None:
    for file in st.session_state[f'uploader_{application}']:
        destination = st.session_state.dir_working / 'uploads' / file.name
        if not destination.exists():
            with open(destination, 'wb') as f:
                f.write(file.getbuffer())
            multiliftstate.add_file(destination, application=application)
    _callback_update_required()


def _callback_file_assigner(genome: str) -> None:
    if (added :=
            set(st.session_state[f'file_assigner_{genome}'])
            - set(multiliftstate.files(genome=genome))):
        multiliftstate[added.pop(), 'genome'] = genome
    else:
        removed = set(multiliftstate.files(genome=genome)) \
            - set(st.session_state[f'file_assigner_{genome}'])
        del(multiliftstate[removed.pop(), 'genome'])
    _callback_update_required()


def _callback_update() -> None:
    with container_message:
        st.info('Working ...')
        sleep(3)
        st.info('Still working ...')
        sleep(3)
        st.success('Finished!')
        sleep(1)
    container_message.empty()
    st.session_state.display_viewer = True
    st.session_state.display_update = False


def _prepare_download() -> None:
    file = \
        st.session_state.dir_working / f'{st.session_state.session_id}.tar.gz'
    with tarfile.open(file, 'w:gz') as F:
        F.add(
            st.session_state.dir_working / 'uploads',
            arcname=f'{st.session_state.session_id}/uploads')
        F.add(
            st.session_state.dir_working / 'config.yml',
            arcname=f'{st.session_state.session_id}/config.yml')
    return file


# Streamlit: state configuration ##############################################


st.set_page_config(page_title=__prog_string__, layout='wide')

if 'session_id' not in st.session_state:
    st.session_state.session_id = st.report_thread.get_report_ctx().session_id
    _session_create()
    _session_init()
multiliftstate = st.session_state.multiliftstate

if 'display_viewer' not in st.session_state:
    st.session_state.display_viewer = False

if 'display_update' not in st.session_state:
    st.session_state.display_update = False


# Streamlit: app ##############################################################
# Everything below here renders!


st.title(f'{__prog_string__}')

container_message = st.empty()

container_viewer = st.container()
if st.session_state.display_viewer:
    with container_viewer:
        st.write('The viewer will display in here!')

container_picker = st.expander(
    'Genomes & Files', not st.session_state.display_viewer)
with container_picker:
    container_picker_main_columns = st.columns((1, 1, 1))

    with container_picker_main_columns[0]:
        st.text_input(
            'Enter genome name',
            key='add_genome',
            on_change=_callback_add_genome)
        st.multiselect(
            'Available genomes',
            key='genome_picker',
            options=multiliftstate.genomes,
            default=multiliftstate.genomes,
            on_change=_callback_genome_picker)
        st.selectbox(
            'Reference genome',
            key='reference_picker',
            options=multiliftstate.genomes,
            index=0,
            on_change=_callback_reference_picker)

    with container_picker_main_columns[1]:
        st.file_uploader(
            'Upload sequences ...',
            key='uploader_sequence',
            accept_multiple_files=True,
            help='''If genomes have multiple sequences, they should be
                consistently named in all files. Accepted formats: FASTA,
                GenBank, EMBL, SnapGene''',
            on_change=_callback_file_uploader,
            args=('sequence', ))
        st.file_uploader(
            '... or pre-computed alignments',
            key='uploader_alignment',
            accept_multiple_files=True,
            help='''If genomes have multiple sequences, each file's name
                indicates the name of the sequence and, within, sequence names
                in the alignment should match the genome names input here.
                Accepted MSA formats: FASTA, Clustal, MAF, Nexus, Phylip,
                Stockholm''',
            on_change=_callback_file_uploader,
            args=('alignment', ))

    with container_picker_main_columns[2]:
        st.file_uploader(
            'Upload data',
            key='uploader_data',
            accept_multiple_files=True,
            on_change=_callback_file_uploader,
            args=('data', ))

    if multiliftstate.genomes:
        st.markdown('---')
    for genome in multiliftstate.genomes:
        st.multiselect(
            f'{genome} files',
            key=f'file_assigner_{genome}',
            options=multiliftstate.files(
                genome=('', genome), application=('sequence', 'data')),
            default=multiliftstate.files(genome=genome),
            on_change=_callback_file_assigner,
            args=(genome, ))

    if st.session_state.display_update:
        st.markdown('---')
        st.button(
            'Update',
            on_click=_callback_update)


container_session = st.expander('Session', False)
with container_session:
    container_session_columns = st.columns((1, 1, 1))

    with container_session_columns[0]:
        st.text_input(
            'Session token',
            key='setter_session_id',
            value=st.session_state.session_id,
            help="Resume a previous session by entering it's token",
            on_change=_callback_setter_session_id)
        if not st.session_state.display_update:
            prepare_download = st.button(
                'Prepare session snapshot',
                help='Prepare a .tar.gz file of the current session state')
            if prepare_download:
                with open(_prepare_download(), 'rb') as file:
                    st.download_button(
                        'Download',
                        data=file,
                        file_name=f'{st.session_state.session_id}.tar.gz',
                        mime='application/octet-stream')

st.write(
    '[github.com/dbauerlab/multilift](https://github.com/dbauerlab/multilift)')


###############################################################################

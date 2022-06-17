from collections import defaultdict as dd
from itertools import chain, cycle
from io import StringIO, BytesIO
from pathlib import Path, PurePath
import tarfile

from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord

import streamlit as st

from multilift import __prog__, __prog_string__, __website__
from multilift.st_utils import message, v_space
from multilift.utils import sniff_filetype
from multilift.msa import align, test_aligners
from multilift.liftover import Lifter, liftover


###############################################################################
# Functions ###################################################################
###############################################################################


def refresh_ui(display_level: int=0, clear_message: bool=True) -> None:
    ''' Clear the messaging system and set the display level '''
    if clear_message:
        message()
    state.display_level = display_level


def session_init() -> None:
    ''' Initiate the streamlit session state '''
    state.session_id = \
        st.experimental_get_query_params().get(
            'session_id', [st.scriptrunner.get_script_run_ctx().session_id])[0]
    state.available_aligners = test_aligners()
    # a list of genome names
    state.multilift_genomes = []
    # the current reference genome
    state.multilift_reference = ''
    refresh_ui()


###############################################################################
# Callbacks ###################################################################
###############################################################################


def callback_add_genome() -> None:
    if not state.multilift_genomes:
        state.multilift_reference = state.uiobj_add_genome
    state.multilift_genomes.append(state.uiobj_add_genome)
    state.uiobj_add_genome = ''
    refresh_ui(0)


def callback_del_genome() -> None:
    for genome in set(state.multilift_genomes) - set(state.uiobj_genome_picker):
        state.multilift_genomes.remove(genome)
        try:
            del(state[f'uiobj_uploader_{genome}'])
        except KeyError:
            pass
    refresh_ui(0)


def callback_set_reference() -> None:
    state.multilift_reference = state.uiobj_reference_picker
    refresh_ui(0)


def callback_file_uploader() -> None:
    success = True
    # a dict of sequences {[file, seq.id, mapping, key]: seq, }
    state.multilift_sequences = {}
    for key in ['alignment'] + state.multilift_genomes:
        for file in state[f'uiobj_uploader_{key}']:
            ftype, application = sniff_filetype(file.name)
            if ftype == '':
                message(f'Unknown or unsupported filetype: {file.name}', 3)
                success = False
                break
            if key == 'alignment':
                if 'alignment' not in application:
                    message(f'{file.name} is not an alignment file', 3)
                    success = False
                    break
                with StringIO(file.getvalue().decode('utf-8')) as F:
                    state.multilift_sequences.update(
                        {(file.name, seq.id, None, key): seq
                        for aln in AlignIO.parse(F, ftype)
                        for seq in aln})
            elif set(('alignment', 'sequence')) & set(application):
                if bool(state.uiobj_uploader_alignment):
                    message(f'Sequences provided as well as alignments: {file.name}', 3)
                    success = False
                    break
                with StringIO(file.getvalue().decode('utf-8')) as F:
                    state.multilift_sequences.update(
                        {(file.name, seq.id, None, key): seq
                        for seq in SeqIO.parse(F, ftype)})
        if not success:
            break
    refresh_ui(2 if success else 1, success)


def callback_assign_sequence(mapping: str, current: list[tuple]) -> None:
    new = set(state[f'uiobj_sequence_assigner_{mapping}'])
    current = set(current)
    try:
        k = (new - current).pop()
        seq = state.multilift_sequences.pop(k)
        state.multilift_sequences[(k[0], k[1], mapping, k[3])] = seq
    except KeyError:
        k = (current - new).pop()
        seq = state.multilift_sequences.pop(k)
        state.multilift_sequences[(k[0], k[1], None, k[3])] = seq
    refresh_ui(2)


def callback_run_multilift() -> None:
    state.multilift_download = BytesIO()
    with tarfile.open(fileobj=state.multilift_download, mode='w:gz') as Tar:

        if not bool(state.uiobj_uploader_alignment):
            with container_message, st.spinner('Aligning sequences'):
                alignments = []
                for ref_fname, ref_seqid, _, _ in (
                        k for k in state.multilift_sequences.keys()
                        if k[3] == state.multilift_reference):
                    with StringIO() as F:
                        for k in state.multilift_sequences.keys():
                            fname, seqid, mapping, genome = k
                            if (fname == ref_fname
                                        and seqid == ref_seqid
                                        and genome == state.multilift_reference) \
                                    or mapping == ref_seqid:
                                SeqIO.write(
                                    SeqRecord(
                                        id=genome, description=seqid,
                                        seq=state.multilift_sequences[k].seq),
                                    F, 'fasta')
                        returncode, result = align(F, state.uiobj_aligner)
                        if returncode:
                            message(
                                f'Error making {state.uiobj_aligner} alignment for '
                                    f'[{ref_fname}] {ref_seqid} mappings',
                                3)
                            return
                        alignments.append(
                            {s.id:s for s in SeqIO.parse(result, 'fasta')})
                        tar_data = BytesIO(bytes(result.getvalue(), 'utf-8'))
                        tar_info = tarfile.TarInfo(
                            f'{state.session_id}/alignment/{ref_seqid}.fa')
                        tar_info.size = len(tar_data.getbuffer())
                        Tar.addfile(tar_info, tar_data)

        with container_message, st.spinner('Calculating liftovers'):
            L = Lifter()
            if not bool(state.uiobj_uploader_alignment):
                for aln in alignments:
                    L.add_alignment(aln, state.multilift_reference)
            else:
                for fname in set(k[0] for k in state.multilift_sequences.keys()):
                    L.add_alignment(
                        {k[2]: SeqRecord(id=k[2], description=v.id, seq=v.seq)
                            for k, v in state.multilift_sequences.items()
                            if k[0] == fname},
                        state.multilift_reference)

        with container_message, st.spinner('Performing liftovers'):
            state.multilift_liftovers = dd(dict)
            for genome in state.multilift_genomes:
                for file in state[f'uiobj_uploader_{genome}']:
                    ftype, application = sniff_filetype(file.name)
                    if 'data' not in application:
                        continue
                    try:
                        ext, lift_file = \
                            liftover(
                                StringIO(file.getvalue().decode('utf-8')),
                                ftype, L, genome)
                        tar_data = BytesIO(bytes(lift_file.getvalue(), 'utf-8'))
                        tar_info = tarfile.TarInfo(
                            f'{state.session_id}/liftover/{genome}/{PurePath(file.name).stem}.{ext}')
                        tar_info.size = len(tar_data.getbuffer())
                        Tar.addfile(tar_info, tar_data)
                    except:
                        message(f'Error lifting over {ref_fname}', 3)
                        return

    message('multilift ran sucessfully')
    refresh_ui(3, False)


###############################################################################
# Streamlit state configuration ###############################################
###############################################################################


state = st.session_state
if not 'session_id' in state:
    session_init()


###############################################################################
# Streamlit page layout #######################################################
###############################################################################
# Everything below here renders!

st.set_page_config(page_title=__prog__, layout='wide')

st.title(f'{__prog__}')

st.markdown(
    '<style> footer {visibility: hidden;} </style>',
    unsafe_allow_html=True)

# top-level messaging & errors
container_message = st.empty()
if state.message:
    with container_message:
        [st.success, st.info, st.warning, st.error][state.message_level](state.message)

container_0 = st.expander(
    'Define genomes',
    state.display_level == 0)

container_1 = st.expander(
    'Add / replace files',
    state.display_level == 1)

container_2 = st.expander(
    'Assign sequences',
    state.display_level == 2)

container_3 = st.container()

v_space(2, True)
container_session = st.expander('Session', False)

st.write(f'[{__prog_string__}]({__website__})')

###############################################################################
# Streamlit page functionality ################################################
###############################################################################


with container_0:
    container_0_cols = st.columns((1, 1, 1))

    with container_0_cols[0]:
        st.text_input(
            'Enter genome name',
            key='uiobj_add_genome',
            on_change=callback_add_genome,
            help='Add a genome to the current state')

    if state.multilift_genomes:

        with container_0_cols[1]:
            st.multiselect(
                'Available genomes',
                key='uiobj_genome_picker',
                options=state.multilift_genomes,
                default=state.multilift_genomes,
                on_change=callback_del_genome,
                help='Remove a genome to delete it from the current state')

        with container_0_cols[2]:
            st.selectbox(
                'Reference genome',
                key='uiobj_reference_picker',
                options=state.multilift_genomes,
                index=state.multilift_genomes.index(state.multilift_reference),
                on_change=callback_set_reference,
                help='Set a genome as the reference for the current state')

        if state.display_level == 0 and len(state.multilift_genomes) >= 2:
            st.button(
                'Continue',
                key='uiobj_update',
                on_click=refresh_ui,
                args=(1,))


###############################################################################


if state.display_level >= 1:

    with container_1:
        container_1_cols = st.columns((1, 1, 1))

        for genome, col in zip(state.multilift_genomes, cycle((0, 1, 2))):
            with container_1_cols[col]:
                st.file_uploader(
                    f'{genome}:',
                    key=f'uiobj_uploader_{genome}',
                    accept_multiple_files=True,
                    on_change=refresh_ui,
                    args=(1, False))

        with container_1_cols[0]:
            st.file_uploader(
                'Alignment(s) (optional):',
                key='uiobj_uploader_alignment',
                accept_multiple_files=True,
                on_change=refresh_ui,
                args=(1, False))

        if st.session_state.display_level == 1:
            st.button(
                'Continue',
                key='uiobj_update',
                on_click=callback_file_uploader)


###############################################################################


if state.display_level >= 2:

    with container_2:
        container_2_cols = st.columns((1, 1))

        if bool(state.uiobj_uploader_alignment):
            for genome, col in zip(state.multilift_genomes, cycle((0, 1))):
                with container_2_cols[col]:
                    current = \
                        [k for k in state.multilift_sequences.keys()
                        if k[2] == genome]
                    st.multiselect(
                        f'{genome}',
                        key=f'uiobj_sequence_assigner_{genome}',
                        options=[
                            k for k in state.multilift_sequences.keys()
                            if k[2] in (None, genome)],
                        default=current,
                        format_func=lambda x: f'[{x[0]}] {x[1]}',
                        on_change=callback_assign_sequence,
                        args=(genome, current))

        else:
            for ref_seqid, col in zip(
                    (seqid
                        for _, seqid, _, genome in state.multilift_sequences.keys()
                        if genome == state.multilift_reference),
                    cycle((0, 1))):
                with container_2_cols[col]:
                    current = \
                        [k for k in state.multilift_sequences.keys()
                        if k[2] == ref_seqid]
                    st.multiselect(
                        f'{ref_seqid}',
                        key=f'uiobj_sequence_assigner_{ref_seqid}',
                        options=[
                            k for k in state.multilift_sequences.keys()
                            if k[2] in (None, ref_seqid)
                            and k[3] is not state.multilift_reference],
                        default=current,
                        format_func=lambda x: f'[{x[0]}] {x[1]}',
                        on_change=callback_assign_sequence,
                        args=(ref_seqid, current))

        if st.session_state.display_level >= 2:
            st.button(
                'Run multilift',
                key='uiobj_update',
                on_click=callback_run_multilift)


###############################################################################


if state.display_level >= 3:
    with container_3:
       st.download_button(
            'Download multilift results',
            data=state.multilift_download,
            file_name=f'{state.session_id}.tar.gz',
            mime='application/octet-stream')


###############################################################################


with container_session:
    st.radio(
        'Alignment program',
        key='uiobj_aligner',
        options=state.available_aligners,
        horizontal=True)

###############################################################################

from collections import defaultdict as dd
from itertools import cycle
from io import StringIO, BytesIO
import tarfile
import zipfile

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

import streamlit as st

from multilift import __prog__, __prog_string__, __website__
from multilift.liftover import annotate, Lifter, liftover
from multilift.msa import align, aligner_limits, generate_consensus, test_aligners
from multilift.utils import add_to_archive, basename, create_igv_session, sniff_filetype


###############################################################################
# Functions ###################################################################
###############################################################################

def v_space(n: int=1, hline: bool=False) -> None:
    ''' '''
    for i in range(n):
        st.write('')
    if hline:
        st.write('***')

def message(message: str='', level: int=0) -> None:
    ''' Display a message in the `container_message` area on script re-run.
    0: success, 1: info, 2: warning, 3: error '''
    state.message = message
    state.message_level = level


def refresh_ui(display_level: int=0, clear_message: bool=True) -> None:
    ''' Clear the messaging system and set the display level '''
    if clear_message:
        message()
    state.display_level = display_level


def session_init() -> None:
    ''' Initiate the streamlit session state '''
    state.session_id = \
        st.experimental_get_query_params().get(
            'session_id',
            [st.runtime.scriptrunner.script_run_context.get_script_run_ctx().session_id])[0]
    state.available_aligners = test_aligners()
    state.multilift_genomes = []
    state.multilift_seq_groups = ['Group1', ]
    refresh_ui()


###############################################################################
# Callbacks ###################################################################
###############################################################################


def callback_add_genome() -> None:
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


def callback_file_uploader() -> None:
    clear_message = success = True
    state.multilift_sequences = {}  # {(key, file, seq.id, group): seq, }
    for key in ['alignment'] + state.multilift_genomes:
        for file in state[f'uiobj_uploader_{key}']:
            ftype, application = sniff_filetype(file.name)
            if ftype == '':
                message(f'Unknown or unsupported filetype: {file.name}', 3)
                clear_message = success = False
                break
            if key == 'alignment':
                if 'alignment' not in application:
                    message(f'{file.name} is not an alignment file', 3)
                    clear_message = success = False
                    break
                try:
                    with StringIO(file.getvalue().decode('utf-8')) as F:
                        state.multilift_sequences.update(
                            {(key, file.name, seq.id, None): seq
                            for seq in AlignIO.read(F, ftype)})
                except:
                    message(f'Error reading sequences: {file.name}', 3)
                    clear_message = success = False
                    break
            elif 'sequence' in application:
                if bool(state.uiobj_uploader_alignment):
                    message(f'Sequences provided as well as alignments: {file.name}', 3)
                    clear_message = success = False
                    break
                try:
                    with StringIO(file.getvalue().decode('utf-8')) as F:
                        state.multilift_sequences.update(
                            {(key, file.name, seq.id, None): seq
                            for seq in SeqIO.parse(F, ftype)})
                except:
                    message(f'Error reading sequences: {file.name}', 3)
                    clear_message = success = False
                    break
        if not success:
            break
    # check sequence inputs & subset aligners
    if not bool(state.uiobj_uploader_alignment):
        if not set((k[0] for k in state.multilift_sequences.keys())) == \
                set(state.multilift_genomes):
            message(f'All genomes must have sequences provided', 3)
            clear_message = success = False
        else:
            longest_seq = \
                max(len(s) for s in state.multilift_sequences.values())
            for aligner, limit in aligner_limits.items():
                if longest_seq > limit:
                    try:
                        state.available_aligners.remove(aligner)
                    except KeyError:
                        pass
            if not state.available_aligners:
                message(
                    'No installed aligners are capable of aligning sequences of '
                    'this length. Please upload pre-computed alignments.',
                    2)
                clear_message = success = False
    refresh_ui(2 if success else 1, clear_message)


def callback_sequence_assigner_name(i: int) -> None:
    old_name = state.multilift_seq_groups[i]
    new_name = state[f'uiobj_sequence_assigner_{i}_name']
    if not new_name:
        message('All sequence groups must have names', 3)
        return
    if new_name in state.multilift_seq_groups:
        message('All sequence groups must have different names', 3)
        return
    state.multilift_seq_groups[i] = new_name
    refresh_ui(2)


def callback_addremove_seq_group(add_group: bool=True) -> None:
    if add_group:
        state.multilift_seq_groups.append(
            f'Group{len(state.multilift_seq_groups) + 1}')
    else:
        drop_idx = str(len(state.multilift_seq_groups) - 1)
        state.multilift_seq_groups.pop()
        state.multilift_sequences = {
            (k[0], k[1], k[2], None if k[3] == drop_idx else k[3]): v
            for k, v in state.multilift_sequences.items()}
        del(state[f'uiobj_sequence_assigner_{i}_name'])
        del(state[f'uiobj_sequence_assigner_{i}'])
    refresh_ui(2)


def callback_assign_sequence(group: str, current: list[tuple]) -> None:
    new = set(state[f'uiobj_sequence_assigner_{group}'])
    current = set(current)
    try:
        k = (new - current).pop()
    except KeyError:
        k = (current - new).pop()
        group = None
    seq = state.multilift_sequences.pop(k)
    state.multilift_sequences[(k[0], k[1], k[2], group)] = seq
    refresh_ui(2)


def callback_run_multilift() -> None:
    state.multilift_download = BytesIO()
    maf_alignments = []
    igv_resources = []
    igv_genomes = []
    L = Lifter()

    with \
            tarfile.open(fileobj=state.multilift_download, mode='w:gz') \
            if state.uiobj_download_format == '.tar.gz' else \
            zipfile.ZipFile(state.multilift_download, 'w') as Arc:

        if not bool(state.uiobj_uploader_alignment):
            for i, seq_group in enumerate(state.multilift_seq_groups):
                with container_message, \
                        st.spinner(f'Aligning sequence(s): {seq_group}'):
                    with StringIO() as F:
                        for (genome, fname, seqid, group), seq \
                                in state.multilift_sequences.items():
                            if group == str(i):
                                SeqIO.write(
                                    SeqRecord(
                                        id=f'{genome} {seqid}',
                                        description='',
                                        seq=seq.seq),
                                    F, 'fasta')
                        returncode, result = align(F, state.uiobj_aligner)
                    if returncode:
                        message(
                            f'Error making {state.uiobj_aligner} alignment '
                            f'for sequence group "{seq_group}"',
                            3)
                        return
                    aln = AlignIO.read(result, 'fasta')
                    L.add_alignment(aln, seq_group)
                    # write alignment as fasta
                    add_to_archive(
                        Arc,
                        BytesIO(bytes(result.getvalue(), 'utf-8')),
                        f'{state.session_id}/alignment/{seq_group}.fa',
                        state.uiobj_download_format)
                    # store consensus
                    cons = Seq(generate_consensus(aln))
                    igv_genomes.append(
                        SeqRecord(
                            id=seq_group, description='',
                            seq=cons))
                    # add to maf
                    maf_alignments.append(
                        MultipleSeqAlignment(
                            [SeqRecord(
                                id=f'multilift.{seq_group}', description='',
                                seq=cons)] +
                            [SeqRecord(
                                id=f'{s.id}.{seq_group}', description='',
                                seq=s.seq)
                            for s in aln]))

        else:  # alignments provided
            with container_message, st.spinner('Reading alignment(s)'):
                for fname in set(k[1] for k in state.multilift_sequences.keys()):
                    aln = MultipleSeqAlignment(
                        [SeqRecord(
                            id=f'{k[3]} {k[2]}', description='', seq=v.seq)
                        for k, v in state.multilift_sequences.items()
                        if k[1] == fname])
                    L.add_alignment(aln, basename(fname))
                    # store consensus
                    cons = Seq(generate_consensus(aln))
                    igv_genomes.append(
                        SeqRecord(
                            id=basename(fname), description='',
                            seq=cons))
                    # add to maf
                    maf_alignments.append(
                        MultipleSeqAlignment(
                            [SeqRecord(
                                id=f'multilift.{basename(fname)}',
                                description='',
                                seq=cons)] +
                            [SeqRecord(
                                id=f'{s.id}.{basename(fname)}',
                                description='',
                                seq=s.seq)
                            for s in aln]))

        # write consensus sequences
        with StringIO() as F:
            SeqIO.write(igv_genomes, F, 'fasta')
            add_to_archive(
                Arc,
                BytesIO(bytes(F.getvalue(), 'utf-8')),
                f'{state.session_id}/genome/multilift_{state.session_id}.fa',
                state.uiobj_download_format)
        del(igv_genomes)

        # write maf alignment
        with StringIO() as F:
            AlignIO.write(maf_alignments, F, 'maf')
            add_to_archive(
                Arc,
                BytesIO(bytes(F.getvalue(), 'utf-8')),
                f'{state.session_id}/alignment/multilift_{state.session_id}.maf',
                state.uiobj_download_format)
            igv_resources.append(
                f'alignment/multilift_{state.session_id}.maf')

        # compute maf coordinates
        for idx in range(len(maf_alignments)):
            aln = maf_alignments[idx]
            coords = [
                SeqRecord(
                    id=aln[0].id, description='',
                    seq=Seq('.' * aln.get_alignment_length()))]
            for s in aln[1:]:
                seq = list(s.seq)
                i = 0
                for j, nt in enumerate(s.seq):
                    if nt != '-':
                        seq[j] = '|' if i % 10 == 9 else '.'
                        i += 1
                i = 1
                for j, nt in enumerate(seq[:]):
                    if nt == '|':
                        coord_str = str(i * 10)
                        if '-' not in seq[j:j+len(coord_str)] \
                                and j+len(coord_str) < len(seq):
                            seq[j:j+len(coord_str)] = list(coord_str)
                        i += 1
                coords.append(
                    SeqRecord(
                        id=s.id, description='',
                        seq=Seq(''.join(seq))))
            maf_alignments[idx] = MultipleSeqAlignment(coords)

        # write coordinates file
        with StringIO() as F:
            AlignIO.write(maf_alignments, F, 'maf')
            add_to_archive(
                Arc,
                BytesIO(bytes(F.getvalue(), 'utf-8')),
                f'{state.session_id}/alignment/multilift_{state.session_id}.coords.maf',
                state.uiobj_download_format)
        igv_resources.append(
            f'alignment/multilift_{state.session_id}.coords.maf')
        del(maf_alignments)

        # liftover data files, create annotations
        with container_message, st.spinner('Performing liftovers'):
            for genome in state.multilift_genomes:
                for file in state[f'uiobj_uploader_{genome}']:
                    ftype, application = sniff_filetype(file.name)
                    if 'data' in application:
                        try:
                            new_ext, lift_file = \
                                liftover(
                                    StringIO(file.getvalue().decode('utf-8')),
                                    ftype, L, genome)
                        except Exception as e:
                            message(
                                f'Error lifting over {file.name} for {genome}. '
                                    f'Job failed with: {e}',
                                3)
                            return
                    elif 'annotation' in application:
                        try:
                            new_ext, lift_file = \
                                annotate(
                                    StringIO(file.getvalue().decode('utf-8')),
                                    L,
                                    {k[2]: v
                                        for k, v in state.multilift_sequences.items()
                                        if k[0] == genome},
                                    genome)
                        except Exception as e:
                            message(
                                f'Error applying {file.name} annotations to {genome}. '
                                    f'Job failed with: {e}',
                                3)
                            return
                    else:
                        continue
                    add_to_archive(
                        Arc,
                        BytesIO(bytes(lift_file.getvalue(), 'utf-8')),
                        f'{state.session_id}/liftover/{genome}/{file.name}{new_ext}',
                        state.uiobj_download_format)
                    igv_resources.append(
                        f'liftover/{genome}/{file.name}{new_ext}')

        with container_message, st.spinner('Preparing IGV session'):
            add_to_archive(
                Arc,
                BytesIO(bytes(
                    create_igv_session(state.session_id, igv_resources),
                    'utf-8')),
                f'{state.session_id}/igv_session.xml',
                state.uiobj_download_format)

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
    '<style> '
    'footer {visibility: hidden;} '
    '#MainMenu {visibility: hidden;} '
    '</style>',
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
    container_0_cols = st.columns((1, 2))

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

        if state.display_level == 0 and len(state.multilift_genomes) >= 2:
            st.button(
                'Continue',
                key='uiobj_update',
                on_click=refresh_ui,
                args=(1,))


###############################################################################


if state.display_level >= 1:

    with container_1:

        for key, col in zip(
                state.multilift_genomes + ['alignment'],
                cycle((0, 1, 2))):
            if col == 0:
                container_1_cols = st.columns((1, 1, 1))
            with container_1_cols[col]:
                st.file_uploader(
                    f'{key}:' if key != 'alignment' else 'Alignment(s) (optional):',
                    key=f'uiobj_uploader_{key}',
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

        unassigned = \
            [k for k in state.multilift_sequences.keys() if k[3] is None]

        if bool(state.uiobj_uploader_alignment):
            container_2_cols = st.columns((1, 1))
            for genome, col in zip(state.multilift_genomes, cycle((0, 1))):
                with container_2_cols[col]:
                    current = \
                        [k for k in state.multilift_sequences.keys()
                        if k[3] == genome]
                    st.multiselect(
                        f'{genome}',
                        key=f'uiobj_sequence_assigner_{genome}',
                        options=unassigned + current,
                        default=current,
                        format_func=lambda x: f'[{x[1]}] {x[2]}',
                        on_change=callback_assign_sequence,
                        args=(genome, current))

        else:
            container_2_cols = st.columns((1, 2))
            for i, seq_group in enumerate(state.multilift_seq_groups):
                with container_2_cols[0]:
                    st.text_input(
                        '',
                        key=f'uiobj_sequence_assigner_{i}_name',
                        value=seq_group,
                        on_change=callback_sequence_assigner_name,
                        args=(i, ))
                with container_2_cols[1]:
                    current = \
                        [k for k in state.multilift_sequences.keys()
                        if k[3] == str(i)]
                    st.multiselect(
                        '',
                        key=f'uiobj_sequence_assigner_{i}',
                        options=unassigned + current,
                        default=current,
                        format_func=lambda x: f'[{x[1]}] {x[2]}',
                        on_change=callback_assign_sequence,
                        args=(str(i), current),
                        help='This value serves as the liftover chromosome name')

        container_2_cols = st.columns(6)
        if not bool(state.uiobj_uploader_alignment):
            with container_2_cols[4]:
                st.button(
                    'Add sequence group',
                    key='uiobj_sequence_group_add',
                    on_click=callback_addremove_seq_group)
            with container_2_cols[5]:
                st.button(
                    'Remove sequence group',
                    key='uiobj_sequence_group_remove',
                    on_click=callback_addremove_seq_group,
                    args=(False, ),
                    disabled=len(state.multilift_seq_groups) == 1)
        if unassigned:
            st.info(
                'Unassigned sequences: ' +
                ', '.join("[{1}] {2}".format(*k) for k in unassigned))
        elif st.session_state.display_level >= 2:
            with container_2_cols[0]:
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
            file_name=f'{state.session_id}{state.uiobj_download_format}',
            mime='application/octet-stream')


###############################################################################


with container_session:

    st.radio(
        'Alignment program',
        key='uiobj_aligner',
        options=state.available_aligners,
        index=\
            state.available_aligners.index('mafft')
            if 'mafft' in state.available_aligners else 0,
        horizontal=True)

    st.radio(
        'Download format',
        key='uiobj_download_format',
        options=['.tar.gz', '.zip'],
        index=0,
        horizontal=True)

    st.write(state)


###############################################################################

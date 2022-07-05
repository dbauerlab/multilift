from collections import Counter
from io import StringIO
import shlex
from subprocess import Popen, PIPE, run

from Bio.Align import MultipleSeqAlignment


aligner_limits = {
    'clustalo': 2_000,
    'kalign': 10_000,
    'muscle': 10_000,
    'mafft': 100_000}

IUPAC_NTS_TO_AMBIG = {
    'AG': 'R', 'CT': 'Y', 'CG': 'S', 'AT': 'W', 'GT': 'K', 'AC': 'M',
    'CGT': 'B', 'AGT': 'D', 'ACT': 'H', 'ACG': 'V',
    'ACGT': 'N'}
IUPAC_AMBIG_TO_NTS = {v: k for k, v in IUPAC_NTS_TO_AMBIG.items()}
IUPAC_AMBIG_TO_NTS.update({
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'U': 'T'})


def consensus_base(column: str, iupac: bool=False) -> str:
    ''' Decide the consensus at a base position '''
    column = ''.join(b for b in column if b != '-')
    if not column:  # Return gap only
        return '-'
    nts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    for base, count in Counter(column).items():
        for nt in IUPAC_AMBIG_TO_NTS[base]:
            nts[nt] += count
    nts = \
        sorted(
            [(k, v) for k, v in nts.items()],
            key=lambda x: x[1], reverse=True)
    if not iupac:
        nts = [nt for nt in nts if nt[1] == nts[0][1]]
        if len(nts) == 1:
            return nts[0][0]  # Return majority
        return 'N'  # Return majority indecision
    nts = [nt for nt in nts if nt[1] >= (nts[0][1] * 0.5)]
    if len(nts) == 1:
        return nts[0][0]  # Return majority
    return IUPAC_NTS_TO_AMBIG[''.join(sorted([nt[0] for nt in nts]))]  # Return iupac indecision


def generate_consensus(alignment: MultipleSeqAlignment) -> str:
    return ''.join(
        consensus_base(alignment[:, i].upper())
        for i in range(alignment.get_alignment_length()))


def test_aligners() -> set[str]:
    ''' '''
    aligner_versions = {
        'mafft': 'mafft --version',
        'kalign': 'kalign -v',
        'muscle': 'muscle -version',
        'clustalo': 'clustalo --version'}
    return set((
        k for k, v in aligner_versions.items()
        if run(shlex.split(v), capture_output=True).returncode == 0))


def align(file: StringIO, prog: str='mafft', threads: int=1) -> tuple[int, StringIO]:
    '''  '''
    if prog == 'mafft':
        with Popen(
                shlex.split(f'mafft --6merpair --nuc --nwildcard --thread {threads} -'),
                stdin=PIPE, stdout=PIPE, stderr=PIPE, text=True) as P:
            std_out, std_err = P.communicate(file.getvalue())
            if P.returncode:
                return P.returncode, StringIO(std_err)
            return 0, StringIO(std_out)

    if prog == 'clustalo':
        with Popen(
                shlex.split(f'clustalo -v -i - --threads {threads}'),
                stdin=PIPE, stdout=PIPE, stderr=PIPE, text=True) as P:
            std_out, std_err = P.communicate(file.getvalue())
            if P.returncode:
                return P.returncode, StringIO(std_err)
            return 0, StringIO(std_out)

    if prog == 'kalign':
        with \
                Popen(
                    ['kalign'],
                    stdin=PIPE, stdout=PIPE, stderr=PIPE, text=True) as P, \
                Popen(
                    shlex.split(f"sed -n '/^>/,$p'"),
                    stdin=P.stdout, stdout=PIPE, text=True) as P2:
            _ = P.stdin.write(file.getvalue())
            P.stdin.close()
            std_out, std_err = P2.communicate()
            if P.returncode:
                return P.returncode, StringIO(std_err)
            if P2.returncode:
                return P2.returncode, StringIO(std_err)
            return 0, StringIO(std_out)

    if prog == 'muscle':
        with Popen(
                shlex.split('muscle -maxiters 1 -diags'),
                stdin=PIPE, stdout=PIPE, stderr=PIPE, text=True) as P:
            std_out, std_err = P.communicate(file.getvalue())
            if P.returncode:
                return P.returncode, StringIO(std_err)
            return 0, StringIO(std_out)




###############################################################################

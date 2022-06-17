from io import StringIO
import shlex
from subprocess import Popen, PIPE, run


def test_aligners() -> list[str]:
    ''' '''
    aligner_versions = {
        'mafft': 'mafft --version',
        'kalign': 'kalign -v',
        'muscle': 'muscle -version',
        'clustalo': 'clustalo --version'}
    return [
        k for k, v in aligner_versions.items()
        if run(shlex.split(v), capture_output=True).returncode == 0]


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

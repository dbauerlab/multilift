import bz2
import gzip
import hashlib
import logging
from pathlib import PurePath

from Bio import SeqIO, AlignIO

from multilift import __prog__


# Globals ######################################################################


__all__ = ['SeqFile', 'AlnFile']

logger = logging.getLogger(__prog__)


################################################################################


class _SeqIOGeneric():

    def __init__(self, file: str) -> None:
        self.file = PurePath(file)
        self._exts = self.file.suffixes
        if '.gz' in self.exts:
            self._fileobj = gzip.open(self.file, 'rt')
            self.compressed = True
        elif '.bz2' in self.exts:
            self._fileobj = bz2.open(self.file, 'rt')
            self.compressed = True
        else:
            self._fileobj = open(self.file)
            self.compressed = False

    @property
    def md5(self) -> str:
        ''' The MD5 checksum hex digest of the file '''
        h = hashlib.md5()
        with open(self.file, 'rb') as f:
            for chunk in iter(lambda: f.read(4096), b''):
                h.update(chunk)
        return h.hexdigest()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        pass


class SeqFile(_SeqIOGeneric):
    ''' An extension of _SeqIOGeneric that wraps Bio.SeqIO, guesses the correct
    parser from the file extension, and provides an interator '''

    def __init__(self, file: str) -> None:
        super().__init__()
        if any(e in self._exts for e in ('.fa', '.fasta')):
            self.format = 'fasta'
        elif '.embl' in self._exts:
            self.format = 'embl'
        elif any(e in self._exts for e in ('.gb', '.genbank')):
            self.format = 'genbank'
        if '.dna' in self._exts:
            self.format = 'snapgene'
        else:
            logger.error(
                f'Cannot determine filetype for the extension(s) '
                f'{"".join(self._exts)}')

    def __iter__(self):
        for record in SeqIO.parse(self.file, self.format):
            yield record


class AlnFile(_SeqIOGeneric):
    ''' An extension of _SeqIOGeneric that wraps Bio.AlignIO, guesses the
    correct parser from the file extension, and provides an interator '''

    def __init__(self, file: str) -> None:
        super().__init__()
        if any(e in self._exts for e in ('.fa', '.fasta')):
            self.format = 'fasta'
        elif '.clustal' in self._exts:
            self.format = 'clustal'
        elif '.maf' in self._exts:
            self.format = 'maf'
        elif any(e in self._exts for e in ('.nex', '.nexus')):
            self.format = 'nexus'
        elif any(e in self._exts for e in ('.phy', '.phylip')):
            self.format = 'phylip'
        elif any(e in self._exts for e in ('.sto', 'sth', '.stockholm')):
            self.format = 'stockholm'
        else:
            logger.error(
                f'Cannot determine filetype for the extension(s) '
                f'{"".join(self._exts)}')

    def __iter__(self):
        for record in AlignIO.parse(self.file, self.format):
            yield record

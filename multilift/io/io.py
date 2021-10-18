import bz2
from collections.abc import Iterator
from functools import cached_property
import gzip
import hashlib
import logging
from pathlib import Path

from Bio import AlignIO, Entrez, SeqIO, SeqRecord

from multilift import __prog__


# Globals #####################################################################


__all__ = ['SeqFile', 'AlnFile', 'fetch']

logger = logging.getLogger(__prog__)


# Classes #####################################################################


class _BioIOGeneric():
    ''' PRIVATE. Base class for the SeqFile and AlnFile parser classes
    containing the shared functionality '''

    def __init__(self, file: str) -> None:
        self.file = Path(file)
        if not self.file.is_file():
            logger.error(f'Cannot find file: {self.file}')
        # invert extension list to respect precedence
        self._exts = [ext.lower() for ext in self.file.suffixes][::-1]
        if '.gz' in self._exts:
            self._fileobj = gzip.open(self.file, 'rt')
            self.compressed = True
        elif '.bz2' in self._exts:
            self._fileobj = bz2.open(self.file, 'rt')
            self.compressed = True
        else:
            self._fileobj = open(self.file)
            self.compressed = False

    @cached_property
    def md5(self) -> str:
        ''' The MD5 checksum hex digest of the file '''
        h = hashlib.md5()
        with open(self.file, 'rb') as f:
            for chunk in iter(lambda: f.read(4096), b''):
                h.update(chunk)
        return h.hexdigest()

    def __repr__(self) -> str:
        return \
            f'{self.__class__.__name__}(' \
            f'file={self.file}, ' \
            f'format={self.format}, ' \
            f'md5={self.md5})'

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self._fileobj.close()


class SeqFile(_BioIOGeneric):
    ''' An extension of _BioIOGeneric that wraps Bio.SeqIO, guesses the correct
    parser from the file extension, and provides an interator '''

    def __init__(self, file: str) -> None:
        super().__init__(file)
        self.format = self._guess_format()

    def _guess_format(self) -> str:
        ''' PRIVATE. Guess the filetype or raise an error via the log '''
        for ext in self._exts:
            if ext in ('.fa', '.faa', '.fasta'):
                return 'fasta'
            if ext == '.embl':
                return 'embl'
            if ext in ('.gb', '.genbank'):
                return 'genbank'
            if ext == '.dna':
                return 'snapgene'
        logger.error(
            f'Cannot determine filetype for the extension(s) '
            f'{"".join(self.file.suffixes)}')

    def __iter__(self) -> Iterator[SeqRecord]:
        ''' Iterate and yield Bio.SeqRecord objects '''
        self._fileobj.seek(0)  # make sure we're at the start of the file
        for record in SeqIO.parse(self._fileobj, self.format):
            yield record


class AlnFile(_BioIOGeneric):
    ''' An extension of _BioIOGeneric that wraps Bio.AlignIO, guesses the
    correct parser from the file extension, and provides an interator '''

    def __init__(self, file: str) -> None:
        super().__init__(file)
        self.format = self._guess_format()

    def _guess_format(self) -> str:
        ''' PRIVATE. Guess the filetype or raise an error via the log '''
        for ext in self._exts:
            if ext in ('.fa', '.faa', '.fasta', '.aln'):
                return 'fasta'
            if ext == '.clustal':
                return 'clustal'
            if ext == '.maf':
                return 'maf'
            if ext in ('.nex', '.nexus'):
                return 'nexus'
            if ext in ('.phy', '.phylip'):
                return 'phylip'
            if ext in ('.sto', 'sth', '.stockholm'):
                return 'stockholm'
        logger.error(
            f'Cannot determine filetype for the extension(s) '
            f'{"".join(exts)}')

    def __iter__(self) -> Iterator[SeqRecord]:
        ''' Iterate and yield Bio.SeqRecord objects '''
        self._fileobj.seek(0)  # make sure we're at the start of the file
        for record in AlignIO.read(self._fileobj, self.format):
            yield record


# Functions ###################################################################


def fetch(accession: str, email: str) -> SeqRecord:
    ''' Fetch a single nucleotide record from NCBI as GenBank and return an
    annotated Bio.SeqRecord '''
    Entrez.email = email
    handle = Entrez.efetch(
        db='nucleotide', id=accession, rettype='gb', retmode='text')
    return SeqIO.read(handle, 'genbank')


###############################################################################

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
        if '.gz' in self.file.suffixes:
            self._fileobj = gzip.open(self.file, 'rt')
            self.compressed = True
        elif '.bz2' in self.file.suffixes:
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
        exts = self.file.suffixes
        if any(e in exts for e in ('.fa', '.faa', '.fasta')):
            return 'fasta'
        if '.embl' in exts:
            return 'embl'
        if any(e in exts for e in ('.gb', '.genbank')):
            return 'genbank'
        if '.dna' in exts:
            return 'snapgene'
        logger.error(
            f'Cannot determine filetype for the extension(s) '
            f'{"".join(exts)}')

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
        exts = self.file.suffixes
        if any(e in exts for e in ('.fa', '.fasta')):
            return 'fasta'
        if '.clustal' in exts:
            return 'clustal'
        if '.maf' in exts:
            return 'maf'
        if any(e in exts for e in ('.nex', '.nexus')):
            return 'nexus'
        if any(e in exts for e in ('.phy', '.phylip')):
            return 'phylip'
        if any(e in exts for e in ('.sto', 'sth', '.stockholm')):
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

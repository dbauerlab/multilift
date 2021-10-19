from collections.abc import Iterator
from functools import cached_property
import logging
from pathlib import Path

from Bio import AlignIO, Entrez, SeqIO, SeqRecord

from multilift import __prog__
from multilift.utils import \
    open_helper, file_hash, filetype_associations, supported_filetypes


# Globals #####################################################################


logger = logging.getLogger(__prog__)


# Classes #####################################################################


class _BioIOGeneric():
    ''' PRIVATE. Base class for the SeqFile and AlnFile parser classes
    containing the shared functionality '''

    def __init__(self, file: str, filetype: str) -> None:
        self.file = Path(file)
        self.filetype = filetype
        if not self.file.is_file():
            logger.error(f'Cannot find file: {self.file}')

    def __repr__(self) -> str:
        return \
            f'{self.__class__.__name__}(' \
            f'file={self.file}, ' \
            f'filetype={self.filetype}, ' \
            f'md5={self.md5})'

    @cached_property
    def md5(self) -> str:
        return file_hash(self.file)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        pass


class SeqFile(_BioIOGeneric):
    ''' An extension of _BioIOGeneric that wraps Bio.SeqIO, guesses the correct
    parser from the file extension, and provides an interator '''

    def __init__(self, file: str, filetype: str) -> None:
        super().__init__(file, filetype)

    def __iter__(self) -> Iterator[SeqRecord]:
        ''' Iterate and yield Bio.SeqRecord objects '''
        with open_helper(self.file) as F:
            for record in SeqIO.parse(F, self.filetype):
                yield record


class AlnFile(_BioIOGeneric):
    ''' An extension of _BioIOGeneric that wraps Bio.AlignIO, guesses the
    correct parser from the file extension, and provides an interator '''

    def __init__(self, file: str) -> None:
        super().__init__(file, filetype)

    def __iter__(self) -> Iterator[SeqRecord]:
        ''' Iterate and yield Bio.SeqRecord objects '''
        with open_helper(self.file) as F:
            for record in AlignIO.parse(F, self.filetype):
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

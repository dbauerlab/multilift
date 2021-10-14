from bisect import bisect_left
import logging

from Bio import SeqRecord

from multilift import __prog__


# Globals #####################################################################


__all__ = ['Lifter']

logger = logging.getLogger(__prog__)


# Classes #####################################################################


class Lifter():
    ''' A class that ingests Bio.SeqRecord objects from MSAs and returns
    lifted-over positions when called '''

    def __init__(self) -> None:
        self._gaps = {}

    def add_liftover(self, genome: str, chrom: str, record: SeqRecord) -> None:
        ''' Calculate the liftover table for a sequence by recording the
        location (and hence number) of alignment gaps relative to the ungapped
        sequence string '''
        key = f'{genome}_{chrom}'
        self._gaps[key] = []
        last_nt = -1
        for i, nt in enumerate(record.seq):
            if nt != '-':
                last_nt = i
            else:
                self._gaps[key].append(last_nt)

    def __call__(self, genome: str, chrom: str, pos: int) -> int:
        ''' Liftover `genome_chrom:pos` into alignment-space '''
        key = f'{genome}_{chrom}'
        try:
            return pos + bisect_left(self._gaps[key], pos)
        except KeyError:
            logger.error(f'No liftover exists for "{key}"')


###############################################################################

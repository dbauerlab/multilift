from bisect import bisect_left, bisect_right
import logging

from Bio import SeqRecord

from multilift import __prog__


# Globals #####################################################################


__all__ = ['Lifter', 'liftover']

logger = logging.getLogger(__prog__)


# Classes #####################################################################


class Lifter():
    ''' A class that ingests Bio.SeqRecord objects from MSAs and returns
    positions lifted over to reference-space when called.

    Multilift uses Virulign for annotation-guided ORF-aware alignment, which
    makes MSAs by repetitive pairwise alignment to a specific reference. We
    specifically require a reference to be set here, therefore, as arbitrary
    liftovers would be inappropriate under this methodology (it's not true
    mutliple sequence alignment). '''

    def __init__(self) -> None:
        self._refs = {}
        self._lifts = {}

    def add_reference(self, chrom: str, record: SeqRecord) -> None:
        ''' Store a reference against which liftovers can be calculated '''
        self._refs[chrom] = record

    def add_liftover(
            self, chrom: str, record: SeqRecord, genome: str='') -> None:
        ''' Calculate a liftover table by recording the location (and hence
        number) of indels relative to the reference.
        Fills `self._lifts[key] = [[insertions], [deletions]]` for use with
        array bisection algorithms - see __call__().
        '''
        if chrom not in self._refs:
            logger.error(
                f'A "{chrom} reference must be added before liftovers '
                'can be calculated')
        key = f'{genome}_{chrom}' if genome else f'{chrom}_{chrom}'
        self._lifts[key] = [[], []]
        last_nt = -1
        for ref_nt, lift_nt in zip(self._refs[chrom], record.seq):
            if '-' == ref_nt == lift_nt:   # alignment gap for both sequences
                pass
            elif ref_nt == '-':  # insertion relative to the reference
                last_nt += 1
                self._lifts[key][0].append(last_nt)
            elif lift_nt == '-':  # deletion relative to the reference
                self._lifts[key][1].append(last_nt)
            else:  # sequence for both
                last_nt += 1

    def __call__(self, chrom: str, pos: int, genome: str='') -> int:
        ''' Liftover `chrom:pos` into reference-space for a given `genome` '''
        key = f'{genome}_{chrom}' if genome else f'{chrom}_{chrom}'
        try:
            return max(
                0,
                pos - bisect_right(self._lifts[key][0], pos)
                    + bisect_left(self._lifts[key][1], pos)
            )
        except KeyError:
            logger.error(f'No liftover has been calculated for "{key}"')

    def __repr__(self) -> str:
        return \
            'Lifter(' \
            f'n_references={len(self._refs)}, ' \
            f'n_sequences={len(self._lifts)})'


# Functions ###################################################################


# TODO: expand this out to different BED* formats
def _liftover_bed() -> None:
    pass


def _liftover_link() -> None:
    pass


def _liftover_gff() -> None:
    pass


def _liftover_wig() -> None:
    pass


def liftover(infile: str, outfile: str, lifter: Lifter) -> None:
    ''' Liftover a file '''
    pass


###############################################################################

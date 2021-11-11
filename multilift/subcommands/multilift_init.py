import argparse
import logging

from multilift import __prog__
from multilift.utils import extensions, ignored_filetypes
from multilift.utils.state import MultiliftState


# Globals #####################################################################


__all__ = ['multilift_init']

logger = logging.getLogger(__prog__)


# Functions ###################################################################


def multilift_init(args: argparse.Namespace) -> None:
    ''' Main routine for the `multilift init` subcommand '''

    with MultiliftState(args.state, overwrite=False) as S:
        S.add_genome(args.reference, True)
        for lift in args.liftover:
            S.add_genome(lift)
        for genome in S.genomes:
            for file in getattr(args, genome):
                if any(ext in ignored_filetypes for ext in extensions(file)):
                    logger.info(f'Ignoring {file}')
                    continue
                S.add_file(file, genome=genome)
        for file in args.alignment:
            if any(ext in ignored_filetypes for ext in extensions(file)):
                logger.info(f'Ignoring {file}')
                continue
            S.add_file(file, application='alignment')


###############################################################################

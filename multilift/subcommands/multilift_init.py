import argparse
import logging
from pathlib import Path

from multilift import __prog__
from multilift.utils import basename, file_hash, guess_filetype
from multilift.utils.state import MultiliftState

import yaml


# Globals #####################################################################


__all__ = ['multilift_init']

logger = logging.getLogger(__prog__)


# Functions ###################################################################


def multilift_init(args: argparse.Namespace) -> None:
    ''' Main routine for the `multilift init` subcommand '''

    with MultiliftState(args.state, overwrite=True) as S:
        S.add_genome(args.reference, True)
        for lift in args.liftovers:
            S.add_genome(lift)
        for genome in S.genomes:
            for file in getattr(args, genome):
                filetype, applications = guess_filetype(file)
                if filetype is None:
                    logger.info(f'Ignoring {file}')
                    continue
                S.add_file(genome, file, filetype=filetype)


###############################################################################

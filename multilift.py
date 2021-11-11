#!/usr/bin/env python3

import logging
import sys

from multilift import __prog__, __prog_string__, parse_args, __file__ as multilift_location


# Globals #####################################################################


logger = logging.getLogger(__prog__)


###############################################################################


if __name__ == "__main__":

    args, remainder = parse_args()
    logger.info(__prog_string__)

    if args.subcommand == 'app':
        from pathlib import Path
        import streamlit.cli as stcli
        multilift_app = \
            Path(multilift_location).resolve().parent / 'subcommands' / 'multilift_app.py'
        command = ['streamlit', 'run', str(multilift_app)] + remainder
        command += ['--', '--cache', args.cache]
        sys.argv = command
        sys.exit(stcli.main())
    elif args.subcommand == 'init':
        from multilift.subcommands.multilift_init import multilift_init
        multilift_init(args)
    elif args.subcommand == 'lift':
        from multilift.subcommands.multilift_lift import multilift_lift
        multilift_lift(args)

    logger.info('Finished!')

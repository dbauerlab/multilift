#!/usr/bin/env python3

import logging
import sys

from multilift import __prog__, __prog_string__, parse_args


# Globals #####################################################################


logger = logging.getLogger(__prog__)


###############################################################################


if __name__ == "__main__":

    args, remainder = parse_args()
    logger.info(__prog_string__)

    if args.subcommand == 'app':
        import streamlit.cli as stcli
        from multilift.subcommands.multilift_app import multilift_app_path
        sys.argv = ['streamlit', 'run', multilift_app_path]
        sys.exit(stcli.main())
    elif args.subcommand == 'init':
        from multilift.subcommands.multilift_init import multilift_init
        multilift_init(args)
    elif args.subcommand == 'lift':
        from multilift.subcommands.multilift_lift import multilift_lift
        multilift_lift(args)
    elif args.subcommand == 'fetch':
        from multilift.subcommands.multilift_fetch import multilift_fetch
        multilift_fetch(args)

    logger.info('Finished!')

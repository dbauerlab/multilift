#!/usr/bin/env python3

import logging
from pathlib import PurePath
import sys

# import streamlit.cli as stcli later (i.e. where necessary) as it's slow

from multilift import __prog__, __prog_string__, parse_args


# Globals ######################################################################


logger = logging.getLogger(__prog__)


################################################################################


if __name__ == "__main__":

    args, remainder = parse_args()
    logger.info(__prog_string__)

    if args.subcommand == 'app':
        import streamlit.cli as stcli
        sys.argv = [
            'streamlit', 'run',
            f'{PurePath(PurePath(__file__).parent, "multilift_app.py")}'
        ] + remainder
        sys.exit(stcli.main())
    elif args.subcommand == 'ini':
        pass
    elif args.subcommand == 'lift':
        pass

    print(args)

    logger.info('Finished!')

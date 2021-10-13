#!/usr/bin/env python3

import logging
from pathlib import PurePath
import sys

import streamlit.cli as stcli

from multilift import __prog_string__, parse_args


if __name__ == "__main__":

    args = parse_args()
    logging.info(__prog_string__)

    if args.subcommand == 'app':
        # TODO: allow unparsed args and forward them to streamlit
        sys.argv = [
            'streamlit', 'run',
            f'{PurePath(PurePath(__file__).parent, "multilift_app.py")}']
        sys.exit(stcli.main())
    elif args.subcommand == 'ini':
        pass
    elif args.subcommand == 'lift':
        pass

    logging.info('Finished!')

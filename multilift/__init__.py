__author__ = 'George Young, Anob Chakrabarti'
__email__ = ''

__prog__ = 'multilift'
__version__ = '0.1'
__status__ = 'Development'
__prog_string__ = f'{__prog__} v{__version__} ({__status__})'

__license__ = 'MIT license'


import argparse
from collections import deque
import logging
from multiprocessing import cpu_count


# Logging #####################################################################


class DieOnError(logging.StreamHandler):
    ''' Extends logging.StreamHandler to write error messages through the
    logging system before the program dies '''

    def __init__(self) -> None:
        super().__init__()
        self.history = deque(maxlen=6)

    def emit(self, logrecord: logging.LogRecord) -> None:
        ''' Push error message to the logging system and die if appropriate '''
        if logrecord.levelno in (logging.ERROR, logging.CRITICAL):
            logrecord.msg = f'[ERROR] {logrecord.msg}'
        elif logrecord.levelno == logging.WARNING:
            logrecord.msg = f'[WARNING] {logrecord.msg}'
        self.history.appendleft(logrecord.msg)
        if (n_identical := sum(m == logrecord.msg for m in self.history)) > 5:
            return
        super().emit(logrecord)
        if n_identical == 5:
            super().emit(logging.makeLogRecord({
                'msg': '... further instances of this message will be suppressed',
                'level': logging.INFO
            }))
        if logrecord.levelno in (logging.ERROR, logging.CRITICAL):
            raise SystemExit(-1)


logging_handler = DieOnError()
logging_handler.setFormatter(logging.Formatter(
    '%(asctime)s\t%(message)s', '%y%m%d %H:%M:%S'))
logging.getLogger(__prog__).addHandler(logging_handler)
logging.getLogger(__prog__).setLevel(logging.INFO)


# Command line parsing ########################################################


class LoggingArgumentParser(argparse.ArgumentParser):
    ''' Extends ArgumentParser to feed its error messages through the logging
    system '''
    def error(self, message: str) -> None:
        logging.getLogger(__prog__).error(message)


# top level parser
parser = LoggingArgumentParser(prog=__prog__)
# subparsers for different run modes
subparsers = parser.add_subparsers(help='commands', dest='subcommand')

# generic flags
parser.add_argument(
    '-@', '--threads', type=int, default=cpu_count(),
    help='number of processes to use (defaults max available [=%(default)s])')
parser.add_argument(
    '-v', '--version', action='version', version=__version__)
parser.add_argument(
    '-q', '--quiet', action='store_true',
    help='silence reporting of progress to stderr')

# multilift app
subparser_app = subparsers.add_parser(
    'app',
    help='launch the multilift web interface')
subparser_app.set_defaults(quiet=True)  # Running streamlit so silence stderr
subparser_app.add_argument(
    '-s', '--state', type=str,
    help='restore session state with this multilift configuration file')

# multilift init
subparser_init = subparsers.add_parser(
    'init',
    help='initialise a multilift configuration file')
subparser_init.add_argument(
    '-r', '--reference', type=str, nargs=1, required=True,
    help='name for the reference genome/sequence (no spaces)')
subparser_init.add_argument(
    '-l', '--liftovers', type=str, nargs='+', required=True,
    help='name(s) for the genome/sequence to liftover (no spaces)')
subparser_init.add_argument(
    '-a', '--alignment', type=str, nargs='+',
    help='pre-prepared alignment(s) to use for liftover calculations')

# multilift lift
subparser_lift = subparsers.add_parser(
    'lift',
    help='perform liftover')
subparser_lift.add_argument(
    '-s', '--state', type=str, required=True,
    help='perform liftover according to this multilift configuration file')

# multilift fetch
subparser_fetch = subparsers.add_parser(
    'fetch',
    help='fetch an annotated GenBank from NCBI')
subparser_fetch.add_argument(
    'accession', type=str, nargs='*',
    help='accession(s) of the nucleotide records to fetch')
subparser_fetch.add_argument(
    '--email', type=str, required=True,
    help='user email address (passed to NCBI to ensure fair usage of the API)')


def parse_args() -> (argparse.Namespace, list):
    ''' Validate the command line inputs '''
    args, unknown_args = parser.parse_known_args()
    if args.quiet:
        logging.getLogger(__prog__).setLevel(logging.ERROR)
    if args.subcommand == 'init':
        for extra_arg in args.reference + args.liftovers:
            subparser_init.add_argument(f'--{extra_arg}', type=str, nargs='+')
    if args.subcommand != 'app':
        args, unknown_args = parser.parse_args(), []
    # TODO: perform any further argument validation in here
    return args, unknown_args

################################################################################

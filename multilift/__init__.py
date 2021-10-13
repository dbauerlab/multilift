__author__ = 'George Young, Anob Chakrabarti'
__email__ = ''

__prog__ = 'multilift'
__version__ = '0.1'
__status__ = 'Development'
__prog_string__ = f'{__prog__} v{__version__} ({__status__})'

__license__ = 'MIT license'


import argparse
import logging
from multiprocessing import cpu_count


# Logging #####################################################################


class DieOnError(logging.StreamHandler):
    ''' Extends logging.StreamHandler to write error messages through the
    logging system before the program dies '''

    def emit(self, logrecord: logging.LogRecord) -> None:
        ''' Push error message to the logging system and die if appropriate '''
        if logrecord.levelno in (logging.ERROR, logging.CRITICAL):
            logrecord.msg = f'[ERROR] {logrecord.msg}'
        super().emit(logrecord)
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
    '-@', '--threads',
    type=int, default=cpu_count(),
    help='number of processes to use (defaults to number of processor cores \
    available [=%(default)s])')
parser.add_argument(
    '--version',
    action='version', version=__version__)
parser.add_argument(
    '-q', '--quiet',
    action='store_true',
    help='silence reporting of progress to stderr')

# multilift app
subparser_app = subparsers.add_parser(
    'app',
    help='Lauch the multilift web interface')
subparser_app.set_defaults(quiet=True)  # Running streamlit so silence stderr
subparser_app.add_argument(
    '--state', type=str,
    help='restore session state with this multilift configuration file')

# multilift ini
subparser_ini = subparsers.add_parser(
    'ini',
    help='Create a multilift configuration file')

# multilift lift
subparser_lift = subparsers.add_parser(
    'lift',
    help='Perform liftover')
subparser_lift.add_argument(
    '--state', type=str,
    help='perform liftover according to this multilift configuration file')


def parse_args() -> (argparse.Namespace, list):
    ''' Validate the command line inputs '''
    args, remainder = parser.parse_known_args()
    if args.quiet:
        logging.getLogger(__prog__).setLevel(logging.ERROR)
    # TODO: perform any argument validation in here
    return args, remainder

################################################################################

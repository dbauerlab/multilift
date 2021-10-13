__author__ = 'George Young, Anob Chakrabarti'
__email__ = ''

__prog__ = 'multilift'
__version__ = '0.1'
__status__ = 'Development'
__prog_string__ = f'{__prog__} v{__version__} ({__status__})'

__license__ = 'MIT license'


from argparse import ArgumentParser, Namespace
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


logging.basicConfig(
    format='%(asctime)s\t%(message)s',
    datefmt='%y%m%d %H:%M:%S',
    level=logging.INFO,
    handlers=[DieOnError()])


# Command line parsing ########################################################


class LoggingArgumentParser(ArgumentParser):
    ''' Extends ArgumentParser to feed its error messages through the logging
    system '''
    def error(self, message: str) -> None:
        logging.error(message)


# top level parser
parser = LoggingArgumentParser(prog=__prog__)
# subparsers for different run modes
subparsers = parser.add_subparsers(help='commands')

# generic / pan-subcommand flags
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


def parse_args() -> Namespace:
    ''' Validate the command line inputs '''
    args = parser.parse_args()
    if args.quiet:
        logging.getLogger().setLevel(logging.ERROR)
    # TODO: perform any argument validation in here
    return args

################################################################################

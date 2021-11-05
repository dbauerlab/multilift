import bz2
from functools import partial
import gzip
from hashlib import md5
import logging
from pathlib import Path, PurePath
from typing import Callable, Union

from multilift import __prog__


# Globals #####################################################################


__all__ = [
    'supported_filetypes', 'filetype_associations', 'ignored_filetypes',
    'basename', 'extensions', 'file_hash', 'guess_filetype', 'open_helper']

Pathish = Union[Path, PurePath, str]

logger = logging.getLogger(__prog__)

supported_filetypes = {
    'sequence': {
        'fasta': ('.fa', '.faa', '.fasta'),
        'embl': ('.embl', ),
        'genbank': ('.gb', '.genbank'),
        'snapgene': ('.dna', )
    },
    'alignment': {
        'fasta': ('.fa', '.faa', '.fasta', '.aln'),
        'clustal': ('.clustal', ),
        'maf': ('.maf', ),
        'nexus': ('.nex', '.nexus'),
        'phylip': ('.phy', '.phylip'),
        'stockholm': ('.sto', '.sth', '.stockholm')
    },
    'data': {
        'bed': ('.bed', ),
        'link': ('.link', ),
        'interact': ('inter', 'interact'),
        'bedgraph': ('.bg', '.bedgraph'),
        'wiggle': ('.wig', ),
        'narrowpeak': (),
        'broadpeak': (),
        'gappedpeak': (),
        'beddetail': ()
    },
    'annotation': {
        'bed': ('.bed', ),
        'general': ('.gtf', '.gff', '.gff2', '.gff3')
    }
}

ignored_filetypes = ('.fai', '.gzi', '.bam', '.bai', '.csi')

filetype_associations = {
    ext: filetype
    for assocs in supported_filetypes.values()
    for filetype, exts in assocs.items()
    for ext in exts}

track_types = {
    'wiggle_0': 'wiggle',
    'narrowPeak': None,
    'broadPeak': None,
    'gappedPeak': None,
    'bedDetail': None,
    'bedGraph': None}


# Functions ###################################################################


def basename(file: Pathish) -> str:
    ''' Return the basename (all extensions stripped) of a file '''
    return PurePath(file).name.partition('.')[0]


def extensions(file: Pathish, precedence_order=True) -> list[str]:
    ''' Return the extensions of a file, converted to lowercase, and optionally
    inverted to reflect precedence '''
    file = Path(file)
    exts = [ext.lower() for ext in file.suffixes]
    return exts[::-1] if precedence_order else exts


def file_hash(file: Pathish) -> str:
    ''' The MD5 checksum hex digest of the file '''
    h = md5()
    with open(file, 'rb') as F:
        for chunk in iter(lambda: F.read(4096), b''):
            h.update(chunk)
    return h.hexdigest()


def guess_filetype(file: Pathish) -> (str, [str]):
    ''' Attempt to guess the filetype from the extension, returning the
    filetype and a list of potential file applications '''
    ft = ''
    file = Path(file)
    for ext in extensions(file):
        if (ft := filetype_associations.get(ext, '')):
            break
    if not ft:
        logger.warning(f'Cannot determine filetype: {file}')
        return ft, list(supported_filetypes.keys())
    return ft, [k for k, v in supported_filetypes.items() if ft in v]


def open_helper(file: Pathish, mode='r', create: bool=True) -> Callable:
    ''' Return the appropriate function to open a file, transparently handling
    gzip or bz2 compression according to the file extension.
    mode='r'    open for reading        must be present         pointer @ 0
    mode='a'    open for read/write     create if needed        pointer @ EOF
    mode='w'    open for writing        create and truncate     pointer @ 0 '''
    file = Path(file)
    if not (dir_path := Path(PurePath(file).parent)).is_dir():
        if create:
            dir_path.mkdir(parents=True)
        else:
            logger.error(f'Cannot open file. No directory: {dir_path}')
    if '.gz' in (exts := extensions(file)):
        creator = partial(gzip.open, mode='wt')
        opener = partial(gzip.open, mode=mode + 't')
    elif '.bz2' in exts:
        creator = partial(bz2.open, mode='wt')
        opener = partial(bz2.open, mode=mode + 't')
    else:
        creator = partial(open, mode='w')
        opener = partial(open, mode=mode)
    if mode == 'r' and not file.is_file():
        if create:
            with creator(file):
                pass
        else:
            logger.error(f'Cannot open file for reading: {file}')
    return opener(file)


###############################################################################

import bz2
from functools import partial
import gzip
from hashlib import md5
import logging
from pathlib import Path, PurePath
from typing import Callable

from multilift import __prog__


# Globals #####################################################################


__all__ = [
    'supported_filetypes', 'filetype_associations',
    'basename', 'extensions', 'file_hash', 'guess_filetype', 'open_helper']

logger = logging.getLogger(__prog__)

supported_filetypes = {}
supported_filetypes['sequence'] = {
    'fasta': ('.fa', '.faa', '.fasta'),
    'embl': ('.embl', ),
    'genbank': ('.gb', '.genbank'),
    'snapgene': ('.dna', )
}
supported_filetypes['alignment'] = {
    'fasta': ('.fa', '.faa', '.fasta', '.aln'),
    'clustal': ('.clustal', ),
    'maf': ('.maf', ),
    'nexus': ('.nex', '.nexus'),
    'phylip': ('.phy', '.phylip'),
    'stockholm': ('.sto', '.sth', '.stockholm')
}
supported_filetypes['annotation'] = {
    'bed': ('.bed', ),
    'general': ('.gtf', '.gff', '.gff2', '.gff3')
}
supported_filetypes['data'] = {
    'bed': ('.bed', ),
    'link': ('.link', ),
    'interact': ('inter', 'interact'),
    'bedgraph': ('.bg', '.bedgraph'),
    'wiggle': ('.wig', ),
    'narrowpeak': (),
    'broadpeak': (),
    'gappedpeak': (),
    'beddetail': ()
}

filetype_associations = {
    ext: filetype
    for assocs in supported_filetypes.values()
    for filetype, exts in assocs.items()
    for ext in exts
}

track_types = {
    'wiggle_0': 'wiggle',
    'narrowPeak': None,
    'broadPeak': None,
    'gappedPeak': None,
    'bedDetail': None,
    'bedGraph': None
}


# Functions ###################################################################


def basename(file) -> str:
    ''' Return the basename (all extensions stripped) of a file '''
    return PurePath(file).name.partition('.')[0]


def extensions(file, precedence_order=True) -> list[str]:
    ''' Return the extensions of a file, converted to lowercase, and optionally
    inverted to reflect precedence '''
    exts = [ext.lower() for ext in Path(file).suffixes]
    return exts[::-1] if precedence_order else exts


def file_hash(file: Path) -> str:
    ''' The MD5 checksum hex digest of the file '''
    h = md5()
    with open(file, 'rb') as F:
        for chunk in iter(lambda: F.read(4096), b''):
            h.update(chunk)
    return h.hexdigest()


def guess_filetype(file: Path) -> (str, [str]):
    ''' Attempt to guess the filetype from the extension, returning the
    filetype and a list of potential file applications '''
    ft = ''
    for ext in extensions(file):
        if (ft := filetype_associations.get(ext, '')):
            break
    if not ft:
        logger.warning(f'Cannot determine filetype for "{file}"')
        return ft, list(supported_filetypes.keys())
    return ft, [k for k, v in supported_filetypes.items() if ft in v]


def open_helper(file, mode='r') -> Callable:
    ''' Return the appropriate function to open a file '''
    exts = extensions(file)
    if '.gz' in exts:
        return partial(gzip.open, mode=mode + 't')
    elif '.bz2' in exts:
        return partial(bz2.open, mode=mode + 't')
    return partial(open, mode=mode)


###############################################################################

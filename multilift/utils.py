from datetime import date
from io import BytesIO
from pathlib import Path, PurePath
import tarfile
import zipfile


def basename(file: str|Path|PurePath) -> str:
    ''' Return the basename (all extensions stripped) of a file '''
    return PurePath(file).name.partition('.')[0]


def extensions(file: str|Path|PurePath, precedence_order=True) -> list[str]:
    ''' Return the extensions of a file, converted to lowercase, and optionally
    inverted to reflect precedence '''
    exts = [ext.lower() for ext in Path(file).suffixes]
    return exts[::-1] if precedence_order else exts


def sniff_filetype(file: str|Path|PurePath) -> tuple[str, tuple[str]]:
    ''' Attempt to guess the filetype from the extension, returning the
    filetype and a list of potential file applications '''
    for ext in extensions(file):
        match ext:
            case '.fa' | '.fasta' | '.fas' | '.faa' | '.fas' | '.fna' | '.seq':
                return 'fasta', ('sequence', 'alignment' )
            case '.em' | '.emb' | '.embl':
                return 'embl', ('sequence', )
            case '.gb' | '.gen' | '.gbk' | '.genbank':
                return 'genbank', ('sequence', )
            case '.dna':
                return 'snapgene', ('sequence', )
            case '.afa':
                return 'fasta', ('alignment', )
            case '.clustal' | '.aln':
                return 'clustal', ('alignment', )
            case '.maf':
                return 'maf', ('alignment', )
            case '.nex' | '.nxs' | '.nexus':
                return 'nexus', ('alignment', )
            case '.ph' | '.phy' | '.phylip':
                return 'phylip', ('alignment', )
            case '.sto' | '.sth' | '.stockholm':
                return 'stockholm', ('alignment', )
            case '.bed':
                return 'bed', ('data', )
            case '.bedpe':
                return 'link', ('data', )
            case 'inter' | 'interact':
                return 'interact', ('data', )
            case '.bg' | '.bedgraph':
                return 'bedgraph', ('data', )
            case '.wig':
                return 'wiggle', ('data', )
            case '.gtf' | '.gff' | '.gff2' | '.gff3':
                return 'gtf', ('data', )
            case '.dot' | '.db' | '.bracket' | '.dbn':
                return 'dotbracket', ('data', )
            case '.bp':
                return 'basepair', ('data', )
            case '.dp':
                return 'dotplot', ('data', )
            case '.yaml':
                return 'annotation', ('annotation', )
    return '', ()


def create_igv_session(session_id: str, igv_resources: list[str]) -> str:
    ''' '''
    xml_string = '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
    xml_string += f'<Session genome="genome/multilift_{session_id}.fa" hasGeneTrack="false" hasSequenceTrack="false" locus="All" version="8">\n'
    xml_string += '    <Resources>\n'
    for r in igv_resources:
        xml_string += f'        <Resource path="{r}"/>\n'
    xml_string += '    </Resources>\n'
    xml_string += '</Session>'
    return xml_string


def add_to_archive(
    arc: BytesIO, data: BytesIO, name: str, arc_type: str='.tar.gz') -> None:
    ''' '''
    if arc_type == '.tar.gz':
        info = tarfile.TarInfo(name)
        info.size = len(data.getbuffer())
        arc.addfile(info, data)
    elif arc_type == '.zip':
        info = zipfile.ZipInfo(name, date.timetuple(date.today())[:6])
        info.file_size = len(data.getbuffer())
        arc.writestr(info, data.getbuffer())

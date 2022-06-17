from bisect import bisect_left, bisect_right
from collections import defaultdict as dd
from difflib import get_close_matches
from io import StringIO
import re



# Globals #####################################################################


header_re = re.compile(r'(\S+)=(\w+?|"[\S ]+?")')


# Classes #####################################################################


class Lifter():
    ''' A class that ingests Bio.SeqRecord objects from MSAs and returns
    positions lifted over to reference-space when called '''

    def __init__(self) -> None:
        self.liftovers = dd(dict)

    def add_alignment(self, alignment: dict, ref_genome: str):
        '''  '''
        ref = alignment.pop(ref_genome)
        _, _, ref_seqid = ref.description.partition(' ')
        for s in alignment.values():
            genome = s.id
            _, _, seqid = s.description.partition(' ')
            self.liftovers[genome][seqid] = (ref_seqid, [[], []])
            last_nt = -1
            for ref_nt, lift_nt in zip(ref.seq, s.seq):
                if '-' == ref_nt == lift_nt:  # alignment gap for both
                    pass
                elif ref_nt == '-':  # insertion relative to the reference
                    last_nt += 1
                    self.liftovers[genome][seqid][1][0].append(last_nt)
                elif lift_nt == '-':  # deletion relative to the reference
                    self.liftovers[genome][seqid][1][1].append(last_nt)
                else:  # sequence for both
                    last_nt += 1

    def __call__(self, genome: str, seqid: str, pos: int) -> (str, int):
        ''' Liftover `seqid:pos` into reference-space for a given `genome` '''
        if not seqid in self.liftovers[genome]:
            matches = get_close_matches(seqid, self.liftovers[genome].keys(), 1)
            if matches and \
                    (matches[0].startswith(seqid) or seqid.startswith(matches[0])):
                seqid = matches[0]
            else:
                raise ValueError(
                    f'No liftover has been calculated for {genome} {seqid}')
        return (
            self.liftovers[genome][seqid][0],
            max(
                0,
                pos - bisect_right(self.liftovers[genome][seqid][1][0], pos)
                    + bisect_left(self.liftovers[genome][seqid][1][1], pos)))


    def __repr__(self) -> str:
        return 'Lifter()'


# Functions ###################################################################


def parse_header(s: str) -> dict:
    ''' Parse space-separated key=value pairs from a `track` line '''
    return {k: v.strip('"') for k, v in header_re.findall(s)}


def format_header(meta: dict, prefix: str='') -> str:
    ''' Format a dictionary to space-separated key=value pairs '''
    meta = {k: f'"{v}"' if ' ' in v else v for k, v in meta.items()}
    header = f'{prefix} ' if prefix else ''
    header += ' '.join(f'{k}={v}' for k, v in meta.items())
    return header


def _liftover_bed(
        infile: StringIO, lifter: Lifter, genome: str) -> tuple[str, StringIO]:
    '''
    PRIVATE. Liftover for BED* format data

    Format spec: https://genome.ucsc.edu/FAQ/FAQformat.html#format1

    !!! BED-based files are 0-based half-open ... [0, x)
    '''
    outfile = StringIO()
    for line in infile:
        if not (line := line.strip('\r\n')) or \
                any(line.startswith(s) for s in ('#', 'browser')):
            continue
        elif line.startswith('track'):
            print(line, file=outfile)
        else:
            line = line.split('\t')
            line[0], line[1] = lifter(genome, seqid := line[0], int(line[1]))
            _, line[2] = lifter(genome, seqid, int(line[2]))
            if len(line) > 6:
                _, line[6] = lifter(genome, seqid, int(line[6]))
                _, line[7] = lifter(genome, seqid, int(line[7]))
            # TODO: Also need to do line[11] and line[12] adjustment!
            print('\t'.join(str(x) for x in line), file=outfile)
    return 'bed', outfile


def _liftover_bedgraph(
        infile: StringIO, lifter: Lifter, genome: str) -> tuple[str, StringIO]:
    '''
    PRIVATE. Liftover for BED* format data

    Format spec: https://genome.ucsc.edu/goldenpath/help/bedgraph.html

    !!! BED-based files are 0-based half-open ... [0, x)
    '''
    outfile = StringIO()
    for line in infile:
        if not (line := line.strip('\r\n')) or \
                any(line.startswith(s) for s in ('#', 'browser')):
            continue
        elif line.startswith('track'):
            print(line, file=outfile)
        else:
            line = line.split('\t')
            line[0], line[1] = lifter(genome, seqid := line[0], int(line[1]))
            _, line[2] = lifter(genome, seqid, int(line[2]))
            print('\t'.join(str(x) for x in line), file=outfile)
    return 'bedgraph', outfile


def _liftover_interact(
        infile: StringIO, lifter: Lifter, genome: str) -> tuple[str, StringIO]:
    '''
    PRIVATE. Liftover for interact (BED5+13) format data

    Format spec: https://genome.ucsc.edu/goldenPath/help/interact.html

    !!! BED-based files are 0-based half-open ... [0, x)
    '''
    outfile = StringIO()
    for line in infile:
        if not (line := line.strip('\r\n')) or \
                any(line.startswith(s) for s in ('#', 'browser')):
            continue
        elif line.startswith('track'):
            print(line, file=outfile)
        else:
            line = line.split('\t')
            line[0], line[1] = lifter(genome, seqid := line[0], int(line[1]))
            _, line[2] = lifter(genome, seqid, int(line[2]))
            line[8], line[9] = lifter(genome, seqid := line[8], int(line[9]))
            _, line[10] = lifter(genome, seqid, int(line[10]))
            line[13], line[14] = lifter(genome, seqid := line[13], int(line[14]))
            _, line[15] = lifter(genome, seqid, int(line[15]))
            print('\t'.join(str(x) for x in line), file=outfile)
    return 'interact', outfile


def _liftover_link(
        infile: StringIO, lifter: Lifter, genome: str) -> tuple[str, StringIO]:
    '''
    PRIVATE. Liftover for linked (BED3+BED3+n) files

    !!! BED-based files are 0-based half-open ... [0, x)
    '''
    outfile = StringIO()
    for line in infile:
        if not (line := line.strip('\r\n')) or \
                any(line.startswith(s) for s in ('#', 'browser')):
            continue
        elif line.startswith('track'):
            print(line, file=outfile)
        else:
            line = line.split('\t')
            line[0], line[1] = lifter(genome, seqid := line[0], int(line[1]))
            _, line[2] = lifter(genome, seqid, int(line[2]))
            line[3], line[4] = lifter(genome, seqid := line[0], int(line[4]))
            _, line[5] = lifter(genome, seqid, int(line[5]))
            print('\t'.join(str(x) for x in line), file=outfile)
    return 'link', outfile


def _liftover_gxf(
        infile: StringIO, lifter: Lifter, genome: str) -> tuple[str, StringIO]:
    '''
    PRIVATE. Liftover for GFF and GTF format data

    Format spec: http://www.ensembl.org/info/website/upload/gff.html

    !!! GFF and GTF are 1-based fully-closed ... [1, x]
    '''
    outfile = StringIO()
    for line in infile:
        if not (line := line.strip('\r\n')) or \
                any(line.startswith(s) for s in ('#', 'browser')):
            continue
        elif line.startswith('track'):
            print(line, file=outfile)
        else:
            line = line.split('\t')
            line[0], line[3] = lifter(genome, seqid := line[0], int(line[3])-1)
            line[3] += 1
            _, line[4] = lifter(genome, seqid, int(line[4])-1)
            line[4] += 1
            print('\t'.join(str(x) for x in line), file=outfile)
    return 'gtf', outfile


def _liftover_wiggle(
        infile: StringIO, lifter: Lifter, genome: str) -> tuple[str, StringIO]:
    '''
    PRIVATE. Liftover for wiggle format data.

    Format spec: https://genome.ucsc.edu/goldenPath/help/wiggle.html

    There is no means of guaranteeing the size of data spans when converting
    wiggle format data, so we convert to bedgraph.

    !!! .wig format is 1-based fully-closed ... [1, x]
    !!! .bedgraph format is 0-based half-open ... [0, x)
    '''
    outfile = StringIO()
    for line in infile:
        if not (line := line.strip('\r\n')) or \
                any(line.startswith(s) for s in ('#', 'browser')):
            continue
        elif line.startswith('track'):
            track_meta = parse_header()
            track_meta['type'] = 'bedGraph'
            print(format_header(track_meta, 'track'), file=outfile)
        elif any(line.startswith(s) for s in ('variableStep', 'fixedStep')):
            if line.startswith('variableStep'):
                wig_meta = {'wig_type': 'variableStep'}
                wig_meta.update(parse_header(line))
            elif line.startswith('fixedStep'):
                wig_meta = {'wig_type': 'fixedStep'}
                wig_meta.update(parse_header(line))
                wig_meta['start'] = int(wig_meta['start'])
                wig_meta['step'] = int(wig_meta['step'])
                current_step = -1
            wig_meta['span'] = int(wig_meta.get('span', 1))
        else:  # a data line
            if wig_meta['wig_type'] == 'variableStep':
                start, value = line.split()
                start = int(start)
            elif wig_meta['wig_type'] == 'fixedStep':
                current_step += 1
                start = wig_meta['start'] + (current_step * wig_meta['step'])
                value = line
            # -1 converts to 0-based
            chrom, start = lifter(genome, wig_meta['chrom'], start - 1)
            # +1 converts to half-open
            _, stop = lifter(
                genome, wig_meta['chrom'], start + wig_meta['span'] + 1)
            print(
                f'{chrom}\t{str(start)}\t{str(stop)}\t{value}',
                file=outfile)
    return 'bedgraph', outfile


def liftover(infile: StringIO, ftype: str, lifter: Lifter, genome: str):
    match ftype:
        case 'bed':
            return _liftover_bed(infile, lifter, genome)
        case 'bedgraph':
            return _liftover_bedgraph(infile, lifter, genome)
        case 'link':
            return _liftover_link(infile, lifter, genome)
        case 'interact':
            return _liftover_interact(infile, lifter, genome)
        case 'wiggle':
            return _liftover_wiggle(infile, lifter, genome)
        case 'gtf':
            return _liftover_gxf(infile, lifter, genome)

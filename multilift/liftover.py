from bisect import bisect_left, bisect_right
from collections import defaultdict as dd, deque
from difflib import get_close_matches
from io import StringIO
import re



# Globals #####################################################################


header_re = re.compile(r'(\S+)=("[\S ]+"|\w+)')


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
                    f'No liftover has been calculated for sequence ID {seqid} for genome {genome}')
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
    return '', outfile


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
    return '', outfile


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
    return '', outfile


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
    return '', outfile


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
    return '', outfile


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
            track_meta = parse_header(line)
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
                f'{chrom}\t{start}\t{stop}\t{value}',
                file=outfile)
    return '.bedgraph', outfile


def _liftover_dotbracket(
        infile: StringIO, lifter: Lifter, genome: str) -> tuple[str, StringIO]:
    '''
    PRIVATE. Liftover for dot bracket format data.

    Format spec: https://software.broadinstitute.org/software/igv/RNAsecStructure

    As brackets have to be exactly paired, there's no straightforward way to
    liftover dot bracket formats directly. We instead parse the data and output
    as .bp (base pair) format.

    !!! .bp format is 1-based fully-closed ... [1, x]
    '''
    outfile = StringIO()
    seqid = infile.readline().strip('\r\n')[1:].partition(' ')[0]
    infile.readline()  # ignore the sequence line
    deques = {'(': deque(), '[': deque(), '{': deque(), '<': deque()}
    brackets = {')': '(', ']': '[', '}': '{', '>': '<'}
    pairs = []
    for i, c in enumerate(infile.readline().strip('\r\n').partition(' ')[0]):
        try:
            deques[c].append(i)
        except KeyError:
            try:
                pairs.append(
                    (c,
                    lifter(genome, seqid, deques[brackets[c]].pop())[1],
                    lifter(genome, seqid, i)[1]))
            except KeyError:
                # unpaired base
                pass
    print('color:\t0\t0\t0', file=outfile)
    ref_seqid = lifter(genome, seqid, 0)[0]
    pairs.sort(key=lambda x: (x[1], -x[2]))
    current = [(p := pairs[0])[0], p[1], p[1], p[2], p[2]]
    for p in pairs[1:]:
        if (p[0] == current[0]) \
                and (p[1] == current[2] + 1) \
                and (p[2] == current[3] - 1):
            current[2] += 1
            current[3] -= 1
        else:
            print(
                f'{ref_seqid}\t'
                f'{current[1]+1}\t{current[2]+1}\t'
                f'{current[3]+1}\t{current[4]+1}\t'
                '0',
                file=outfile)
            current = [p[0], p[1], p[1], p[2], p[2]]
    print(
        f'{ref_seqid}\t'
        f'{current[1]+1}\t{current[2]+1}\t'
        f'{current[3]+1}\t{current[4]+1}\t'
        '0',
        file=outfile)
    return '.bp', outfile


def _liftover_basepair(
        infile: StringIO, lifter: Lifter, genome: str) -> tuple[str, StringIO]:
    '''
    PRIVATE. Liftover for RNA base pairing format data.

    Format spec: https://software.broadinstitute.org/software/igv/RNAsecStructure

    !!! .bp format is 1-based fully-closed ... [1, x]
    '''
    outfile = StringIO()
    for line in infile:
        if not (line := line.strip('\r\n')) or \
                any(line.startswith(s) for s in ('#', )):
            continue
        elif line.startswith('color'):
            print(line, file=outfile)
        else:
            line = line.split('\t')
            line[1:5] = [int(x)-1 for x in line[1:5]]
            pairs = [
                (lifter(genome, seqid, p[0])[1], lifter(genome, seqid, p[1])[1])
                for p in zip(
                    range(line[1], line[2]+1), range(line[4], line[3]-1, -1))]
            ref_seqid = lifter(genome, line[0], 0)[0]
            current = [(p := pairs[0])[0], p[0], p[1], p[1]]
            for p in pairs[1:]:
                if (p[0] == current[1] + 1) and (p[1] == current[2] - 1):
                    current[1] += 1
                    current[2] -= 1
                else:
                    print(
                        f'{ref_seqid}\t'
                        f'{current[0]+1}\t{current[1]+1}\t'
                        f'{current[2]+1}\t{current[3]+1}\t'
                        f'{line[5]}',
                        file=outfile)
                    current = [p[0], p[0], p[1], p[1]]
            print(
                f'{ref_seqid}\t'
                f'{current[0]+1}\t{current[1]+1}\t'
                f'{current[2]+1}\t{current[3]+1}\t'
                f'{line[5]}',
                file=outfile)
    return '', outfile


def _liftover_dotplot(
        infile: StringIO, lifter: Lifter, genome: str) -> tuple[str, StringIO]:
    '''
    PRIVATE. Liftover for RNA base pairing probability data.

    Format spec: https://software.broadinstitute.org/software/igv/RNAsecStructure

    !!! .bp format is 1-based fully-closed ... [1, x]
    '''
    if len(lifter.liftovers[genome]) > 1:
        raise ValueError(
            'Cannot lift over dotplot (.dp) data where more than one sequence \
            has been provided for the genome. Data should be converted to \
            base pair (.bp) format.')
    seqid = list(lifter.liftovers[genome].keys())[0]
    infile.readline()  # length of sequence
    infile.readline()  # column headers
    prob_pairs = [[], [], []]  # 3 for the different probability ranges
    for line in infile:
        if not (line := line.strip('\r\n')) or \
                any(line.startswith(s) for s in ('#', )):
            continue
        left, right, prob = line.split('\t')
        left = int(left) - 1
        right = int(right) - 1
        prob = 10 ** -float(prob)
        if 0.1 <= prob < 0.3:
            prob_pairs[0].append(
                (lifter(genome, seqid, left)[1],
                lifter(genome, seqid, right)[1]))
        elif 0.3 <= prob < 0.8:
            prob_pairs[1].append(
                (lifter(genome, seqid, left)[1],
                lifter(genome, seqid, right)[1]))
        elif 0.8 <= prob:
            prob_pairs[2].append(
                (lifter(genome, seqid, left)[1],
                lifter(genome, seqid, right)[1]))
    prob_pairs = \
        [sorted(pairs, key=lambda x: (x[0], -x[1])) for pairs in prob_pairs]
    ref_seqid = lifter(genome, seqid, 0)[0]
    outfile = StringIO()
    print(
        'color:\t255\t218\t125\tPairing probability 10-30%\n'
        'color:\t113\t195\t209\tPairing probability 30-80%\n'
        'color:\t51\t114\t38\tPairing probability >80%',
        file=outfile)
    for prob_level, pairs in enumerate(prob_pairs):
        current = [(p := pairs[0])[0], p[0], p[1], p[1]]
        for p in pairs[1:]:
            if (p[0] == current[1] + 1) and (p[1] == current[2] - 1):
                current[1] += 1
                current[2] -= 1
            else:
                print(
                    f'{ref_seqid}\t'
                    f'{current[0]+1}\t{current[1]+1}\t'
                    f'{current[2]+1}\t{current[3]+1}\t'
                    f'{prob_level}',
                    file=outfile)
                current = [p[0], p[0], p[1], p[1]]
        print(
            f'{ref_seqid}\t'
            f'{current[0]+1}\t{current[1]+1}\t'
            f'{current[2]+1}\t{current[3]+1}\t'
            f'{prob_level}',
            file=outfile)
    return '.bp', outfile


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
        case 'dotbracket':
            return _liftover_dotbracket(infile, lifter, genome)
        case 'basepair':
            return _liftover_basepair(infile, lifter, genome)
        case 'dotplot':
            return _liftover_dotplot(infile, lifter, genome)

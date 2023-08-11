import tempfile
import subprocess
import gzip
import os
from typing import List, Iterable

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import motifs
from search_n2 import wrapper
from html_output import to_html
import enum
import itertools

TFMP = False
try:
    from pytfmpval import tfmp
    TFMP = True
except ImportError:
    print("tfmpval not found and can't be used")


# coolors pallete
class Pallete(object):
    def __init__(self):
        self.p = {
            'vivid_sky_blue': {'DEFAULT': '#0ad2ff', 100: '#002b35', 200: '#00576a', 300: '#00829f', 400: '#00add4', 500: '#0ad2ff', 600: '#3bdbff', 700: '#6ce4ff', 800: '#9dedff', 900: '#cef6ff'},
            'neon_blue': {'DEFAULT': '#2962ff', 100: '#00103b', 200: '#002076', 300: '#002fb1', 400: '#003fed', 500: '#2962ff', 600: '#5481ff', 700: '#7ea1ff', 800: '#a9c0ff', 900: '#d4e0ff'},
            'electric_violet': {'DEFAULT': '#9500ff', 100: '#1e0033', 200: '#3b0066', 300: '#590099', 400: '#7700cc', 500: '#9500ff', 600: '#aa33ff', 700: '#bf66ff', 800: '#d499ff', 900: '#eaccff'},
            'folly': {'DEFAULT': '#ff0059', 100: '#330012', 200: '#660024', 300: '#990036', 400: '#cc0047', 500: '#ff0059', 600: '#ff337a', 700: '#ff669c', 800: '#ff99bd', 900: '#ffccde'},
            'dark_orange_(web)': {'DEFAULT': '#ff8c00', 100: '#331c00', 200: '#663800', 300: '#995400', 400: '#cc7000', 500: '#ff8c00', 600: '#ffa333', 700: '#ffba66', 800: '#ffd199', 900: '#ffe8cc'},
            'lime': {'DEFAULT': '#b4e600', 100: '#242e00', 200: '#485c00', 300: '#6c8a00', 400: '#90b800', 500: '#b4e600', 600: '#ceff1f', 700: '#dbff57', 800: '#e7ff8f', 900: '#f3ffc7'},
            'aquamarine': {'DEFAULT': '#0fffdb', 100: '#00362e', 200: '#006c5c', 300: '#00a28a', 400: '#00d8b8', 500: '#0fffdb', 600: '#3fffe2', 700: '#6fffe9', 800: '#9ffff1', 900: '#cffff8'}
            }
        self.i = 0
        self.colors = list(self.p.keys())

    def get_color(self, shade='DEFAULT'):
        _i = self.i
        self.next_color()
        return self.p[self.colors[_i]][shade]

    def next_color(self):
        if self.i == len(self.colors) - 1:
            self.i = 0
        else:
            self.i += 1


class FilterBy(enum.Enum):
    TFMP = enum.auto()
    SCORE = enum.auto()


class Filter(object):
    def __init__(self, filter_kind='score', score_thresh=0, pval_thresh=0.05):
        if filter_kind == 'score':
            self.filterby = FilterBy.SCORE
        elif filter_kind == 'pval':
            self.filterby = FilterBy.TFMP
        else:
            raise ValueError('Invalid ft kind. Only "score" and "pval" allowed.')
        self.filterThreshScore = score_thresh
        self.filterThreshPval = pval_thresh


class Diagram(object):
    """Diagram class
    inspired by MAST diagrams
    If we have motifs named a, b, c and we expect that they appear in order of b-a-c
    with 6-8 spacer in between the first two and 7 between the second two
    we would write [b]6-8[a]7[c]
    we would write [+b]6-8[-a]7[+c]
    Motif names cannot contain [] characters
    Motif names cannot start with '+' or '-' characters
    """
    def __init__(self, diagram: str):
        self.d = diagram

        br_pairs = self._find_brackets(diagram)
        self.models, self.strands = self._find_model_names(diagram, br_pairs)
        self.dists = self._find_distances(diagram, br_pairs)
        self.matches = []

        # joined motif matrix for p-value computation from combined scores
        self.tfmp = None

    def to_df(self):
        df = pd.DataFrame(r.pretty_dict() for r in self.matches)
        return df

    @staticmethod
    def _find_model_names(diagram: str, bracket_pairs: List[tuple]) -> (List[str], List[str]):
        model_names = []
        model_strands = []
        for l, r in bracket_pairs:
            m = diagram[l + 1:r]
            s = '+'
            if m.startswith('-'):
                s = '-'
                m = m[1:]

            model_names.append(m)
            model_strands.append(s)
        return model_names, model_strands

    @staticmethod
    def _find_distances(diagram: str, bracket_pairs: List[tuple]) -> List[tuple]:
        dists = []
        for i in range(len(bracket_pairs)-1):
            s = diagram[bracket_pairs[i][1]+1:bracket_pairs[i+1][0]]
            spacer_range = s.split('-')
            if len(spacer_range) == 1:
                # we have only one distance
                d = int(spacer_range[0])+1  # convert (correct) from spacer length to distance
                dists.append((d, d))
            elif len(spacer_range) == 2:
                # we have a range
                dists.append([int(d)+1 for d in spacer_range])  # convert (correct) from spacer length to distance
            else:
                raise ValueError("Invalid spacer range specification: {}".format(s))
        return dists

    @staticmethod
    def _find_brackets(diagram: str) -> List[tuple]:
        output = []
        stack = []
        for i, t in enumerate(diagram):
            if t == '[':
                stack.append(i)
            elif t == ']':
                if len(stack) == 0:
                    raise ValueError('Unbalanced diagram: unexpected closing bracket at {}'.format(i))
                output.append((stack.pop(), i))

            if len(stack) > 1:
                raise ValueError('Unbalanced diagram: unexpected opening bracket at {}'.format(i))
        if len(stack) > 0:
            raise ValueError('Unbalanced diagram: too few closing brackets.')
        return output

    def site_combinations(self, df):
        # access the sequence name, so we don't have to take care of it latter
        sn = df.sequence_name.iloc[0]
        df['se'] = list(zip(df.start, df.stop))
        df.set_index(['motif_id', 'strand', 'se'], drop=False, inplace=True)

        def _process_strand(orientation, models, model_strands, dists):
            compatible_tupes = wrapper([
                df.loc[(df.motif_id == m) & (df.strand == ms), 'se'].tolist() for m, ms in zip(models, model_strands)
            ], dists)

            # maybe just provide list of per motif-strand tuples
            for one_csite in compatible_tupes:
                positions = []
                matched_sequences = []
                scores = []
                pvals = []
                for mi in zip(models, model_strands, one_csite):
                    r = df.xs(mi)
                    positions.append(mi[2])
                    matched_sequences.append(r.matched_sequence)
                    scores.append(r['score'])
                    pvals.append(r['p-value'])

                yield CombinedMatch(
                    self,
                    sequence_name=sn,
                    orientation=orientation,
                    scores=scores,
                    matched_sequences=matched_sequences,
                    sites=positions,
                    pvals=pvals,
                    strands=model_strands
                )

        # make process strand lookup combinations in + and -
        # but take into account that the diagram may contain both orientations (the 1 and -1)
        # this does not interfere with the specification but rather provides more flexibility
        for c in _process_strand('fw', self.models, self.strands, self.dists):
            yield c
        for c in _process_strand('rc', self.models[::-1], [reverse_strand(s) for s in self.strands[::-1]], self.dists[::-1]):
            yield c

    def to_gff(self, completeGFF=False):
        # todo: add code to make complete gff record (headers and footers)
        return '\n'.join(itertools.chain.from_iterable(m.to_gff(id=str(i+1)) for i, m in enumerate(self.matches)))

    def to_interact(self):
        return '\n'.join(itertools.chain.from_iterable(m.to_interact() for m in self.matches))

    def to_vis_tracks(self):
        pal = Pallete()
        for i, m in enumerate(self.matches):
            c = pal.get_color()
            gff = m.to_gff(site_id=str(i+1), color=c)
            interact = m.to_interact(site_id=str(i+1), color=c)
            yield '\n'.join(gff), '\n'.join(interact)


def reverse_strand(strand):
    if strand == '+':
        return '-'
    elif strand == '-':
        return '+'
    else:
        raise ValueError(f'Strand may only be +/- but {strand} given.')


class CombinedMatch(object):
    def __init__(self, diagram: Diagram, sequence_name: str, orientation: str, scores: list[float], matched_sequences: list[str], sites: list[tuple[int]], pvals: list[str], strands: list[str]):
        self.sequence_name = sequence_name
        self.diagram = diagram.d
        self._diagram_names = diagram.models

        self.orientation = orientation  # orientation switches order and strand of the matches
        self.scores = scores
        self.matched_sequences = matched_sequences
        self.sites = sites
        self.pvals = pvals
        self.strands = strands

        self.pval = float('nan')   # prealocate for tfmp_pval
        self.score = sum(scores)

    def to_dict(self) -> dict:
        return dict(vars(self))

    def to_dict_reverse(self) -> dict:
        d = self.to_dict()
        d['scores'] = d['scores'][::-1]
        d['matched_sequences'] = d['matched_sequences'][::-1]
        d['sites'] = d['sites'][::-1]
        d['pvals'] = d['pvals'][::-1]
        d['strands'] = d['strands'][::-1]
        # do not reverse the _diagram names, it is needed for correct visualization of sites in igv
        # d['_diagram_names'] = d['_diagram_names'][::-1]
        return d

    def pretty_dict(self, sep='==') -> dict:
        """Make nicer output"""

        if self.orientation == 'fw':
            d = self.to_dict()
            d['dists'] = sep.join(str(self.sites[i + 1][0] - self.sites[i][1]-1) for i in range(len(self.sites) - 1))
        elif self.orientation == 'rc':
            d = self.to_dict_reverse()
            # compute dists from the reported site organization, then reverse
            # need to check dists.
            d['dists'] = sep.join([str(self.sites[i + 1][0] - self.sites[i][1]-1) for i in range(len(self.sites) - 1)][::-1])
        else:
            raise ValueError(f'Invalid diagram orientation. Only "fw" or "rc" accepted, but {self.orientation} recieved.')

        d['matched_sequences'] = sep.join(d['matched_sequences'])
        d['sites'] = sep.join(f'({i[0]}-{i[1]})' for i in d['sites'])
        d['pvals'] = sep.join(d['pvals'])

        # do not output _diagram_names to
        del d['_diagram_names']
        return d

    def to_gff(self, site_id='', color='blue') -> list[str]:
        rows = []
        fstr = '{sn}\tfimo\t{tstr}\t{start}\t{end}\t{sc}\t{st}\t.\tName={name};ID={site_id};pval={pval};sequence={ms};color={color};'

        if self.orientation == 'fw':
            d = self.to_dict()
        else:
            d = self.to_dict_reverse()

        for i in range(len(d['sites'])):
            dn = d['strands'][i]+d['_diagram_names'][i]
            gffr = fstr.format(
                sn=d['sequence_name'],
                tstr='nucleotide_motif',
                start=d['sites'][i][0],
                end=d['sites'][i][1],
                sc=d['scores'][i],
                st=d['strands'][i],
                name=dn,
                site_id=site_id+f'_{dn}_{i}',
                pval=d['pvals'][i],
                ms=d['matched_sequences'][i],
                color=color,
            )
            rows.append(gffr)
        return rows

    def to_interact(self, site_id='.', color='blue') -> list[str]:
        # interact file is in pairs
        # we thus need n-1 rows
        # make span for whole interact
        # start - 1 correction to render block correctly in igv.js
        lines = []

        if self.orientation == 'fw':
            d = self.to_dict()
        else:
            d = self.to_dict_reverse()

        sn = d['sequence_name']
        sites = d['sites']
        score = d['score']
        dn = d['_diagram_names']
        strands = d['strands']
        ma = max(list(itertools.chain.from_iterable(sites)))
        mi = min(list(itertools.chain.from_iterable(sites)))
        span = f'{sn}\t{mi}\t{ma}\t{site_id}\t{round(score)}\t{score}\t.\t{color}'
        for i in range(len(sites)-1):
            ss = sites[i][0]-1
            se = sites[i][1]
            ts = sites[i+1][0]-1
            te = sites[i+1][1]
            r = f'{span}\t{sn}\t{ss}\t{se}\t{dn[i]}\t{strands[i]}\t{sn}\t{ts}\t{te}\t{dn[i+1]}\t{strands[i+1]}'
            lines.append(r)
        return lines


def fimo_wrapper(motifs: str, fasta: Iterable[SeqRecord], options=None):
    for seq in fasta:
        yield run_fimo(motifs, seq, options)


def run_fimo(motifs: str, fasta: SeqRecord, options=None) -> str:
    """Run FIMO sequence by sequence
    Although this may seem inefficient, it produces sequence by sequence output and eliminates the need for sorting the
    output later on or reading the huge output all at once.
    It also allows for multiple threads to be used. (Not implemented yet, see fimo_wrapper)

    FIMO is assumed to run in text mode. Enforce it here.
    """
    options = [] if options is None else options

    # Enforce the text mode
    if '--text' not in options:
        options += ['--text', ]

    fd, fp = tempfile.mkstemp()
    with os.fdopen(fd, 'w') as f:
        SeqIO.write([fasta], f, 'fasta')

    r = subprocess.run(['fimo'] + options + [motifs, fp], capture_output=True, text=True)
    os.remove(fp)
    return r.stdout


def analyze_output(fo: str, diagrams: List[Diagram], ft: Filter):
    df = pd.DataFrame().from_records(i.split('\t') for i in fo.split('\n') if i)
    df.columns = df.iloc[0]
    df = df[1:]
    df.loc[:, ['start', 'stop']] = df.loc[:, ['start', 'stop']].astype(int)
    df.score = df.score.astype(float)

    seqnames = df.sequence_name.unique()
    if len(seqnames) != 1:
        raise ValueError(f'Only one sequence per search allowed; problem with: {seqnames}')

    for d in diagrams:
        for comb in d.site_combinations(df):
            if ft.filterby is FilterBy.SCORE:
                if comb.score >= ft.filterThreshScore:
                    d.matches.append(comb)

            elif ft.filterby is FilterBy.TFMP:
                pval = tfmp.score2pval(d.tfmp, comb.score)

                # ============= p-value float error handling =========================================================
                # the user is allowed to specify the pval_thresh of 1
                #  since the tfmp.score2pval returns floats slightly larger than 1 for very low scores
                #  (probably as a result of float errors), we set p-value as 1 directly to avoid downstream issues.
                if pval > 1:
                    pval = float(1)         # needed due to floating point errors
                # ====================================================================================================

                if pval <= ft.filterThreshPval:
                    comb.pval = pval
                    d.matches.append(comb)
            else:
                raise NotImplementedError('Filter option not implemented. Please contact developers.')

    return diagrams


def join_PFMs_to_dict(pfms, alphabet='ACGT'):
    """Join PFMs dict representation in order in which they are given.
    They must share the same alphabet (specified).
    """
    op = {a: [] for a in alphabet}
    for p in pfms:
        for a in alphabet:
            op[a] += p.counts[a]

            if a not in p.alphabet:
                raise ValueError(f'Alphabet mismatch. Expected {alphabet} but got {p.alphabet}.')
    return op


def join_motif_background_frequencies(pfms, alphabet='ACGT'):
    bg = {}
    for a in alphabet:
        bgs = [p.background[a] for p in pfms]
        if len(set(bgs)) == 1:
            bg[a] = bgs[0]
        else:
            raise ValueError(f'Background frequency mismatch. Expected same freq. for all motifs but got {a}: {bgs}.')
    return bg


def main(
    motif_file='',
    fasta='',
    diag_match_output_file='',
    fimo_options=('--thresh', '0.05'),
    diagrams=(),
    vis_track_prefix='',
    filter_score=0,
):
    """

    :param motif_file: File with motifs.
    :param fasta: Fasta file to search.
    :param diag_match_output_file: Output file .xlsx|.tsv|.csv.
    :param fimo_options: Fimo options. Must be compatible with "--text" option.
    :param diagrams: list of diagram strings
    :param vis_track_prefix: path prefix for .gff and .interact files
    :return:
    """
    fimo_options = list(fimo_options)
    filter_options = Filter(score_thresh=filter_score)

    # todo - We don't need to enforce alphabet if we use the score cuttoff
    # todo - Infer alphabet from pfms when score cutoff so it can be used in further functions
    alphabet = 'ACGT'       # this is limitation of the tfmpval package
    D = [Diagram(d) for d in diagrams]

    # implement for combined score first
    # oh, this must be done for each diagram separately
    # todo: implement merging the motifs specified in the diagram
    M = dict()
    with open(motif_file) as f:
        for m in motifs.parse(f, 'MINIMAL'):
            if filter_options.filterby is FilterBy.TFMP and set(m.alphabet) != set(alphabet):
                raise ValueError(f'Invalid alphabet {m.alphabet} for {m.name}. Only {alphabet} is supported.')
            M[m.name] = m

    # otherwise we need user input to determine how they want to handle the situation
    for d in D:
        # note that variable spaces are ignored in the scoring procedure

        # ======== Prepare the joined motif =========
        diag_motif = join_PFMs_to_dict([M[i] for i in d.models])
        diag_background = join_motif_background_frequencies([M[i] for i in d.models])

        # set up the matrix for p-value computation from combined scores
        if filter_options.filterby is FilterBy.TFMP:
            d.tfmp = tfmp.read_matrix(
                ' '.join(' '.join(str(i) for i in diag_motif[k]) for k in alphabet),
                bg=[diag_background[k] for k in alphabet]
            )

    # this accepts input in the JASPAR motif format
    # this maybe the conversion code (c) holds Bio.motifs.Motif ' '.join(' '.join(str(i) for i in c[k]) for k in 'ACGT')
    # checked with the Bio.motifs.Motif.format() for pfm and JASPAR and it appears to be correct

    if fasta.endswith('.gz'):
        fasta_handle = gzip.open(fasta, 'rt')
    else:
        fasta_handle = open(fasta, 'r')
    input_seqs = [s for s in SeqIO.parse(fasta_handle, 'fasta')]

    for seq_results in fimo_wrapper(motif_file, input_seqs, options=fimo_options):
        _ = analyze_output(seq_results, D, filter_options)

    all_sites = pd.concat([d.to_df() for d in D], ignore_index=True)
    all_sites.sort_values('score', ascending=False, inplace=True)

    # output section
    # todo: output some genome browser compatible format (gff for example)
    # todo: allow more outputs at once

    if vis_track_prefix != '':
        with open(vis_track_prefix + '.gff', 'w') as f1, open(vis_track_prefix + '.interact', 'w') as f2:
            for d in D:
                for i, j in d.to_vis_tracks():
                    f1.write(i+'\n')
                    f2.write(j+'\n')

        to_html(vis_track_prefix, fasta)

    if diag_match_output_file.endswith('xlsx'):
        all_sites.to_excel(diag_match_output_file, index=False)
    elif diag_match_output_file.endswith('tsv'):
        all_sites.to_csv(diag_match_output_file, sep='\t', index=False)
    else:
        # default to csv
        all_sites.to_csv(diag_match_output_file, index=False)


if __name__ == '__main__':
    from tests import *
    # test1()
    # test1_rc()
    # test2()
    # test2v2()
    # test2v3()
    # test2v4()
    # test3v1()
    pass




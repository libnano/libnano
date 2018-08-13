import sys
from pprint import pprint
from typing import (
    List,
    Dict,
    NamedTuple
)
from contextlib import contextmanager

import click

# update path as required
try:
    import libnano
except:
    # B. libnano not installed so update the path
    from os.path import (
        dirname,
        abspath
    )
    LIBNANO_PATH = dirname(dirname(dirname(abspath(__file__))))
    print(LIBNANO_PATH)
    sys.path.append(LIBNANO_PATH)

from libnano.ensemblrest import (
    lookUpID,
    lookUpSymbolList,
    overlap,
    LookUp,
    permittedSequences,
    Transcript
)
import libnano.ensemblrest as ger

from libnano.padlock import (
    generatePadlocks,
    writePadlocksToCSV,
    DEFAULT_PADLOCK_CONFIG
)

import libnano.barcode_tools as bt
BARCODE_LIST = bt.getBarcodeSet('34_3hd_5mer_v00')

THIS_FILE: str = __file__

GeneIds = NamedTuple('GeneIds', [('mouse', str), ('human', str)])

ALL_GENES: Dict[str, str] = {
'TUBB3':    GeneIds('ENSMUSG00000062380', 'ENSG00000258947'),
'SLC18A2':  GeneIds('ENSMUSG00000025094', 'ENSG00000165646'),
'SLC17A8':  GeneIds('ENSMUSG00000019935', 'ENSG00000179520'),
'GAD1':     GeneIds('ENSMUSG00000070880', 'ENSG00000128683'),
'GFAP':     GeneIds('ENSMUSG00000020932', 'ENSG00000131095'), # look at off target
'MBP':      GeneIds('ENSMUSG00000041607', 'ENSG00000197971'),
'TMEM119':  GeneIds('ENSMUSG00000054675', 'ENSG00000183160'),
'REST':     GeneIds('ENSMUSG00000029249', 'ENSG00000084093'),
'FOXO1':    GeneIds('ENSMUSG00000044167', 'ENSG00000150907'),
'FOXO3':    GeneIds('ENSMUSG00000048756', 'ENSG00000118689'),
'DISC1':    GeneIds('ENSMUSG00000043051', 'ENSG00000162946'),
'DRD2':     GeneIds('ENSMUSG00000032259', 'ENSG00000149295'),
'SST':      GeneIds('ENSMUSG00000004366', 'ENSG00000157005'),
'CALB1':    GeneIds('ENSMUSG00000028222', 'ENSG00000104327'),  # look at off target
'SNAP25':   GeneIds('ENSMUSG00000027273', 'ENSG00000132639'),
'SYN1':     GeneIds('ENSMUSG00000037217', 'ENSG00000008056'),
'PSEN2':    GeneIds('ENSMUSG00000010609', 'ENSG00000143801'),
'DAXX':     GeneIds('ENSMUSG00000002307', 'ENSG00000204209'),
'CDK5R1':   GeneIds('ENSMUSG00000048895', 'ENSG00000176749'),
'SLC22A7':  GeneIds('ENSMUSG00000067144', 'ENSG00000137204'),
'APP':      GeneIds('ENSMUSG00000022892', 'ENSG00000142192')
}

ALL_SYMBOLS: List[str] = sorted(ALL_GENES.keys())

@contextmanager
def disable_cache():
    ger.USE_CACHE = False
    yield
    print("cache on")
    ger.USE_CACHE = True

def useCache(is_on: bool):
    ger.USE_CACHE = is_on
# end def

def listExons(transcript_id: str):
    res: dict = lookUpID(transcript_id)
    msg = '%-10s%8d'
    exon = res.get('Exon')
    if exon is None:
        raise ValueError('{} is not an Ensembl Transcript ID'.format(transcript_id))
    for exon in res['Exon']:
        start, end  = exon['start'], exon['end']
        print(msg % (exon['id'],
                    end - start
                    ))
# end def

def listTranscriptsIdxs(species: str,
                        symbols: List[str],
                        filename: str = None):
    '''List the start and end indices of a transcript.
    This help to compare transcript of a gene to determine what location to
    target

    Correct for minimum start index to allow for easier visual comparison.

    Args:
        species:
        symbols:
    '''
    print("Listing transcript indices")
    res: dict = lookUpSymbolList(species, symbols)
    msg: str = '%-12s%-20s%8s%8s%8s%8s%8s%8s'
    print(msg % ('Name', 'EID', 'start', 'end', 'delta', 'elength', 'canon', 'strand'))
    msg: str = '%-12s%-20s%8d%8d%8d%8d%8s%8s'
    for symbol in symbols:
        item: dict = res[symbol]
        transcripts: dict = item['Transcript']
        # start_min: float = 1e20
        end_canon: int = 0
        for transcript in transcripts:
            start, end  = transcript['start'], transcript['end']
            if transcript['is_canonical'] == 1:
                end_canon = end
                start_canon: int = start
        if item['strand'] == 1:
            for transcript in transcripts:
                start, end  = transcript['start'], transcript['end']
                is_canonical = 'Yes' if transcript['is_canonical'] == 1 else 'No'
                elength: int = 0
                for exon in transcript['Exon']:
                    elength += abs(exon['end'] - exon['start'])
                print(msg % (
                                transcript['display_name'],
                                transcript['id'],
                                start - start_canon,
                                end - start_canon,
                                end - end_canon,
                                elength,
                                is_canonical,
                                'fwd'))
        else:
            for transcript in transcripts:
                start, end  = transcript['start'], transcript['end']
                is_canonical: str = 'Yes' if transcript['is_canonical'] == 1 else 'No'
                elength: int = 0
                for exon in transcript['Exon']:
                    elength += abs(exon['end'] - exon['start'])
                print(msg % (
                                transcript['display_name'],
                                transcript['id'],
                                end_canon - start,
                                end_canon - end,
                                start - start_canon,
                                elength,
                                is_canonical,
                                'rev'))
# end def

def listDetails(species: str,
                symbols: List[str],
                filename: str = None):
    '''List the details of a group of symbols for a species

    Currently::

        Symbol,
        # of protein coding transcripts,
        strand fwd/rev,
        gene length,
        length of mRNA of the canonical transcript
        # of variant in the 3' UTR

    Args:
        species:
        symbols:
        filename: Default is ``None``. This is unused
    '''
    print("Listing Details")
    res: dict = lookUpSymbolList(species, symbols)
    msg: str = '%-8s%6s%8s%8s%8s%8s'
    print(msg % ('Name', '#PCTs', 'strand', 'length', 'elength', '#vars'))
    msg: str = '%-8s%6d%8s%8d%8d%8d'
    for symbol in symbols:
        try:
            item: dict = res[symbol]
        except:
            print(res)
            raise
        length = item['end'] - item['start']

        elength: int = 0
        variant_count: int = -1
        for transcript in item['Transcript']:
            if transcript['is_canonical'] == 1:
                three_prime_utr_id: str = transcript['Exon'][-1]['id']
                variant_count = len(overlap(three_prime_utr_id))
                for exon in transcript['Exon']:
                    elength += exon['end'] - exon['start']
                break

        to_print = ( item['display_name'],
                    len(item['Transcript']),
                    'fwd' if item['strand'] == 1 else 'rev',
                    length,
                    elength,
                    variant_count)
        print(msg % to_print)
# end def

def designPadlocks( species: str,
                    symbols: List[str],
                    exon_index: int = -1,
                    arm_length: int = 20,
                    three_prime_delta: int = 10,
                    barcodes: List[str] = None,
                    do_print: bool = True,
                    filename: str = None):
    '''Get the canon transcript and design a set of padlocks that fall in the
    3' UTR of the transcript by default to be compatible with polyA targeted
    RT of the mRNA

    Args:
        species:
        symbols:
        exon_index: default is -1 (AKA the 3' UTR)
        arm_length: padlock arm length
        three_prime_delta: distance to heed from the three prime end
        filename: Default is ``None``. This is used to write a CSV of the results
    '''
    print("Designing Padlocks")
    out: dict = {symbol: [] for symbol in symbols}

    p_params: dict = DEFAULT_PADLOCK_CONFIG()

    if barcodes is None:
        barcodes = BARCODE_LIST

    arm_length_2X: int = 2*arm_length
    res: dict = lookUpSymbolList(species, symbols)
    barcode_idx: int = 0
    for symbol in symbols:
        lookup: LookUp = LookUp(res[symbol])
        for prospective_transcript in lookup.transcripts:
            if prospective_transcript['is_canonical'] == 1:
                canon_transcript = Transcript(prospective_transcript)
                break

        canon_transcript_id: str = canon_transcript.id
        exon_id: str = canon_transcript.Exon[exon_index]['id']

        p_seqs: List[List[Tuple[int, str]]] = permittedSequences(canon_transcript)
        exon_segments: List[Tuple[int, str]] = p_seqs[exon_index]
        # print(exon_segments)
        max_idx: int = sum([len(x[1]) for x in exon_segments])

        strand_dir: str = 'fwd' if lookup.is_fwd else 'rev'
        if do_print:
            print(symbol, len(exon_segments), strand_dir)
            msg = '%8s%8s%8s%8s%10s'
            print(msg % ('start', 'maxidx', 'seglen', 'addlen', 'g_index'))
            msg = '%8d%8d%8d%8d%10d'
        candidate_segments: List[Tuple[int, str]] = []
        idx: int = 0
        bound_idx: int = max_idx - three_prime_delta
        for g_index, segment in exon_segments:
            len_segment: int = len(segment)
            # print(len_segment, arm_length_2X)
            if len_segment > arm_length_2X:
                if idx + len_segment > bound_idx:
                    # Of the form:: len_segment - (idx + len_segment - bound_idx)
                    adjust_bound: int = bound_idx - idx
                    if adjust_bound < arm_length_2X:
                        continue
                    # print("Do adjust", adjust_bound,
                    #                     len_segment,
                    #                     bound_idx,
                    #                     idx + len_segment, idx)
                    if lookup.is_rev:
                        g_index = g_index - (len_segment - adjust_bound)
                    add_segment: str = segment[:adjust_bound]
                else:
                    add_segment: str = segment
                if do_print:
                    print(msg %     (idx, max_idx, len_segment,
                                    len(add_segment), g_index) )
                candidate_segments.append((g_index, add_segment))
            idx += len_segment
        # end for
        barcode: str = barcodes[barcode_idx]
        barcode_idx += 1

        pads: list = out[symbol]
        for g_index, segment in candidate_segments:
            pads += generatePadlocks(segment,
                                    canon_transcript_id,
                                    exon_id,
                                    strand_dir,
                                    [barcode],
                                    genome_idx=g_index,
                                    params=p_params)
            if do_print:
                pprint([x.padlock_seq for x in pads])

        '''Retry for No hits with higher allowable GC max content'''
        if len(pads) == 0:
            print(">> bumping GC content for %s" % symbol)
            bumped_params: dict = DEFAULT_PADLOCK_CONFIG()
            bumped_params['arm_gc_max'] = 0.65
            for g_index, segment in candidate_segments:
                pads += generatePadlocks(segment,
                                        canon_transcript_id,
                                        exon_id,
                                        strand_dir,
                                        [barcode],
                                        genome_idx=g_index,
                                        params=bumped_params)
                if do_print:
                    pprint([x.padlock_seq for x in pads])
        cannon_transcript = None    # clear the transcript
    # end for
    if filename is not None:
        writePadlocksToCSV(out, filename)
# end def

output_help = '''Kind of function to run:
    d: listDetails
    e: listExons
    i: listTranscriptsIdxs
    p: designPadlocks
'''
barcode_help = '''if genes are mapped to specific barcodes, only used when genes
                is specifically set. use the form:

                \b
                SYMBOL0 BARCODE0 SYMBOL1 BARCODE1 ...'''

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--output', '-o', default='d', help=output_help)
@click.option('--genes', '-g',
                default=None,
                type=str,
                help='Symbols of genes in quotes, like \"SST MBP\"')
@click.option('--dobarcode', '-b',
                is_flag=True,
                help=barcode_help)
@click.option('--ids', '-i',
                default=None,
                type=str,
                help='Symbols of genes in quotes, like \"SST MBP\"')
@click.option('--species', '-s', default='mouse', help='species, mouse or human')
@click.option('--filename', '-f', default=None, help='file name to write output')
def cli(genes, output, species, ids, dobarcode, filename):
    '''This script will provide information on a list of genes for a give species
    '''
    barcodes: List[str]
    if genes is None:
        genes = ALL_SYMBOLS
        barcodes = None
    else:
        genes = genes.split()
        if dobarcode:
            barcodes =  [x[i] for x in range(1, len(genes), 2)]
            genes =     [x[i] for x in range(0, len(genes), 2)]
    if ids is not None:
        ids = ids.split()
    if output == 'd':
        func = listDetails
    elif output == 'i':
        func = listTranscriptsIdxs
    elif output == 'p':
        designPadlocks(species, genes, barcodes=barcodes, filename=filename)
        return
    elif output == 'e':
        if ids is not None:
            for eid in ids:
                listExons(eid)
            return
        else:
            raise ValueError("Need to use the `--ids` option to specify a transcript")
    else:
        raise ValueError("unknown function command {}".format(output))

    if species not in ('mouse', 'human'):
        raise ValueError("unknown specieds {}".format(species))
    # print("going to run")
    func(species, genes, filename=filename)
# end def

if __name__ == '__main__':
    cli()
    # genes = ['GFAP', 'GAD1', 'MBP', 'CALB1']
    # barcodes = ['ACAGC', 'ACCAA', 'AATAC', 'CTAAG']
    # listDetails('mouse', genes)
    # designPadlocks('mouse', genes, barcodes=barcodes, filename='winners_v03.csv')
    # listTranscriptsIdxs('mouse', genes)
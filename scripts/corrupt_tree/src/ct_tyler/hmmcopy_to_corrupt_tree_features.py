from argparse import ArgumentParser
import pandas as pd


def get_args():
    desc = 'Generate corrupt-tree input from hmmcopy output'
    p = ArgumentParser(description=desc)

    p.add_argument('hmmcopy', help='hmmcopy copy number tsv file')
    p.add_argument('output', help='generated corrupt-tree input file')

    return p.parse_args()


def generate_bin_edge_ids(chrom, chrom_cns):
    starts = chrom_cns['start'].iloc[:-1].astype(str).values
    ends = chrom_cns['end'].iloc[:-1].astype(str).values
    bin_edge_ids = chrom + '_' + starts + '_' + ends
    #  starts = chrom_cns['start'].iloc[1:].astype(str).values
    #  ends = chrom_cns['end'].iloc[:-1].astype(str).values
    #  bin_edge_ids = chrom + '_' + ends + '_' + starts
    return bin_edge_ids


def convert_cn_to_breakpoints(cn):
    cn.sort_values(by=['chr', 'start', 'end'], inplace=True)
    
    
    pieces = []
    for chrom, group in cn.groupby('chr'):
        bin_edge_ids = generate_bin_edge_ids(chrom, group)
        bin_gap = (
            group['start'].iloc[1:].values - group['end'].iloc[:-1].values
        )

        group.set_index(['chr', 'start', 'end', 'width'], inplace=True)

        chrom_brkpt = (group.diff().iloc[1:] != 0).astype(int)
        chrom_brkpt['loci'] = bin_edge_ids
        chrom_brkpt.set_index('loci', inplace=True)

        #  chrom_brkpt = chrom_brkpt[bin_gap == 1]

        pieces.append(chrom_brkpt)

    return pd.concat(pieces)


def generate_corrupt_tree_input(cn):
    brkpt = convert_cn_to_breakpoints(cn)

    #  all_similar = brkpt.apply(lambda x: (x == x[0]).all(), axis=1)
    #  brkpt = brkpt[~all_similar]

    brkpt = brkpt.drop_duplicates()
    brkpt = brkpt.reset_index().melt(
        'loci', var_name='cells', value_name='tipInclusionProbabilities')

    return brkpt[['cells', 'loci', 'tipInclusionProbabilities']]

if __name__ == '__main__':
    argv = get_args()
    cn = pd.read_csv(argv.hmmcopy, header=0, dtype={'chr': str})
#     cn = pd.read_csv(argv.hmmcopy, sep='\t', dtype={'chr': str})
    brkpt = generate_corrupt_tree_input(cn)
    brkpt.to_csv(argv.output, index=False)

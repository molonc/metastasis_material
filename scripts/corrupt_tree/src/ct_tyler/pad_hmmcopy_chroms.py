from argparse import ArgumentParser
import pandas as pd
from scipy.stats import mode


def get_args():
    desc = 'Generate corrupt-tree input from hmmcopy output'
    p = ArgumentParser(description=desc)

    p.add_argument('input', help='hmmcopy copy number tsv file')
    p.add_argument('output', help='generated corrupt-tree input file')

    return p.parse_args()


def generate_chrom_pad(chrom, start, end, width, cells, cell_ploidy):
    pad = pd.DataFrame(
        columns=['chr', 'start', 'end', 'width'] + cells,
        data=[[chrom, start, end, width] + cell_ploidy]
    )
    print(pad.shape)
    return pad


def pad_chrom_cn(chrom_cn, cell_ploidies):
    cells = [
        x for x in chrom_cn.columns
        if x not in ['chr', 'start', 'end', 'width']
    ]

    chrom = chrom_cn.iloc[0]['chr']
    width = chrom_cn.iloc[0]['width']
    end_max = chrom_cn['end'].max()

    ploidy = cell_ploidies.loc[cells]['ploidy'].tolist()
    print('Add left padding...')
    # Tyler: -1 involve lots of errors
    # left_pad = generate_chrom_pad(chrom, -1, 0, width, cells, ploidy)
    left_pad = generate_chrom_pad(chrom, 0, 1, width, cells, ploidy)
    print('Add right padding...')
    right_pad = generate_chrom_pad(
        chrom, end_max + 1, end_max + 2, width, cells, ploidy)

    return pd.concat([left_pad, chrom_cn, right_pad])
    # return pd.concat([chrom_cn, right_pad])


def pad_hmmcopy_chroms(cn):
    print('Sorting')
    cn.sort_values(by=['chr', 'start', 'end'], inplace=True)

    cell_ploidy = (
        pd.melt(
            cn, id_vars=['chr', 'start', 'end', 'width'],
            var_name='cell_id', value_name='state'
        ).\
        groupby('cell_id')['state'].\
        agg(lambda x: mode(x)[0][0]).\
        to_frame().\
        rename(columns={'state': 'ploidy'})
    )
    print('Padding...')
    pieces = []
    for chrom, group in cn.groupby('chr'):
        pieces.append(pad_chrom_cn(group, cell_ploidy))
            
    return pd.concat(pieces, axis=0)

if __name__ == '__main__':
    argv = get_args()

    # cn = pd.read_csv(argv.input, sep='\t', dtype={'chr': str})
    # Hoa
    cn = pd.read_csv(argv.input, header=0, dtype={'chr': str})
    print('Filtered mtx: {0}{1}'.format(cn.shape[0],cn.shape[1]))
    print('Start padding')
    padded_cn = pad_hmmcopy_chroms(cn)
    print('Save data to file:')
    # padded_cn.to_csv(argv.output, sep='\t', index=False)
    padded_cn.to_csv(argv.output, index=False)

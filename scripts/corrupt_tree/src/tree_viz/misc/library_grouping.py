import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--library_ids', required=True)
    parser.add_argument('--groupings', required=True)
    parser.add_argument('--mainsites', required=True)
    parser.add_argument('--pdxids', required=True)
    parser.add_argument('--output', '-o', default="library_groupings.csv")

    args = vars(parser.parse_args())

    library_ids = args['library_ids'].split(",")
    groupings = args['groupings'].split(",")
    mainsites = args['mainsites'].split(",")
    pdxids = args['pdxids'].split(",")
    df = pd.DataFrame(list(zip(library_ids, groupings, mainsites,pdxids)), columns = ['library_labels', 'grouping','mainsite','pdxid'])

    df.to_csv(args['output'], index=None, header=True)


if __name__ == '__main__':
    main()
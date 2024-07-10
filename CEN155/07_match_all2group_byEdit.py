import pandas as pd
import argparse
import Levenshtein


def match_monomer_group(monomerID, merge, outfile):
    df_know = pd.read_table(merge, sep="\t", header=None, names=["monomerID", "chr", "group", "seq"])
    df_big_monomer = pd.read_table(monomerID, sep="\t", header=None, names=["monomerID", "seq"])
    o = open(outfile, 'w')
    for index, row in df_big_monomer.iterrows():
        monomer_ID = row["monomerID"]
        monomer_seq = row["seq"]
        min_value = 9999
        min_group = "None"
        for i, row_dict in df_know.iterrows():
            know_group = row_dict["group"]
            monomer_s = row_dict["seq"]
            Edit = Levenshtein.distance(str(monomer_seq), str(monomer_s))
            if Edit < min_value:
                min_group = know_group
                min_value = Edit
            else:
                continue
        o.write(str(monomer_ID) + '\t' + str(min_group) + '\t' + str(min_value) + '\t' + str(monomer_seq) + '\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i1', '--input1', help='Mgroup', dest='input1', required=True)
    parser.add_argument('-i2', '--input2', help='kmer', dest='input2', required=True)
    parser.add_argument('-o', '--output', help='output', dest='output', required=True)
    args = parser.parse_args()
    match_monomer_group(args.input1, args.input2, args.output)

if __name__ == '__main__':
    main()
import argparse
import pandas as pd
import numpy as np


def get_centromere_region(prefix, outfile):
    infile = str(prefix) + "_consensus.blast"
    df = pd.read_table(infile, sep='\t', header=None)
    df.columns = ["consensus", "Chr", "iden", "Length", "gap", "mismatch", "s1", "e1", "s2", "e2", "evalue", "score"]
    Chr_list = df.drop_duplicates(subset=['Chr'], keep='last', ignore_index=True)['Chr'].tolist()
    o = open(outfile, "a")
    for chr in Chr_list:

        df_chr = df[(df['Chr'] == str(chr)) & (df['Length'] >= 100)]
        for index, row in df_chr.iterrows():
            s2 = df_chr.loc[index, 's2']
            e2 = df_chr.loc[index, 'e2']
            if s2 > e2:
                df_chr.loc[index, 's2'] = int(e2)
                df_chr.loc[index, 'e2'] = int(s2)
        df_chr_sort = df_chr.sort_values(by='s2').reset_index()
        circle = 0
        line_num = 0
        next_line = 1
        while next_line < df_chr_sort.shape[0]:
            if (int(df_chr_sort.iloc[next_line, 10]) - int(df_chr_sort.iloc[line_num, 9]) + 1) <= (
                    int(sum((df_chr_sort.iloc[line_num:(next_line + 1), 4]))) + 150000):
                circle += 1
                line_num += 1
                next_line += 1
            else:
                start_line0 = line_num - circle
                end_line0 = next_line - 1
                if end_line0 - start_line0 > 1:
                    start0 = df_chr_sort.loc[start_line0, "s2"]
                    end0 = df_chr_sort.loc[end_line0, "e2"]
                    length0 = end0 - start0 + 1
                    o.write(str(prefix) + "\t" + str(chr) + "\t" + str(int(start0)) + "\t" + str(int(end0)) + "\t" + str(int(length0)) + "\n")
                    circle = 0
                    line_num += 1
                    next_line += 1
                else:
                    circle = 0
                    line_num += 1
                    next_line += 1
        start_line = line_num - circle
        end_line = next_line - 1
        start = df_chr_sort.loc[start_line, "s2"]
        end = df_chr_sort.loc[end_line, "e2"]
        length = end - start + 1
        o.write(str(prefix) + "\t" + str(chr) + "\t" + str(int(start)) + "\t" + str(int(end)) + "\t" + str(int(length)) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='centromere_regions', dest='input', required=True)
    parser.add_argument('-o', '--output', help='output', dest='output', required=True)
    args = parser.parse_args()
    get_centromere_region(args.input, args.output)

if __name__ == '__main__':
    main()

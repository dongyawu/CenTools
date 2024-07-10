import argparse
import pandas as pd

def get_count2pro(infile, materiallist, outfile):
    df_m = pd.read_table(materiallist, header=None, names=["ID"])
    m_list = df_m["ID"].to_list()
    df_all = pd.read_table(infile, sep="\t", index_col=0)
    df_all["sum"] = df_all.sum(axis=1)  ##  索引 0 对应第一行的总和，索引 1 对应第二行的总和....
    for index, row in df_all.iterrows():
        sum = int(row["sum"]) + 0.000001
        for m in m_list:
            df_all.loc[index, m] = round(int(row[m])/sum, 3)
    df_all.to_csv(outfile, index=False, sep="\t", header=True)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i1', '--input1', help='infile', dest='input1', required=True)
    parser.add_argument('-i2', '--input2', help='materiallist', dest='input2', required=True)
    parser.add_argument('-o', '--output', help='output', dest='output', required=True)
    args = parser.parse_args()
    get_count2pro(args.input1, args.input2, args.output)
if __name__ == '__main__':
    main()
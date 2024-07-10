import argparse
import pandas as pd

def get_count2pro(materiallist, outfile1, outfile2):
    df_m = pd.read_table(materiallist, header=None, names=["ID"])
    m_list = df_m["ID"].to_list()
    chr_list = ["Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12"]
    # chr_list = ["Chr01", "Chr02", "Chr03"]
    dir = "your file path/"
    df_first = pd.read_table(dir + chr_list[0] + "/" + chr_list[0] + "_61m_counttable.txt", sep="\t", index_col=0)
    for chr in chr_list[1:]:
        df_add = pd.read_table(dir + chr + "/" + chr + "_61m_counttable.txt", sep="\t", index_col=0)
        df_first = df_first.add(df_add, fill_value=0).astype(int)
    df_first["sum"] = df_first.sum(axis=1)
    df_first.to_csv(outfile1, sep="\t", header=True)
    df_merge = df_first.copy(deep=True)
    for index, row in df_merge.iterrows():
        sum = int(row["sum"]) + 0.000001
        for m in m_list:
            df_merge.loc[index, m] = round(int(row[m])/sum, 3)
    df_merge.to_csv(outfile2, sep="\t", header=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='infile', dest='input', required=True)
    parser.add_argument('-o1', '--output1', help='counttable', dest='output1', required=True)
    parser.add_argument('-o2', '--output2', help='proptable', dest='output2', required=True)
    args = parser.parse_args()
    get_count2pro(args.input, args.output1, args.output2)
if __name__ == '__main__':
    main()

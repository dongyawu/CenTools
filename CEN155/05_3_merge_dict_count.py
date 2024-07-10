import argparse
import pandas as pd

def merge_dict_count(materiallist, chr, outfile):
    df_m = pd.read_table(materiallist, header=None, names=["ID"])
    m_list = df_m["ID"].to_list()
    df_merge = pd.read_table(str(chr) + "_" + str(m_list[0]) + "_count.txt", header=None, names=["ID", "13-65"])
    for m in m_list[1:]:
        df = pd.read_table(str(chr) + "_" + str(m) + "_count.txt", header=None, names=["ID", m])
        df_merge = pd.merge(df_merge, df, on="ID", how="outer")
    df_merge.to_csv(outfile, index=False, sep="\t", header=True)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i1', '--input1', help='materiallist', dest='input1', required=True)
    parser.add_argument('-i2', '--input2', help='chr', dest='input2', required=True)
    parser.add_argument('-o', '--output', help='output', dest='output', required=True)
    args = parser.parse_args()
    merge_dict_count(args.input1, args.input2, args.output)
if __name__ == '__main__':
    main()

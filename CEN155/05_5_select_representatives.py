import argparse
import pandas as pd

def select_representatives(tablefile, monomerfile, outfile1, outfile2):
    df_monomer = pd.read_table(monomerfile, sep="\t", header=None, index_col=0)  ### 让其索引为M1 M2列
    df_all = pd.read_table(tablefile, sep="\t", index_col=0)
    df_select1 = df_all[df_all["sum"] >= 31]
    df_select3 = df_all[(df_all["sum"] >= 25) & (df_all["sum"] <= 30)]
    df_select2 = df_all[(df_all["sum"] >= 25) & (df_all["sum"] <= 30)].iloc[:, :-1]   ## 去掉最后一列
    df_list = [df_select1]
    for i in range(0, df_select2.shape[0]):
        row_to_count = df_select2.iloc[i]
        max_value = row_to_count.max()
        if max_value >= 0.2:
            df_add_i = pd.DataFrame(df_select3.iloc[i, :])
            df_add_i_T = df_add_i.T
            df_list.append(df_add_i_T)
    df_merge = pd.concat(df_list, axis=0)
    df_merge["sum"] = df_merge["sum"].astype(int)
    df_merge.to_csv(outfile1, header=True, sep="\t")
    monomer_ID_list = []
    for index, row in df_merge.iterrows():
        monomer_ID_list.append(index)   ##获得所有的索引
    selected_rows = df_monomer.loc[monomer_ID_list]
    selected_rows.to_csv(outfile2, sep="\t")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i1', '--input1', help='table', dest='input1', required=True)
    parser.add_argument('-i2', '--input2', help='monomer', dest='input2', required=True)
    parser.add_argument('-o1', '--output1', help='outfile1', dest='output1', required=True)
    parser.add_argument('-o2', '--output2', help='outfile2', dest='output2', required=True)
    args = parser.parse_args()
    select_representatives(args.input1, args.input2, args.output1, args.output2)
if __name__ == '__main__':
    main()

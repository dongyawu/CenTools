import argparse
import pandas as pd

def dict_count(monomerfile, prefix, outfile):
    dict_monomer_name = {}
    dict_monomer_count = {}
    o = open(outfile, 'w')
    df_monomer = pd.read_table(monomerfile, sep="\t", header=None)
    df_monomer.columns = ["ID", "seq"]
    for index, row in df_monomer.iterrows():
        ID = row["ID"]
        seq = row["seq"]
        dict_monomer_name[ID] = seq
        dict_monomer_count[ID] = 0
    df_repeat = pd.read_table("./allrepeats_CENtype/CENtype_repeats_" + str(prefix) + ".csv", sep=",")
    for index, row in df_monomer.iterrows():
        ID = row["ID"]
        df_CENT = df_repeat[(df_repeat["seq.name"].str.contains('Chr01')) & (df_repeat['class'].str.contains('CEN')) & (df_repeat['seq'] == str(dict_monomer_name[ID]))]
        dict_monomer_count[ID] += df_CENT.shape[0]
    for index, row in df_monomer.iterrows():
        ID = row["ID"]
        count = dict_monomer_count[ID]
        o.write(str(ID) + "\t" + str(count) + "\n")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i1', '--input1', help='monomer', dest='input1', required=True)
    parser.add_argument('-i2', '--input2', help='prefix', dest='input2', required=True)
    parser.add_argument('-o', '--output', help='output', dest='output', required=True)
    args = parser.parse_args()
    dict_count(args.input1, args.input2, args.output)
if __name__ == '__main__':
    main()

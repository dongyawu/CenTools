import pandas as pd
import argparse

def monomer_anno_add(prefix, match2group, outfile):
    df_m_anno = pd.read_table(str(prefix) + "_monomer_anno_final.txt", sep="\t") ##Chr	Start	End	Strand	MonomerID	Group	Sample	Type
    df_match = pd.read_table(match2group, sep="\t", header=None, names=["monomerID", "Group", "Editscore", "seq"])

    for index, row in df_m_anno.iterrows():
        mID = row["MonomerID"]
        if row["Group"] == "group_rare":
            df_new_group = df_match[df_match["monomerID"] == mID].reset_index(drop=True)
            if df_new_group.shape[0] != 0:
                if df_new_group.loc[0, "Editscore"] < 10:
                    new_group = "group" + str(df_new_group.loc[0, "Group"])
                    df_m_anno.at[index, "Group"] = new_group
                else:
                    continue
            else:
                continue
    df_new = df_m_anno[(df_m_anno["Group"] != "group_long") & (df_m_anno["Group"] != "group_short")]
    df_new.to_csv(outfile, sep="\t", header=True, index=False)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i1', '--input1', help='prefix', dest='input1', required=True)
    parser.add_argument('-i2', '--input2', help='match2group', dest='input2', required=True)
    parser.add_argument('-o', '--output', help='output', dest='output', required=True)
    args = parser.parse_args()
    monomer_anno_add(args.input1, args.input2, args.output)

if __name__ == '__main__':
    main()

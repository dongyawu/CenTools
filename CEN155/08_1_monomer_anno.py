import pandas as pd
import argparse

def monomer_annotation(prefix, allrepeatfile, outfile, out):    ### allrepeatfile = "TRASH_allrepeat_unalign_70m_rmdup.monomerID"
    df_cent = pd.read_table("centromere_position.bed", sep='\t')  ## Chr 	start	end
    df_all = pd.read_table("all.repeats.from." + str(prefix) + ".fasta.csv", sep=',', low_memory=False)
    # start,end,width,seq,strand,class,region.name,seq.name,edit.distance,repetitiveness

    df_all = df_all.applymap(str)
    df_all = df_all.applymap(remove_quotes)
    df_all = df_all.applymap(custom_conversion)

    df_monomer = pd.read_table(allrepeatfile, sep="\t", header=None)
    df_monomer.columns = ["ID", "seq"]  ## 这个monomer文件是提前制作好的，主要就是把序列的前面加上了一列代表monomerID的标识
    dict_monomer_name = {}
    for index, row in df_monomer.iterrows():  ### 分别把这些加入两个字典
        ID = row["ID"]
        seq = row["seq"]
        dict_monomer_name[seq] = ID

    df_group_type = merge(allrepeatfile, out)
    df_group_type.columns = ["monomer", "ID", "group", "seq"]

    df_mtype = pd.read_table("material.type", sep="\t", header=None, names=["material", "type"])
    o = open(outfile, 'a')
    Chr_list = ["Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12"]
    
    ########## get index dataframe ##########
    for each_chr in Chr_list:
        df_centromere = df_cent[df_cent['Chr'] == str(each_chr) + "_" + str(prefix)].reset_index(drop=True)
        cen_s = df_centromere.loc[0, "start"]
        cen_e = df_centromere.loc[0, "end"]
        df_all["start"] = df_all["start"].astype(float).astype(int)
        df_all["end"] = df_all["end"].astype(float).astype(int)
        df_all_CEN = df_all[(df_all['seq.name'] == str(each_chr) + "_" + str(prefix)) & (df_all['class'].str.contains('CEN')) & (df_all["start"] >= cen_s) & (df_all["end"] <= cen_e)].reset_index(drop=True)
        if df_all_CEN.shape[0] > 0:
            for index, row in df_all_CEN.iterrows():
                Start = row["start"]
                End = row["end"]
                Strand = row["strand"]
                seq = row["seq"]
                if seq in dict_monomer_name:
                    MonomerID = dict_monomer_name[seq]
                else:
                    MonomerID = "else"
                if seq in df_group_type["seq"].values:
                    df_temp = df_group_type[df_group_type["seq"] == seq].reset_index(drop=True)
                    Subfamily = df_temp.loc[0, "group"]
                else:
                    Subfamily = "100"
                df_temp2 = df_mtype[df_mtype["material"] == prefix].reset_index(drop=True)
                Type = df_temp2.loc[0, "type"]
                o.write(str(each_chr) + "_" + str(prefix) + "\t" + str(Start) + "\t" + str(End) + "\t" + str(Strand) + "\t" + str(MonomerID) + "\tgroup" + str(Subfamily) + "\t" + str(prefix) + "\t" + str(Type) + "\n")


def merge(allrepeatfile, out):
    df_subfamily = pd.read_table("monomer_subfamily.txt", sep="\t", header=None, names=["ID", "group"])
    df_ID = pd.read_table("ID_monomer6624.list", sep="\t", header=None, names=["monomer", "ID"])
    df_seq = pd.read_table(allrepeatfile, sep="\t", header=None, names=["monomer", "seq"])
    df_merge1 = df_ID.merge(df_subfamily, on="ID", how="right")
    df_merge2 = df_merge1.merge(df_seq, on="monomer", how="left")
    df_merge2.to_csv(out, sep="\t", index=False, header=False)
    return df_merge2

def remove_quotes(text):
    return text.replace('"', '')

def custom_conversion(cell_value):
    try:
        return int(cell_value)  # 尝试将字符串转换为整数
    except ValueError:
        return cell_value  # 如果无法转换，保留原字符串


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i1', '--input1', help='prefix', dest='input1', required=True)
    parser.add_argument('-i2', '--input2', help='allrepeatfile', dest='input2', required=True)
    parser.add_argument('-o', '--output', help='output', dest='output', required=True)
    parser.add_argument('-o2', '--output2', help='output2', dest='output2', required=True)
    args = parser.parse_args()
    monomer_annotation(args.input1, args.input2, args.output, args.output2)

if __name__ == '__main__':
    main()

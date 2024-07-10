import argparse
import pandas as pd

def zoom_in_centromeric_region_get_window_chip_100k(prefix, each_chr):
    df_mtype = pd.read_table("mtype.txt", sep="\t", header=None, names=["Material", "Type"])
    mtype = {}
    for index_r_type, row_type in df_mtype.iterrows():
        Material = row_type["Material"]
        Typee = row_type["Type"]
        mtype[Material] = Typee
    Type = mtype[prefix]
    chip_infile = "/public2/labmember/xielj/T2Trice_centromere/06_Methylation/04_ChIP/window_ChIP2/" + str(prefix) + "_log2ratio_2k.bdg"
    df_chippeak_region = pd.read_table("CENH3_position.txt", sep="\t")  ##Chr	CENH3S	CENH3E
    df_chr_chippeak = df_chippeak_region[df_chippeak_region["Chr"] == str(each_chr) + "_" + str(prefix)].reset_index(drop=True)
    outfile_chip = str(each_chr) + "_log2ratio_window2k_range500k.txt"
    if df_chr_chippeak.shape[0] != 0:
        chip_s = int(df_chr_chippeak.loc[0, "CENH3S"])
        chip_e = int(df_chr_chippeak.loc[0, "CENH3E"])
        with open(chip_infile, 'r') as inf:
            with open(outfile_chip, 'a') as outf2:
                for lines in inf:
                    text = lines.strip().split("\t")
                    chr = text[0]
                    win_s = int(text[1])
                    win_e = int(text[2])
                    prop = float(text[3])
                    type = text[4]

                    if chr == str(each_chr) + "_" + str(prefix):
                        if (win_s >= chip_s - 500000) and (win_e <= chip_e + 500000):
                            if prop < 0:
                                outf2.write(str(Type) + " " + str(prefix) + " " + str(each_chr) + ": " + str(chip_s-500000) + "-" + str(chip_e+500000) + "\t" + str(win_s-chip_s+500000) + "\t" + str(win_e-chip_s+500000) + "\t0\t" + str(type) + "\n")
                            else:
                                outf2.write(str(Type) + " " + str(prefix) + " " + str(each_chr) + ": " + str(chip_s-500000) + "-" + str(chip_e+500000) + "\t" + str(win_s-chip_s+500000) + "\t" + str(win_e-chip_s+500000) + "\t" + str(prop) + "\t" + str(type) + "\n")
        #### black box
        s_box = 0
        e_box = 1000000 + (chip_e - chip_s)
        with open(outfile_chip, 'a') as outf2:
            outf2.write(str(Type) + " " + str(prefix) + " " + str(each_chr) + ": " + str(chip_s-500000) + "-" + str(chip_e+500000) + "\t" + str(s_box) + "\t" + str(e_box) + "\t0\tbox\n")
        with open(str(each_chr)+ "_temp.txt", 'a') as outtemp:
            outtemp.write('"' + str(Type) + " " + str(prefix) + " " + str(each_chr) + ": " + str(chip_s-500000) + "-" + str(chip_e+500000) + '", ')

def annotation_tracks(prefix, each_chr):
    output = str(each_chr) + "_AnnoTracks_window2k_range500k.txt"
    o = open(output, 'a')
    df_mtype = pd.read_table("mtype.txt", sep="\t", header=None, names=["Material", "Type"])
    mtype = {}
    for index_r_type, row_type in df_mtype.iterrows():
        Material = row_type["Material"]
        Typee = row_type["Type"]
        mtype[Material] = Typee
    Type = mtype[prefix]
    df_chippeak_region = pd.read_table("CENH3_position.txt", sep="\t")  ##Chr	CENH3S	CENH3E
    df_chr_chippeak = df_chippeak_region[df_chippeak_region["Chr"] == str(each_chr) + "_" + str(prefix)].reset_index(drop=True)
    if df_chr_chippeak.shape[0] != 0:
        chip_s = int(df_chr_chippeak.loc[0, "CENH3S"])
        chip_e = int(df_chr_chippeak.loc[0, "CENH3E"])

        Dir = "/public2/labmember/xielj/T2Trice_centromere/10_synteny_centromere/3_SynPanCent/"
        df_newID = pd.read_table(str(Dir) + "00_AllAnnotation/" + str(prefix) + "_full_annotation.rmdup.txt", sep="\t", header=None, names=["chr", "start", "end", "strand", "id", "group"])
        df_newID_chr = df_newID[df_newID["chr"] == str(each_chr) + "_" + str(prefix)].reset_index(drop=True)
        for index, row in df_newID_chr.iterrows():
            s = int(row["start"])
            e = int(row["end"])
            group = row["group"]
            id = row["id"]

            if (s >= chip_s - 500000) and (e <= chip_e + 500000):
                if str(group).find("+") != -1: ##monomer
                    group_r1 = str(id).split(".")[0]
                    group_r2 = str(group_r1).split("_")[3]
                    o.write(str(Type) + " " + str(prefix) + " " + str(each_chr) + ": " + str(chip_s-500000) + "-" + str(chip_e+500000) + "\t" + str(id) + "\t" + str(s-chip_s+500000) + "\t" + str(e-chip_s+500000) + "\t-1\t" + str(group_r2) + "\n")
                elif str(group).find("group") != -1:
                    group_r1 = str(id).split(".")[0]
                    group_r2 = str(group_r1).split("_")[3]
                    o.write(str(Type) + " " + str(prefix) + " " + str(each_chr) + ": " + str(chip_s-500000) + "-" + str(chip_e+500000) + "\t" + str(id) + "\t" + str(s-chip_s+500000) + "\t" + str(e-chip_s+500000) + "\t-1\t" + str(group_r2) + "\n")
                elif group == "Interval":
                    o.write(str(Type) + " " + str(prefix) + " " + str(each_chr) + ": " + str(chip_s-500000) + "-" + str(chip_e+500000) + "\t" + str(id) + "\t" + str(s-chip_s+500000) + "\t" + str(e-chip_s+500000) + "\t-5\t" + str(group) + "\n")
                else:
                    pass
        df_TE = pd.read_table("/public2/labmember/xielj/T2Trice_centromere/00_data/CENH3_ChIPseq/ChIPpeak/Omit_nonCEN155_ChIP/anno_gff.bed", sep="\t", header=None, names=["Chr", "start", "end", "anno"])
        df_TE_chr = df_TE[df_TE["Chr"] == str(each_chr) + "_" + str(prefix)].reset_index(drop=True)
        for indexx, roww in df_TE_chr.iterrows():
            s = int(roww["start"])
            e = int(roww["end"])
            anno = roww["anno"]
            if (s >= chip_s - 500000) and (e <= chip_e + 500000):
                if str(anno).find("SZ-22") != -1:
                    groupnew = "sz22"
                    o.write(str(Type) + " " + str(prefix) + " " + str(each_chr) + ": " + str(chip_s-500000) + "-" + str(chip_e+500000) + "\t" + str(id) + "\t" + str(s - chip_s+500000) + "\t" + str(e - chip_s+500000) + "\t-2\t" + str(groupnew) + "\n")
                elif str(anno).find("CRM") != -1:
                    groupnew = "CRM"
                    o.write(str(Type) + " " + str(prefix) + " " + str(each_chr) + ": " + str(chip_s-500000) + "-" + str(chip_e+500000) + "\t" + str(id) + "\t" + str(s - chip_s+500000) + "\t" + str(e - chip_s+500000) + "\t-3\t" + str(groupnew) + "\n")
                elif str(anno).find("RIRE7") != -1:
                    groupnew = "RIRE7"
                    o.write(str(Type) + " " + str(prefix) + " " + str(each_chr) + ": " + str(chip_s-500000) + "-" + str(chip_e+500000) + "\t" + str(id) + "\t" + str(s - chip_s+500000) + "\t" + str(e - chip_s+500000) + "\t-4\t" + str(groupnew) + "\n")
                else:
                    pass

        s_box = 0
        e_box = 1000000 + (chip_e - chip_s)
        o.write(str(Type) + " " + str(prefix) + " " + str(each_chr) + ": " + str(chip_s-500000) + "-" + str(chip_e+500000) + "\tBOXid1\t" + str(s_box) + "\t" + str(e_box) + "\t-1\tbox\n")
        o.write(str(Type) + " " + str(prefix) + " " + str(each_chr) + ": " + str(chip_s-500000) + "-" + str(chip_e+500000) + "\tBOXid2\t" + str(s_box) + "\t" + str(e_box) + "\t-2\tbox\n")
        o.write(str(Type) + " " + str(prefix) + " " + str(each_chr) + ": " + str(chip_s-500000) + "-" + str(chip_e+500000) + "\tBOXid3\t" + str(s_box) + "\t" + str(e_box) + "\t-3\tbox\n")
        o.write(str(Type) + " " + str(prefix) + " " + str(each_chr) + ": " + str(chip_s-500000) + "-" + str(chip_e+500000) + "\tBOXid2\t" + str(s_box) + "\t" + str(e_box) + "\t-4\tbox\n")
        o.write(str(Type) + " " + str(prefix) + " " + str(each_chr) + ": " + str(chip_s-500000) + "-" + str(chip_e+500000) + "\tBOXid3\t" + str(s_box) + "\t" + str(e_box) + "\t-5\tbox\n")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--prefix', help='prefix', dest='prefix', required=True)
    parser.add_argument('-c', '--chr', help='chr', dest='chr', required=True)
    args = parser.parse_args()
    zoom_in_centromeric_region_get_window_chip_100k(args.prefix, args.chr)
    annotation_tracks(args.prefix, args.chr)

if __name__ == '__main__':
    main()
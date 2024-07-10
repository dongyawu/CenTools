import os
import pandas as pd
import argparse

def get_met_ONT_prop(prefix, chr):
    Dir = "/your file path/"
    ONTmet_type = str(Dir) + str(chr) + "_" + str(prefix) + "_ONTmet_type_filter.txt"
    

    #### the length of choromosome ####
    df_chr_length = pd.read_table(str(prefix) + ".fasta.fai", sep='\t', header=None)
    df_chr_length.columns = ["chr", "length", "a", "b", "c"]
    MAX = df_chr_length.loc[df_chr_length["chr"] == str(chr) + "_" + str(prefix), "length"].iloc[0]

    outfile = str(prefix) + "_window_ONTmet_ChIPLOG2_10k.txt"
    outf = open(outfile, 'a')
    window = 10000
    s = 0
    e = 0 + window

    CHH_all = CHG_all = CG_all = 1
    CHH_met = CHG_met = CG_met = 0
    while s < MAX:
        with open(ONTmet_type, 'r') as inf:
            for lines in inf:
                line = lines.strip()
                text = line.split("\t")
                Pos = int(text[1])
                Prop = float(text[2])
                Type = text[3]
                if Pos <= e and Pos >= s:
                    if Type == "CHH":
                        CHH_all += 1
                        if Prop > 0:
                            CHH_met += 1
                    elif Type == "CHG":
                        CHG_all += 1
                        if Prop > 0:
                            CHG_met += 1
                    elif Type == "CG":
                        CG_all += 1
                        if Prop > 0:
                            CG_met += 1
                    else:
                        pass

        outf.write(str(chr) + "_" + str(prefix) + "\t" + str(s) + "\t" + str(e) + "\t" + str(CHH_met / CHH_all) + "\tCHH\n")
        outf.write(str(chr) + "_" + str(prefix) + "\t" + str(s) + "\t" + str(e) + "\t" + str(CHG_met / CHG_all) + "\tCHG\n")
        outf.write(str(chr) + "_" + str(prefix) + "\t" + str(s) + "\t" + str(e) + "\t" + str(CG_met / CG_all) + "\tCG\n")
        CHH_all = CHG_all = CG_all = 1
        CHH_met = CHG_met = CG_met = 0
        s += window
        e += window



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i1', '--input1', help='prefix', dest='input1', required=True)
    parser.add_argument('-i2', '--input2', help='chr', dest='input2', required=True)
    args = parser.parse_args()
    get_met_ONT_prop(args.input1, args.input2)

if __name__ == '__main__':
    main()
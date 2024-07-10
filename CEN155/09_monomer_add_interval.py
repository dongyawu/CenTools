import pandas as pd
import argparse

def check_region_gaps(prefix, each_chr):
    Dir = "/your file path/"
    infile = str(Dir) + str(each_chr) + "_" + str(prefix) + "_full_annotation.txt"
    df_all = pd.read_table(infile, sep="\t") ## "Chr", "Start", "End", "Strand", "ID", "Group"
    df_all_m = df_all[(df_all["Chr"] == str(each_chr) + "_" + str(prefix))].reset_index(drop=True)
    print(df_all_m.shape[0])
    line_total = df_all_m.shape[0] - 1
    outfile = str(Dir) + str(each_chr) + "_" + str(prefix) + "_region_omit.txt"
    with open(outfile, 'w') as outf:
        for index, row in df_all_m.iterrows():
            Chr = row["Chr"]
            ID_front = row["ID"]
            End = row["End"]
            if index > 0 and index < line_total:
                nest_index = int(index) + 1
                Start_next = df_all.loc[nest_index, "Start"]
                ID_next = df_all.loc[nest_index, "ID"]
                if int(End) + 1 == int(Start_next):
                    pass
                else:
                    length = int(Start_next) - int(End)
                    outf.write(str(Chr) + "\t" + str(End) + "\t" + str(Start_next) + "\t" + str(ID_front) + "\t" + str(ID_next) + "\t" + str(length) + "\n")

def monomer_anno_add(prefix, each_chr, outfile):
    Dir = "/your file path/"
    infile = str(Dir) + str(each_chr) + "_" + str(prefix) + "_full_annotation.txt"
    interval = str(Dir) + str(each_chr) + "_" + str(prefix) + "_region_omit.txt"
    tempfile = str(each_chr) + "_" + str(prefix) + "_temp_annotation.txt"
    os.system("cat %s %s > %s" % (infile,interval, tempfile))
    os.system("less %s | sort -k 2 | sort -k 1 >> %s" % (tempfile,outfile))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--prefix', help='prefix', dest='prefix', required=True)
    parser.add_argument('-c', '--Chr', help='Chr', dest='Chr', required=True)
    parser.add_argument('-o', '--output', help='output', dest='output', required=True)
    args = parser.parse_args()
    check_region_gaps(args.prefix, args.Chr)
    monomer_anno_add(args.prefix, args.Chr, args.output)

if __name__ == '__main__':
    main()

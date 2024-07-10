import argparse
import os
import re
import pandas as pd

def get_new_repeat(dir1, dir2):
    for file in os.listdir(dir1):
        material = re.search(r"all.repeats.from.(.+).fasta.csv", file).group(1)
        outfile = dir2 + "/CENtype_repeats_" + material + ".csv"
        repeat_file = dir1 + "/" + file
        df_repeat = pd.read_table(repeat_file, sep=',', low_memory=False)
        df_repeat = df_repeat.applymap(str)
        df_repeat = df_repeat.applymap(remove_quotes)
        df_repeat = df_repeat.applymap(custom_conversion)
        df_CENT = df_repeat[df_repeat['class'].str.contains('CEN') & (df_repeat['width'] >= 140) & (df_repeat['width'] <= 170)]
        df_CENT.to_csv(outfile, header=True, index=False)

def remove_quotes(text):
    return text.replace('"', '')
def custom_conversion(cell_value):
    try:
        return int(cell_value)
    except ValueError:
        return cell_value

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='monomer', dest='input', required=True)
    parser.add_argument('-o', '--output', help='output', dest='output', required=True)
    args = parser.parse_args()
    get_new_repeat(args.input, args.output)
if __name__ == '__main__':
    main()
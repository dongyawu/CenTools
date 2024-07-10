import pandas as pd
import argparse

def get_consensus(infile, outfile):
    df_Summary = pd.read_table(infile, sep=",") #name start end ave.score most.freq.value.N consensus.primary consensus.count width fasta.name class consensus.secondary repeats.identified
    consensus_repeat_max_num = df_Summary[(df_Summary['most.freq.value.N'] > 140) & (df_Summary['most.freq.value.N'] < 170)]['repeats.identified'].max()

    repeat_length = int(df_Summary[df_Summary['repeats.identified'] == consensus_repeat_max_num]['most.freq.value.N'].iloc[0])
    seq = df_Summary.loc[df_Summary['repeats.identified'] == consensus_repeat_max_num]['consensus.secondary'].iloc[0]
    o = open(outfile, 'w')
    o.write("name,length,seq\n")
    o.write("CEN" + str(repeat_length) + "," + str(repeat_length) + "," + str(seq))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='summary_repeats', dest='input', required=True)
    parser.add_argument('-o', '--output', help='output', dest='output', required=True)
    args = parser.parse_args()
    get_consensus(args.input, args.output)

if __name__ == '__main__':
    main()

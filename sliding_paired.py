#!/usr/bin/env python

"""A script to detect regions >60 kb no spannned by paired-end BESs"""

from sys import argv, exit
import pandas as pd

#original dataframe
ori_df = pd.read_csv(argv[1], sep="\t", header=None, names=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])

window = 10000
step = 1000
start = 1
cycle = 1

scf_ls = ori_df.seqid.unique()

for scf in scf_ls:
    size = int(scf.split("|")[1].replace("size", ""))
    spanning_remind = 1
    start_remind = 0
    end_remind = 0

    for i in range(0, size + 1, step):
        end = start + window
        """Paired-end BESs spanning, i.e. start position in gff3 is less than start of sliding window,
        and end position in gff3 is greater than start of sliding window"""
        enough_paired = ori_df[(ori_df.attributes.str.contains("Paired_BES")) & (ori_df.seqid == scf)]

        #Checks if there are enough paired-end BESs to perform the analysis
        if len(enough_paired["start"]) == 0:
            break

        spanning_df = ori_df[(ori_df.start < start) & (ori_df.end > end) & (ori_df.attributes.str.contains("Paired_BES")) & (ori_df.seqid == scf)]

        spanning_nu = len(spanning_df["start"])

        if spanning_nu == 0 and spanning_remind != 0:
            start_remind = start

        #off recording
        if spanning_nu != 0 and spanning_remind == 0:
            end_remind = end
            misa_size = end_remind - start_remind

            #Cut off to delete small regions, check value
            if misa_size > 60000:
                print("%s\t%d\t%d\tMisa%d" % (scf, start_remind, end_remind, cycle))
                cycle += 1

        #After the cycle
        spanning_remind = spanning_nu

        start += step

        #until sliding window end is no greater than scaffold end
        ### DEBUG: with this structure, it does not do the last window
        if end > size:
            start = 1
            end_remind = end

            misa_size = end_remind - start_remind

            if misa_size > 60000:
                print("%s\t%d\t%d\tMisa%d" % (scf, start_remind, end_remind, cycle))
                cycle += 1

            break

#!/home/cintia/anaconda3/bin/python

"""A script to detect the regions with more than 2 consistent unpaired-end BESs"""

from sys import argv, exit
import pandas as pd

#original dataframe
ori_df = pd.read_csv(argv[1], sep="\t", header=None, names=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])

window = 100000
step = 10000
start = 1
end = start + window

scf_ls = ori_df.seqid.unique()

for scf in scf_ls:
    size = int(scf.split("|")[1].replace("size", ""))

    #empty misassembly dictionary
    misa_dict = {}

    orientation = ["+", "-"]

    for orient in orientation:
        #sliding window
        for i in range(0, size + 1, step):
            #A dataframe with the aligned unpaired-end BESs contained in the window
            containing = ori_df[(ori_df.start > start) & (ori_df.end < end) & (ori_df.attributes.str.contains("Unpaired_BES")) & (ori_df.seqid == scf) & (ori_df.strand == orient)]

            #List of partner BESs
            points_ls = [partner_bes.split(";")[3].replace("Points=", "") for partner_bes in containing.attributes]
            #List of unique partner BESs
            uq_points_ls = set(points_ls)

            for item in uq_points_ls:
                count = points_ls.count(item)
                if count > 2 and item not in misa_dict:
                    misa_dict[item] = [start, end, min(containing["start"]), max(containing["end"]), orient]
                elif count > 2 and item in misa_dict and end > misa_dict[item][1] and misa_dict[item][4] == orient:
                    misa_dict[item][1], misa_dict[item][3] = end, max(containing["end"])
                elif count > 2 and item in misa_dict and end > misa_dict[item][1] and misa_dict[item][4] != orient:
                    print("%s\t%s\t%s\t%s_bes_%s" % (scf, misa_dict[item][2], misa_dict[item][3], misa_dict[item][4].replace("+", "posi").replace("-", "nega"), item))
                    misa_dict[item] = [start, end, min(containing["start"]), max(containing["end"]), orient]

            #until sliding window end is no greater than scaffold end
            ### DEBUG: with this structure, it does not do the last window
            if end > size:
                start = 1
                end = start + window
                break
            else:
                start += step
                end = start + window

    #Write results
    if len(misa_dict) > 0:
        for item in misa_dict:
            print("%s\t%s\t%s\t%s_bes_%s" % (scf, misa_dict[item][2], misa_dict[item][3], misa_dict[item][4].replace("+", "posi").replace("-", "nega"), item))

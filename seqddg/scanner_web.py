#!/usr/bin/python

#By Sun Jinyuan, 2021

#To generate single point mutations of a given sequence
#Using MMseq2 API for MSA generation


from utilities import parsermodule
#from utilities import HHsearch
from utilities import MMseq2
from utilities import potts
from utilities import feature_maker
from tensorflow.keras.models import load_model
import numpy as np


parser = parsermodule.get_parser()
args = parser.parse_args()
seqfilename = args.sequence
#a3mfilename = args.a3mfile
#path_to_database = args.database
#iter_num = args.iteration_number
#num_threads = args.num_threads


outfilename = seqfilename.split(".")[0]+".scan.txt"


#if a3mfilename == None:
    #a3mfilename = HHsearch.hhsearch(seqfilename, iter_num, path_to_database, num_threads)
a3mfilename = MMseq2.run_mmseqs(seqfilename)


w_out, v_out = potts.potts(a3mfilename)

#feature = feature_maker.feature()
#al_dd = feature.get_al_dd()
#print(seqfilename)
#fasta=seqfilename
#msa = feature_maker.get_MSA()
seq = feature_maker.getSeq(seqfilename)
msa = feature_maker.get_MSA(a3mfilename)


msaddg = load_model("model/SeqDDG.h5", compile = False)
msa_C = -0.5083466

def run(wpm,msa,v_out,w_out,seq):
    static_features = feature_maker.getStatic(wpm)
    statical_features = feature_maker.getStatical(msa,v_out,w_out,seq,wpm)
    feature_list = static_features + statical_features
    return feature_list

def gen_mut_scan(seq):
    i = 0
    mut_list = []
    AA_list = ["Q",
               "W",
               "E",
               "R",
               "T",
               "Y",
               "I",
               "P",
               "A",
               "S",
               "D",
               "F",
               "G",
               "H",
               "K",
               "L",
               "C",
               "V",
               "N",
               "M"]
    for res in seq:
        i += 1
        for AA in AA_list:
            if AA == res:
                continue
            else:
                mut_list.append(res+'_'+str(i)+'_'+AA)
    return mut_list

mut_list = gen_mut_scan(seq)


with open(outfilename,"w+") as of:
    of.write("mutation\tscore\n")
    of.close()
    for mut in mut_list:
        wpm = mut.split("_")
        #wpm = mutation.split("_")
        if wpm[0] and wpm[2] in "ARNDCQEGHILKMFPSTWYV":
            if seq[int(wpm[1])-1] == wpm[0]:
                feature_list = run(wpm,msa,v_out,w_out,seq)
                with open(outfilename,"a+") as of:
                    of.write(mut+"\t"+str(msaddg.predict(np.array([feature_list]))[0][0]-msa_C)+"\n")
                    of.close()
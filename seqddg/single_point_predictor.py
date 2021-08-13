#!/usr/bin/python

#By Sun Jinyuan, 2021

from utilities import parsermodule
from utilities import HHsearch
from utilities import potts
from utilities import feature_maker
from tensorflow.keras.models import load_model
import numpy as np


parser = parsermodule.get_parser()
args = parser.parse_args()
seqfilename = args.sequence
a3mfilename = args.a3mfile
path_to_database = args.database
iter_num = args.iteration_number
num_threads = args.num_threads
mutation = args.mutation
mutation_list = args.mutation_list

outfilename = seqfilename.split(".")[0]+".predicted.txt"
print("output at: "+outfilename)

if a3mfilename == None:
    a3mfilename = HHsearch.hhsearch(seqfilename, iter_num, path_to_database, num_threads)



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

if __name__ == '__main__':
    if mutation != None:
        with open(outfilename,"w+") as of:
            of.write("mutation\tscore\n")
            of.close()
        wpm = mutation.split("_")
        if wpm[0] and wpm[2] in "ARNDCQEGHILKMFPSTWYV":
            if seq[int(wpm[1])-1] == wpm[0]:
                feature_list = run(wpm,msa,v_out,w_out,seq)
                with open(outfilename,"a+") as of:
                    of.write(mutation+"\t"+str(msaddg.predict(np.array([feature_list]))[0][0]-msa_C)+"\n")
                    of.close()
                #print(msaddg.predict(np.array([feature_list]))[0][0]-msa_C)
    else:
        mtfile = open(mutation_list)
        with open(outfilename,"w+") as of:
            of.write("mutation\tscore\n")
            of.close()
        for line in mtfile:
            wpm = line.strip().split("_")
            #wpm = mutation.split("_")
            if wpm[0] and wpm[2] in "ARNDCQEGHILKMFPSTWYV":
                if seq[int(wpm[1])-1] == wpm[0]:
                    feature_list = run(wpm,msa,v_out,w_out,seq)
                    print(line.strip()+"\t"+str(msaddg.predict(np.array([feature_list]))[0][0]-msa_C))
                    with open(outfilename,"a+") as of:
                        of.write(line.strip()+"\t"+str(msaddg.predict(np.array([feature_list]))[0][0]-msa_C)+"\n")
                        of.close()


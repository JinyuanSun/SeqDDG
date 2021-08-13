#!/usr/bin/python

#By Sun Jinyuan, 2021


import pandas as pd
import numpy as np


class feature:
    __a = """AAa H V P IP HT ST GSI F0 F1 F2 C B EC
    A -0.171 −0.677 −0.680 −0.170 0.900 −0.476 −0.350 −0.044 −0.234 −0.269 0.587 −0.099 0.829
    D -0.767 −0.281 −0.417 −0.900 −0.155 −0.635 −0.213 −0.103 0.900 0.014 −0.475 −0.082 0.247
    C 0.508 −0.359 −0.329 −0.114 −0.652 0.476 −0.140 −0.642 −0.773 −0.035 −0.433 0.094 −0.388
    E -0.696 −0.058 −0.241 −0.868 0.900 −0.582 −0.230 0.347 0.480 0.021 −0.900 0.105 0.565
    F 0.646 0.412 0.373 −0.272 0.155 0.318 0.363 −0.863 −0.504 −0.113 −0.673 0.721 0.035
    G -0.342 −0.900 −0.900 −0.179 −0.900 −0.900 −0.900 0.701 0.527 −0.050 0.378 −0.900 0.829
    H −0.271 0.138 0.110 0.195 −0.031 −0.106 0.384 −0.480 −0.186 −0.255 −0.297 0.115 −0.088
    I 0.652 −0.009 −0.066 −0.186 0.155 0.688 0.900 −0.332 −0.662 −0.411 −0.288 0.879 −0.900
    K −0.889 0.163 0.066 0.727 0.279 −0.265 −0.088 0.339 0.844 0.900 −0.375 0.317 0.547
    L 0.596 −0.009 −0.066 −0.186 0.714 −0.053 0.213 −0.590 −0.115 −0.064 −0.288 0.879 0.865
    M 0.337 0.087 0.066 −0.262 0.652 −0.001 0.110 −0.738 −0.900 −0.893 −0.205 0.370 0.724
    N −0.674 −0.243 −0.329 −0.075 −0.403 −0.529 −0.213 0.516 0.242 0.000 −0.166 0.031 0.265
    P 0.055 −0.294 −0.900 −0.010 −0.900 0.106 0.247 0.059 0.868 0.014 0.900 0.487 0.212
    Q −0.464 −0.020 −0.110 −0.276 0.528 −0.371 −0.230 0.870 0.416 −0.319 −0.403 0.192 0.529
    R −0.900 0.466 0.373 0.900 0.528 −0.371 0.105 −0.066 0.416 −0.206 0.430 0.175 −0.106
    S −0.364 −0.544 −0.637 −0.265 −0.466 −0.212 −0.337 0.900 0.575 −0.050 −0.024 −0.300 0.600
    T −0.199 −0.321 −0.417 −0.288 −0.403 0.212 0.402 0.192 0.599 0.028 −0.212 0.323 0.406
    V 0.331 −0.232 −0.285 −0.191 −0.031 0.900 0.677 −0.480 −0.385 −0.120 −0.127 0.896 0.794
    W 0.900 0.900 0.900 −0.209 0.279 0.529 0.479 −0.900 −0.464 −0.900 −0.074 0.900 0.900
    Y 0.188 0.541 0.417 −0.274 −0.155 0.476 0.363 −0.634 −0.361 −0.659 −0.738 0.546 0.582"""

    vloumn = {"A":88.6,"R":173.4,"N":114.1,"D":111.1,"C":108.5,
                "Q":143.8,"E":138.4,"G":60.1,"H":153.2,"I":166.7,
                "L":166.7,"K":168.6,"M":162.9,"F":189.9,"P":112.7,
                "S":89.0,"T":116.1,"W":227.8,"Y":193.6,"V":140.0}

    hydropathy_index = {"R":-2.5,"K":-1.5,"D":-0.9,"Q":-0.85,"N":-0.78,
                          "E":-0.74,"H":-0.4,"S":-0.18,"T":-0.05,"P":0.12,
                          "Y":0.26,"C":0.29,"G":0.48,"A":0.62,"M":0.64,
                          "W":0.81,"L":1.1,"V":1.1, "F":1.2,"I":1.4}

    def get_al_dd(self):

        al_dd = {"A":{},"D":{},"C":{},"E":{},"F":{},
                 "G":{},"H":{},"I":{},"K":{},"L":{},
                 "M":{},"N":{},"P":{},"Q":{},"R":{},
                 "S":{},"T":{},"V":{},"W":{},"Y":{}}

        for line in self.__a.replace("−","-").split("\n"):
            lst = line.split()
            if lst[0] == "AAa":
                name_lst = lst[1:]
            else:
                al_dd[lst[0]] = dict(zip(name_lst,lst[1:]))

        #al_df = pd.DataFrame(al_dd).T
        return al_dd

feature = feature()
al_dd = feature.get_al_dd()
vloumn = feature.vloumn
hydropathy_index = feature.hydropathy_index


alphabet = "ARNDCQEGHILKMFPSTWYV-"
states = len(alphabet)
a2n = {}
for a,n in zip(alphabet,range(states)):
    a2n[a] = n  

def get_MSA(seqa3m):
    import string
    rm_lc = str.maketrans(dict.fromkeys(string.ascii_lowercase))
    ali_dict = {}
    a3mfile = open(seqa3m)
    for line in a3mfile:
        if line[0] == ">":
            head = line.strip()
            ali_dict[head] = ""
            #print(line.strip())
        else:
            line = line.translate(rm_lc).strip()
            lst = []
            for x in line:
                lst.append(x)
            ali_dict[head] = lst
            #print(lst)
    df = pd.DataFrame(ali_dict)
    seq_length = df.shape[0]
    seq_num = df.shape[1]
    df_msa = df.T
    
    alphabet = "ARNDCQEGHILKMFPSTWYV-"
    p_msa = np.zeros([21,seq_length])
    for p in range(seq_length):
        for num,AA in enumerate(alphabet):
            if AA != "-":
                count = df_msa[p].value_counts().get(AA)
                if count == None:
                    p_msa[num][p] = 0
                else:
                    p_msa[num][p] = count/seq_num
                    #print(num,p,count/seq_num)
        p = p + 1
        
    return p_msa
def getSeq(fasta):
    fasta = open(fasta)
    seq = ''
    for line in fasta:
        if line[0] == ">":
            continue
        else:
            seq += line.strip().replace("X","A")
    return seq

def getStatic(wpm):
    dP = 0.5*(-float(al_dd[wpm[0]]['P'])+float(al_dd[wpm[2]]['P']))
    dF1 = -0.5*(-float(al_dd[wpm[0]]['F1'])+float(al_dd[wpm[2]]['F1']))
    dB = 0.5*(-float(al_dd[wpm[0]]['B'])+float(al_dd[wpm[2]]['B']))
    dST = 0.5*(-float(al_dd[wpm[0]]['ST'])+float(al_dd[wpm[2]]['ST']))
    dHT = 0.5*(-float(al_dd[wpm[0]]['HT'])+float(al_dd[wpm[2]]['HT']))
    dV = -0.007*(vloumn[wpm[0]]-vloumn[wpm[2]])
    dH = -0.25*(hydropathy_index[wpm[0]]-hydropathy_index[wpm[2]])
    return [dP, dF1, dB, dST, dHT, dV, dH]
def getStatical(msa,v_out,w_out,seq,mutation):
    wild = mutation[0]
    pos = int(mutation[1])
    mut = mutation[2]
    wild_score = []
    mut_score = []
    i = 0
    wild_v = v_out[pos-1][a2n[wild]]
    mut_v = v_out[pos-1][a2n[mut]]
    for aa in seq:
        wild_score.append(w_out[i][pos-1][:-1].reshape(21,21)[a2n[aa]][a2n[wild]])
        mut_score.append(w_out[i][pos-1][:-1].reshape(21,21)[a2n[aa]][a2n[mut]])
        i += 1
    
    return [-0.1*(wild_v-mut_v),
            -0.5*(msa[a2n[wild]][pos-1]-msa[a2n[mut]][pos-1]),
            -0.5*(max(wild_score) - max(mut_score)),
            0.5*(min(wild_score) - min(mut_score)),
            -0.15*(sum(wild_score) - sum(mut_score)),
            -0.1*(np.linalg.norm((np.array(wild_score)), ord=1)-np.linalg.norm((np.array(mut_score)), ord=1)),
            -0.5*(np.linalg.norm((np.array(wild_score)), ord=2)-np.linalg.norm((np.array(mut_score)), ord=2))]


if __name__ == '__main__':
    feature = feature()
    al_dd = feature.get_al_dd()
    seq = feature.getSeq(fasta)
    msa = feature.get_MSA(seqa3m)
    static_features = feature.getStatic(wpm)
    statical_features = feature.getStatical(msa,v_out,w_out,seq,mutation)
    





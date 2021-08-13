#!/usr/bin/env python3
from requests import get, post
from time import sleep
import os


def read_fasta(seqfilename):
    seq = ''
    with open(seqfilename) as fasta:
        for line in fasta:
            if line[0] == ">":
                continue
            else:
                seq += line.strip()
    return seq

# submit a new job


def run_mmseqs(seqfilename):

    repeat = True
    while repeat:
        ID = os.popen("curl -s -F q=@"+seqfilename+" -F mode=env https://a3m.mmseqs.com/ticket/msa | jq -r '.id'").read()
        status = os.popen("curl -s https://a3m.mmseqs.com/ticket/"+ID+" | jq -r '.status'").read()
        #status = get('https://search.mmseqs.com/api/ticket/' + ticket['id']).json()
        if status['status'] == "ERROR":
            # handle error
            #sys.exit(0)
            print("MMseq2 error!")
            return 0

        # wait a short time between poll requests
        sleep(1)
        repeat = status['status'] != "COMPLETE"


    print("MMseq2 finished!")
    os.popen("curl -s https://a3m.mmseqs.com/result/download/"+ID+" > results.mmseqs2.tar.gz")
    with open("MMseq2_ready","w+") as mmseq2:
        mmseq2.write("ready")
        mmseq2.close()


    os.popen("tar xzf results.mmseqs2.tar.gz")
    os.popen("cat uniref.a3m bfd.mgnify30.metaeuk30.smag30.a3m > tmp.a3m")
    os.popen("tr -d '\000' < tmp.a3m > "+seqfilename+".a3m")
    with open("A3M_ready","w") as a3m:
        a3m.write("ready.")
        a3m.close()
    return seqfilename+".a3m"

if __name__ == '__main__':
    import sys

    seqfilename = sys.argv[1]
    run_mmseqs(seqfilename)



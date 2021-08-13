#!/usr/bin/python

#By Sun Jinyuan, 2021

# use HHblist to search sequences and build a3m file

from os import popen
import subprocess
import time
#subprocess.call('a.exe')


def hhsearch(seqfilename, iter_num, path_to_database, num_threads):
    searchcmd = "hhblits -i " + seqfilename + " -o " + seqfilename + ".hhr -oa3m " + seqfilename + ".a3m -n " + str(
        iter_num) + " -d " + path_to_database + " -cpu " + str(num_threads)
    # print("grep \">\" " + seqfilename + ".a3m|wc -l")
    search = subprocess.Popen(searchcmd,shell=True)
    while search.poll() != 0:
        time.sleep(1)

    #if search.poll() == 0:
    hits_num = popen("grep \">\" " + seqfilename + ".a3m|wc -l").read()
    print("Found " + hits_num + "hits!")
    a3mfilename = seqfilename + ".a3m"
    return a3mfilename
    #else:
        #search.wait(1)


if __name__ == '__main__':
    # print_hi('PyCharm')
    seqfilename = "g.fasta"
    iter_num = "3"
    path_to_database = "/ydata/jsun/database/UniRef30_2020_03"
    num_threads = 8

    hhsearch(seqfilename, iter_num, path_to_database, num_threads)

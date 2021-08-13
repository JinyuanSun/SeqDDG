#!/usr/bin/python

#By Sun Jinyuan, 2021

import argparse
def get_parser():
    parser = argparse.ArgumentParser(description='sequence based ∆∆G prediction')
    parser.add_argument("-s",
                        "--sequence",
                        help="target sequence which must be provided")
    parser.add_argument("-a",
                        "--a3mfile",
                        help="the sequence alignment file in a3m format")
    parser.add_argument("-db",
                        "--database",
                        help="path/to/UniRep30_2020_03 for hhblits")
    parser.add_argument("-ni",
                        '--iteration_number',
                        default=3,
                        help='the iteration number of hhblits search')
    parser.add_argument("-nt",
                        '--num_threads',
                        default=4,
                        help='number of threads used to run hhblist, default is 4')
    parser.add_argument("-m",
                        '--mutation',
                        help='mutation to predict, in form of wild_resnum_mut, like A_24_G')
    parser.add_argument("-ml",
                        '--mutation_list',
                        help='list of mutation to predict, one mutation per line')
            
    return parser

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    seqfilename = args.sequence
    a3mfilename = args.a3mfile
    path_to_database = args.database
    iter_num = args.iteration_number
    num_threads = args.num_threads
    mutation = args.mutation
    mutation_list = args.mutation_list
    print(seqfilename,a3mfilename)


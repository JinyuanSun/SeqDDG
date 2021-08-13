

# SeqDDG Documentation

## Introduction

SeqDDG is a sequence-based protein single-point mutation effect predictor. SeqDDG using multiple sequence alignemnt (MSA) together with Potts model to capture the co-variation of residue pairs. SeqDDG achieved state-of-the-art performance on S830 among all tested sequence based predictor (see the Table below ). SeqDG also showed comparable precision on prediction of ∆Tm (see the Figure below) . Detailed methods and results please refer to the paper.

| Algorithm     | Pearson’s r |
| ------------- | :---------- |
| **SeqDDG***   | **0.56**    |
| PoPMuSiC      | 0.56        |
| FoldX         | 0.51        |
| EASE-MM*      | 0.45        |
| Rosetta_NoMin | 0.33        |
| BoostDDG*     | N. A.       |
| DynaMut       | N. A.       |

*Sequence-based predictor. N. A. for server not returned results.

![image](http://185.201.226.155/Precisions.png)

## Contact

We appreciate bug reports, comments, and suggestions.  

> Email: jysun@im.ac.cn

For technical details please also posts issue in this repo.

## Files and Directories 

Organizations of the seqddg package:

```tex
seqddg/
    __init__.py
    scanner.py
    single_point_predictor.py
    utilities/
        __init__.py
        HHsearch.py
        feature_maker.py
        parsermodule.py
        potts.py
    model/
        SeqDDG.h5
    demo/
        Gb1.fasta
        Gb1.fasta.a3m
```

## Dependency 

SeqDDG is wirtten in python3 and tested under Ubuntu 20.04.

[TensorFlow 2.x](https://www.tensorflow.org/) to run the NN,

[HHblits and UniRef30_2020_03](https://github.com/soedinglab/hh-suite) database for sequence search develpoed by soedinglab.

## Instructions

To run a demo on Gb1 domain:

```bash
python3 scanner.py -s demo/Gb1.fasta -a demo/Gb1.fasta.a3m
```

For detailed help information:

```tex
usage: scanner.py [-h] [-s SEQUENCE] [-a A3MFILE] [-db DATABASE]
                  [-ni ITERATION_NUMBER] [-nt NUM_THREADS] [-m MUTATION]
                  [-ml MUTATION_LIST]

sequence based ∆∆G prediction

optional arguments:
  -h, --help            show this help message and exit
  -s SEQUENCE, --sequence SEQUENCE
                        target sequence which must be provided
  -a A3MFILE, --a3mfile A3MFILE
                        the sequence alignment file in a3m format
  -db DATABASE, --database DATABASE
                        path/to/UniRep30_2020_03 for hhblits
  -ni ITERATION_NUMBER, --iteration_number ITERATION_NUMBER
                        the iteration number of hhblist search
  -nt NUM_THREADS, --num_threads NUM_THREADS
                        number of threads used to run hhblist, default is 4
  -m MUTATION, --mutation MUTATION
                        mutation to predict, in form of wild_resnum_mut, like
                        A_24_G
  -ml MUTATION_LIST, --mutation_list MUTATION_LIST
                        list of mutation to predict, one mutation per line
```

The `scanner.py` and `single_point_predictor.py` now sharing the same parse moduel, for `scanner.py`,  `-m` and `-ml` will be ignored.

An example of `mutation_list` file:

```tex
A_24_G
A_24_M
```

## Reference

Sun J.,Cui YL, Wu B.* SeqDDG: Sequence-based protein point mutation effect predictor. *Unpublished*. 2021.


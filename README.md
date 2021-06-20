

# SeqDDG Documentation

## Introduction

SeqDDG is a sequence-based protein single-point mutation effect predictor. SeqDDG using multiple sequence alignemnt (MSA) together with Potts model to capture the co-variation of residue pairs. SeqDDG achieved state-of-the-art performance on S830 among all tested sequence based predictor. Detailed methods and results please refer to the paper.

## Reference

Sun J., Wu B.* SeqDDG: Sequence-based protein point mutation effect predictor. *Unpublished*. 2021.

## Contact

We appreciate bug reports, comments, and suggestions.  

> Email: jysun@im.ac.cn

For technical details please post issue in this repo.

## Files and Directories 

Organizations of the seqddg package:

```bash
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
        demo.py
```

## Dependency 

**TensorFlow >= 2.0, HHblits-3 and UniRef30_2020_03 database** 

## Instructions

To run a 


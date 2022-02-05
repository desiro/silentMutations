# <samp>SilentMutations</samp> - manual

***

## Usage
```
silentMutations.py -f <in_fasta> -p <out_prefix> -s1 <name>:<frame>:<start>-<end> -s2 <name>:<frame>:<start>-<end> [options]
```

## Version
```
silentMutations.py 1.0.0
```

## Dependencies
```Python v3.9.7```, ```NumPy v1.22.2```, ```ViennaRNA v2.4.18```, ```VARNA v3.93```

## Description
<samp>SilentMutations</samp> generates compensatory codon mutations with similar minimum free energy. It is recommended to not extend the length for each sequence over 30 nt. This could lead to extensive run times and each sequence would probably form intra sequence interactions prior to their interaction with each other. This is not covered with SIM. Example call: python3 silentMutations.py -p example -f example.fa -s1 seq1:2:20-36 -s2 seq2:0:23-44 -cls ssRNA- -r -c -thr=4

## Options

```
--prefix,-p
    output prefix for result files

--fasta,-f
    fasta file with two sequences; will create all possible silent codon combinations for a given range of the two segments

--snip1,-s1
    define the range of the first snip with <name>:<frame>:<start>-<end>

--snip2,-s2
    define the range of the second snip with <name>:<frame>:<start>-<end>

--complement,-c
    creates complements of each strain if set (default: False)

--reverse,-r
    creates reverse of each strain if set (default: False)

--virusclass,-cls
    define virus class; (default: ssRNA+) (choices: ssRNA+,ssRNA-)

--filterperc,-prc
    define minimum percentage for the mutant-WT mfe to be below the WT-WT mfe
    (default: 0.5)

--filterperx,-prx
    optionally define a maximum percentage for the mutant-WT mfe to be above the
    WT-WT mfe; the only purpose of this parameter is to reduce the runtime; it 
    should be used with care (default: 0.0)

--upperdev,-udv
    define upper similarity threshold for the mutant-mutant mfe in contrast to
    the WT-WT mfe (default: 1.05)

--lowerdev,-ldv
    define lower similarity threshold for the mutant-mutant mfe in contrast to
    the WT-WT mfe (default: 0.95)

--mutrange,-mrg
    define percentage of similarity for the mean mfe of both single mutants
    (default: 0.1)

--noncoding1,-nc1
    specify if the first sequence includes a non-coding region (default: False)

--noncoding2,-nc2
    specify if the second sequence includes a non-coding region (default: False)

--mutations,-mut
    define maximum number of allowed mutations per sequence, lower values will
    improve runtime (default: 5)

--ignoreMutations,-imt
    ignores all mutations that don't have at least one sequence with exactly the
    maximum number allowed mutations, useful to reduce runtime with higher mutation
    counts (default: False)

--prioMut1,-pm1
    prioritize to minimize first mutant sequence (default: False)

--prioMut2,-pm2
    prioritize to minimize second mutant sequence (default: False)

--stabilize,-stb
    stabilze interaction instead of destabilization (default: False)

--dangles,-dng
    use dangling ends for foldings (default: 0) (choices: 0,1,2,3)

--noLP,-nlp
    disable lonely pairs for RNAcofold (default: False)

--threads,-thr
    number of threads to use for RNAcofold (default: 1)

--varnabin,-var
    use this VARNA binary; example: VARNAv3-93.jar (default: )

--inkscbin,-ink
    use this Inkscape binary; example: inkscape (default: )

reference
    https://doi.org/10.1016/j.virusres.2018.11.005
```

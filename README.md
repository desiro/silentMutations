# silentMutations

This tool can automatically construct interrupting and restoring silent mutation pairs within coding regions for combinatorial *in-silico* analysis of RNA-RNA interactions. The tool has been written in ```Python 3.6.5``` and relies heavily on the ```RNAcofold``` python site-package of the ```ViennaRNA Package 2.4```.

## Mandatory Prerequisites

* [python 3.6](https://www.python.org/downloads/release/python-365/)
* [numpy](http://www.numpy.org/)
* [viennaRNA 2.4](https://www.tbi.univie.ac.at/RNA/documentation.html#install)

### Conda Installation

Installing everything with Conda:
```
conda install -c bioconda  viennarna=2.4.13
conda install -c conda-forge numpy=1.16.4
conda install -c lb_arrakistx varna=3.93
git clone https://github.com/desiro/silentMutations.git 
```

### Docker Installation

Using Docker Hub:
```
docker run --user $(id -u):$(id -g) --rm -v <your working directory>:/source --workdir /source desiro/sim silentMutations.py -p example -f example.fa -s1 seq1:2:20-36 -s2 seq2:0:23-44 -cls ssRNA- -r -c -thr=4
```
* It is currently not possible to use VARNA within the Docker Hub repository due to some java issues.

### Unix Installation

Installing the ViennaRNA package on Linux:
```
tar -zxvf ViennaRNA-2.4.9.tar.gz
cd ViennaRNA-2.4.9
./configure --with-python3
make
sudo make install
```

Installing the ViennaRNA package on MAC:
```
tar -zxvf ViennaRNA-2.4.9.tar.gz
cd ViennaRNA-2.4.9
./configure --enable-universal-binary --with-python3
make
sudo make install
```

For Windows 10 users, please use the [Ubuntu](https://www.microsoft.com/en-us/p/ubuntu/9nblggh4msv6?cid=msft_web_chart) subsystem and follow the Linux installation steps. 

## Optional Prerequisites

* [VARNA 3.93](http://varna.lri.fr/)
* [Inkscape 0.92](https://inkscape.org/en/)

## Examples

The basic input for ```silentMutations.py``` includes the following parameters:
* ```-p``` - name of the output folder
* ```-f``` - input fasta file, will only read the first two sequences
* ```-s1 <name>:<frame>:<start>-<end>``` - for the first sequence
* ```-s2 <name>:<frame>:<start>-<end>``` - for the second sequence
  
### Basic Example
```
python3 silentMutations.py -p example -f example.fa -s1 seq1:2:20-36 -s2 seq2:0:23-44 -cls ssRNA- -r -c -thr=4
```

### Example with VARNA and Inkscape
```
python3 silentMutations.py -p example -f example.fa -s1 seq1:2:20-36 -s2 seq2:0:23-44 -cls ssRNA- -r -c -thr=4 -var VARNAv3-93.jar -ink inkscape
```

## Authors

* [Daniel Desirò](https://github.com/desiro)

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Reference

Please cite [SilentMutations](https://www.sciencedirect.com/science/article/pii/S016817021830577X) if you find our tool useful.
There is also an open access version at [bioRxiv](https://doi.org/10.1101/424002).

## Workflow overview

![workflow](https://github.com/desiro/silentMutations/blob/master/workflow.png "(a) extract sequences and remove unpaired codons (b) create possible codon permutations (c) keep only mutants with a weak mutant-WT fold mfe (d) keep only double-mutants with similar double-WT fold mfe (e) minimize fold mfe of single-mutants")

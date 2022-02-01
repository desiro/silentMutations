
## <samp>SilentMutations</samp>
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-teal.svg)](https://www.gnu.org/licenses/gpl-3.0)![Python3.6](https://img.shields.io/badge/Language-Python_3.6-darkred.svg)

***

## Description

This tool can automatically construct interrupting and restoring silent mutation pairs within coding regions for combinatorial *in silico* analysis of RNA-RNA interactions. The predictions can be used for *in vitro* and *in vivo* experiments.

***

## Mandatory Prerequisites

* [python 3.6](https://www.python.org/downloads/release/python-365/)
* [numpy](http://www.numpy.org/)
* [viennaRNA 2.4](https://www.tbi.univie.ac.at/RNA/documentation.html#install)

## Optional Prerequisites

* [VARNA 3.93](http://varna.lri.fr/)
* [Inkscape 0.92](https://inkscape.org/en/)

***

### Conda Installation

Installing everything with Conda:
```
conda install -c bioconda  viennarna=2.4.13
conda install -c conda-forge numpy=1.16.4
conda install -c lb_arrakistx varna=3.93
git clone https://github.com/desiro/silentMutations.git 
```

### Other Installation

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

***

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

***

## Authors

* [Daniel Desirò](https://github.com/desiro)

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Reference

Please cite [SilentMutations](https://www.sciencedirect.com/science/article/pii/S016817021830577X) if you find our tool useful.
There is also an open access version at [bioRxiv](https://doi.org/10.1101/424002).

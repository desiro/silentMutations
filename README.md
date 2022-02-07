
# <samp>SilentMutations</samp>
[![License: GPL v3](https://img.shields.io/badge/License-GPL_v3-bd0000.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python v3.9.7](https://img.shields.io/badge/Language-Python_v3-75a8d3.svg)](https://www.python.org/)
[![Conda v4.11.0](https://img.shields.io/badge/Uses-Conda-43b02a.svg)](https://docs.conda.io/en/latest/miniconda.html)

***

## Description

This tool can automatically construct interrupting and restoring silent mutation pairs within coding regions for combinatorial *in silico* analysis of RNA-RNA interactions. The predictions can be used for *in vitro* and *in vivo* experiments.

***

## Mandatory Prerequisites

* [![Python v3.9.7](https://img.shields.io/badge/Python_v3.9.7-75a8d3.svg)](https://www.python.org/downloads/release/python-397/)
* [![NumPy v1.22.2](https://img.shields.io/badge/NumPy_v1.22.2-013243.svg)](http://www.numpy.org/)
* [![ViennaRNA v2.5.0](https://img.shields.io/badge/ViennaRNA_v2.5.0-006795.svg)](https://www.tbi.univie.ac.at/RNA/)

## Optional Prerequisites

* [![Conda v4.11.0](https://img.shields.io/badge/Conda_v4.11.0-43b02a.svg)](https://docs.conda.io/en/latest/miniconda.html)
* [![VARNA v3.93](https://img.shields.io/badge/VARNA_v3.93-ffba27.svg)](http://varna.lri.fr/)

***

### Conda Installation

To run <samp>SilentMutations</samp>, I recommend using Miniconda and following the steps below.

Installing Miniconda:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
conda update conda
conda update python
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Installing Conda packages:
```
conda create --name silentMutations python=3.9.7
conda activate silentMutations
conda install -c bioconda viennarna=2.5.0
conda install -c conda-forge numpy=1.22.2
conda install -c lb_arrakistx varna=3.93
git clone https://github.com/desiro/silentMutations.git
cd silentMutations
```

### Other Installations

Installing the ViennaRNA package on Linux:
```
tar -zxvf ViennaRNA-2.5.0.tar.gz
cd ViennaRNA-2.5.0
./configure --with-python3
make
sudo make install
```

Installing the ViennaRNA package on MAC:
```
tar -zxvf ViennaRNA-2.5.0.tar.gz
cd ViennaRNA-2.5.0
./configure --enable-universal-binary --with-python3
make
sudo make install
```

For Windows 10 users, please use the [Ubuntu 20.04 LTS](https://www.microsoft.com/en-us/p/ubuntu-2004-lts/9n6svws3rx71?cid=msft_web_chart) subsystem and follow the Linux installation steps. 

***

## Examples

The basic input for ```silentMutations.py``` includes the following parameters:
* ```-p``` - name of the output folder
* ```-f``` - input fasta file, will only read the first two sequences
* ```-s1 <name>:<frame>:<start>-<end>``` - for the first sequence
* ```-s2 <name>:<frame>:<start>-<end>``` - for the second sequence
  
### Basic Example
```
python silentMutations.py -p example -f example.fa -s1 seq1:2:20-36 -s2 seq2:0:23-44 -cls ssRNA- -r -c -thr=4
```

### Example with VARNA
```
python silentMutations.py -p example -f example.fa -s1 seq1:2:20-36 -s2 seq2:0:23-44 -cls ssRNA- -r -c -thr=4 -var VARNAv3-93.jar
```

### Options

For more command line options, see the [manual](https://github.com/desiro/silentMutations/blob/master/manual.md).

***

## Authors

* [Daniel Desirò](https://github.com/desiro)

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Reference

Please cite [<samp>SilentMutations</samp>](https://www.sciencedirect.com/science/article/pii/S016817021830577X) if you find our tool useful:
```
D. Desirò, M. Hölzer, B. Ibrahim and M. Marz.
"SilentMutations (SIM): A tool for analyzing long-range RNA–RNA interactions in viral genomes and structured RNAs."
Virus Research, 260:135-141, 2019.
```

There is also an open access version at [bioRxiv](https://doi.org/10.1101/424002).

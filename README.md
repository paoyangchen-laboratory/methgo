# MethGo
**MethGo** is a comprehensive tool for analyzing whole genome bisulfite sequencing data.

## Dependencies
Make sure you've installed [Python 2.7](https://www.python.org/downloads/release/python-279/),
[samtools](http://sourceforge.net/projects/samtools/files/samtools/1.2/),
and [bedtools](http://bedtools.readthedocs.org/en/latest/content/installation.html)
before you install packages listed below:

- [NumPy](http://www.numpy.org/)
- [matplotlib](http://matplotlib.org/)
- [pandas](http://pandas.pydata.org/)
- [pysam](http://pysam.readthedocs.org/)
- [biopython](http://biopython.org/)
- [pyfasta](https://pypi.python.org/pypi/pyfasta/)
- [Cython](http://cython.org/)
- [pybedtools](https://pythonhosted.org/pybedtools/)

## Installation

        $ mkdir ~/.local
        $ cd ~/.local
        $ git clone https://github.com/wwliao/methgo.git
        $ chmod +x methgo/methgo

        # Add "export PATH=$PATH:~/.local/methgo" to ~/.bashrc
        $ source ~/.bashrc

## Setting up environment
0. Downloading the test data

        $ cd ~
        $ curl -O http://paoyang.ipmb.sinica.edu.tw/~gattaca/methgo_demo.tar.gz
        $ tar xvfz methgo_demo.tar.gz
        $ cd methgo_demo

1. Installing virtualenv

        $ curl -O https://pypi.python.org/packages/source/v/virtualenv/virtualenv-12.1.1.tar.gz
        $ tar xvfz virtualenv-12.1.1.tar.gz
        $ python virtualenv-12.1.1/virtualenv.py --no-site-packages --python=python2.7 .venv

2. Activating the virtual environment

        $ source .venv/bin/activate

3. Installing packages

        $ pip install -r ~/.local/methgo/requirements/base.txt
        $ pip install -r ~/.local/methgo/requirements/addition.txt

## Demo how to use it
0. Change to data folder

        $ cd data

1. Running MET module

        $ methgo met genes.gtf genome.fa demo.CGmap

2. Running TXN module

        $ methgo txn -t ~/.local/methgo/scripts/txn/tair10_txn -l ATHB1_binding_site_motif,CCA1_binding_site_motif -c demo.CGmap

3. Running CNV module

        $ methgo cnv genome.fa.fai demo.sorted.bam

4. Running SNP module

        $ methgo snp -g genome.fa demo.sorted.bam

## License
MIT

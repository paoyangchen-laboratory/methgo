# MethGo
**MethGo** is a comprehensive tool for analyzing whole genome bisulfite sequencing data.

## Dependencies
It needs Python 2.7, samtools, bedtools, and the packages listed below:
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
        $ git clone git@github.com:wwliao/methgo.git
        $ chmod +x methgo/methgo

- Add this line to ~/.bashrc file

        export PATH=$PATH:~/.local/methgo

## Setting up environment
1. Installing virtualenv

        $ curl -O https://pypi.python.org/packages/source/v/virtualenv/virtualenv-12.1.1.tar.gz
        $ tar xvfz virtualenv-12.1.1.tar.gz
        $ python virtualenv-12.1.1/virtualenv.py --no-site-packages --python=python2.7 .venv

2. Activating the virtual environment

        $ source .venv/bin/activate

3. Installing packages

        $ pip install -r methgo/requirements.txt

## Demo how to use it
1. Downloading the test data
2. Running SNP module

        $ methgo snp -g genome.fa demo.sorted.bam

3. Running CNV module

        $ methgo cnv genome.fa.fai demo.sorted.bam

4. Running MLEVEL module

        $ methgo mlevel genes.gtf genome.fa demo.CGmap

5. Running TXN module

        $ methgo txn -t ~/.local/methgo/scripts/txn/tair10_txn -l CCA1_binding_site_motif -c demo.CGmap

## License
MIT

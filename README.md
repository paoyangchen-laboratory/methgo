# MethGo: a comprehensive tool for analyzing whole genome bisulfite sequencing data

## Authors
Ming-Ren Yen, Wen-Wei Liao, Evaline Ju, Fei-Man Hsu, Larry Lam, Pao-Yang Chen

## Affiliation
Institute of Plant and Microbial Biology, Academia Sinica, Taipei 11529, Taiwan

## Dependencies
Please install [Python 2.7](https://www.python.org/downloads/release/python-279/),
[samtools](http://sourceforge.net/projects/samtools/files/samtools/1.2/),
and [bedtools](http://bedtools.readthedocs.org/en/latest/content/installation.html).

## Install MethGo

        $ mkdir ~/.local
        $ cd ~/.local
        $ git clone https://github.com/wwliao/methgo.git
        $ chmod +x methgo/methgo
        $ echo 'export PATH=$PATH:~/.local/methgo' >> ~/.bashrc
        $ source ~/.bashrc

## Set up environment with a sample input file
0. Download the a sample input file

        $ cd ~
        $ curl -O http://paoyang.ipmb.sinica.edu.tw/~gattaca/methgo_demo.tar.gz
        $ tar xvfz methgo_demo.tar.gz
        $ cd methgo_demo

1. Install and activate a virtual environment (10 mins)

        $ curl -O https://pypi.python.org/packages/source/v/virtualenv/virtualenv-12.1.1.tar.gz
        $ tar xvfz virtualenv-12.1.1.tar.gz
        $ python virtualenv-12.1.1/virtualenv.py --no-site-packages --python=python2.7 .venv
        $ source .venv/bin/activate
        $ pip install -r ~/.local/methgo/requirements/base.txt
        $ pip install -r ~/.local/methgo/requirements/addition.txt
        
	MethGo requires the following packages:
    - [NumPy](http://www.numpy.org/)
	- [matplotlib](http://matplotlib.org/)
	- [pandas](http://pandas.pydata.org/)
	- [pysam](http://pysam.readthedocs.org/) == 0.8.0
	- [biopython](http://biopython.org/)
	- [pyfasta](https://pypi.python.org/pypi/pyfasta/)
	- [Cython](http://cython.org/)
	- [pybedtools](https://pythonhosted.org/pybedtools/)

## Tutorial
0. Change to data folder

        $ cd ~/methgo_demo/data

1. Run MET module (5 mins)

        $ methgo met genes.gtf genome.fa demo.CGmap
        
    Output figures:
    - Global methylation levels ([demo.bulk.mean.png](http://paoyang.ipmb.sinica.edu.tw/~gattaca/methgo_demo_results/MET/demo.bulk.mean.png))
    - Methylation level distribution ([demo.bulk.hist.png](http://paoyang.ipmb.sinica.edu.tw/~gattaca/methgo_demo_results/MET/demo.bulk.hist.png))
    - Methylation levels in different genomic regions ([demo.feature.CG.png](http://paoyang.ipmb.sinica.edu.tw/~gattaca/methgo_demo_results/MET/demo.feature.CG.png),
    [demo.feature.CHG.png](http://paoyang.ipmb.sinica.edu.tw/~gattaca/methgo_demo_results/MET/demo.feature.CHG.png),
    [demo.feature.CHH.png](http://paoyang.ipmb.sinica.edu.tw/~gattaca/methgo_demo_results/MET/demo.feature.CHH.png))
    - Global methylation level of all chromosomes ([demo.genomewide.png](http://paoyang.ipmb.sinica.edu.tw/~gattaca/methgo_demo_results/MET/demo.genomewide.png))

    Output files:
    - Methylation levels of each gene ([demo.feature.CG.txt](http://paoyang.ipmb.sinica.edu.tw/~gattaca/methgo_demo_results/MET/demo.feature.CG.txt),
    [demo.feature.CHG.txt](http://paoyang.ipmb.sinica.edu.tw/~gattaca/methgo_demo_results/MET/demo.feature.CHG.txt),
    [demo.feature.CHH.txt](http://paoyang.ipmb.sinica.edu.tw/~gattaca/methgo_demo_results/MET/demo.feature.CHH.txt))
   
2. Run TXN module (5 mins)

        $ methgo txn -t ~/.local/methgo/scripts/txn/tair10_txn -l ATHB1_binding_site_motif,CCA1_binding_site_motif -c demo.CGmap
   
    Output figures:
    - Methylation level around transcription binding site ([demo.CG.txn.png](http://paoyang.ipmb.sinica.edu.tw/~gattaca/methgo_demo_results/TXN/demo.CG.txn.png), [demo.CHG.txn.png](http://paoyang.ipmb.sinica.edu.tw/~gattaca/methgo_demo_results/TXN/demo.CHG.txn.png), [demo.CG.txn.png](http://paoyang.ipmb.sinica.edu.tw/~gattaca/methgo_demo_results/TXN/demo.CHH.txn.png))

3. Run CNV module (5 mins)

        $ methgo cnv genome.fa.fai demo.sorted.bam

    Output figure:
    - CNV on all chromosomes ([demo.sorted.cov.png](http://paoyang.ipmb.sinica.edu.tw/~gattaca/methgo_demo_results/CNV/demo.sorted.cov.png))

4. Run SNP module (5 mins)

        $ methgo snp -g genome.fa demo.sorted.bam

    Output files:
    - homozygous SNPs ([demo.sorted.homozygous.txt](http://paoyang.ipmb.sinica.edu.tw/~gattaca/methgo_demo_results/SNP/demo.sorted.homozygous.txt))
    - heterozygous SNPs ([demo.sorted.heterozygous.txt](http://paoyang.ipmb.sinica.edu.tw/~gattaca/methgo_demo_results/SNP/demo.sorted.heterozygous.txt))

## License
MIT

## Acknowledgements This work was supported by a grant from Academia Sinica, and grants from  MOST-103-2313-B-001-003-MY3 and MOST-103-2633-B-001-002 and NHRI-EX104-10324SC.

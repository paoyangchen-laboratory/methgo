Tutorial
========

1. Download the sample input file

  ::

  $ curl -O http://paoyang.ipmb.sinica.edu.tw/~gattaca/methgo_demo.tar.gz
  $ tar xvfz methgo_demo.tar.gz
  $ cd methgo_demo/data

2. Run COV module:

  ::

  $ methgo cov genome.fa demo.CGmap

  .. image:: _static/demo.cov.png

3. Run MET module:

  ::

  $ methgo met genes.gtf genome.fa demo.CGmap

4. Run TXN module:

  ::

  $ methgo txn -t methgo/scripts/txn/tair10_txn -l ATHB1_binding_site_motif,CCA1_binding_site_motif -c demo.CGmap

5. Run SNP module:

  ::

  $ methgo snp -g genome.fa demo.sorted.bam

6. Run CNV module:

  ::

  $ methgo cnv genome.fa.fai demo.sorted.bam

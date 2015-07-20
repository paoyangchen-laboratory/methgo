#!/usr/bin/env python
from __future__ import division
import re
import os
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt

def get_parser():
    """
    Create a parser and add arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='reference genome FASTA file')
    parser.add_argument('cgmap', help='CGmap file')
    return parser

def get_ctxnum(reffile):
    """
    Get the number of CG/CHG/CHH from a reference genome FASTA file
    """
    with open(reffile) as infile:
        fasta = SeqIO.to_dict(SeqIO.parse(infile, 'fasta'))
        for chr in fasta:
            fasta[chr] = str(fasta[chr].seq).upper()
    num_cg = 0
    num_chg = 0
    num_chh = 0
    for chr in fasta:
        num_cg += len([match.start() for match in re.finditer(r'(?=(CG))', fasta[chr])])
        num_cg += len([match.start()-1 for match in re.finditer(r'(?<=(CG))', fasta[chr])])
        num_chg += len([match.start() for match in re.finditer(r'(?=(C[ACT]G))', fasta[chr])])
        num_chg += len([match.start()-1 for match in re.finditer(r'(?<=(C[AGT]G))', fasta[chr])])
        num_chh += len([match.start() for match in re.finditer(r'(?=(C[ACT][ACT]))', fasta[chr])])
        num_chh += len([match.start()-1 for match in re.finditer(r'(?<=([AGT][AGT]G))', fasta[chr])])
    return num_cg, num_chg, num_chh

def read_ctxcov(cgmapfile):
    """
    Read the column of coverage from CGmap file
    """
    cgmap = pd.read_table(cgmapfile, header=None, usecols=[3, 7], names=['ctx', 'cov'])
    cov_cg = cgmap[cgmap['ctx'] == 'CG']['cov']
    cov_chg = cgmap[cgmap['ctx'] == 'CHG']['cov']
    cov_chh = cgmap[cgmap['ctx'] == 'CHH']['cov']
    return cov_cg, cov_chg, cov_chh

def plot_ctxcov(num_cg, num_chg, num_chh, cov_cg, cov_chg, cov_chh):
    """
    Plot the CG/CHG/CHH coverage distribution
    """
    colors = { 'CG': (38/255, 173/255, 84/255),
              'CHG': (44/255, 180/255, 234/255),
              'CHH': (249/255, 42/255, 54/255)}

    plt.switch_backend('Agg')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_position(('outward', 8))
    ax.spines['left'].set_position(('outward', 8))
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    maxval = int(max([cov_cg.mean(), cov_chg.mean(), cov_chh.mean()]) + 2*max([cov_cg.std(), cov_chg.std(), cov_chh.std()]))
    n, bins = np.histogram(cov_cg, bins=np.linspace(0, maxval, maxval+1))
    n = np.cumsum(n[::-1])[::-1]
    ax.plot(np.arange(0.5, maxval, 1), n/num_cg*100, linewidth=2, color=colors['CG'], label='CG')
    n, bins = np.histogram(cov_chg, bins=np.linspace(0, maxval, maxval+1))
    n = np.cumsum(n[::-1])[::-1]
    ax.plot(np.arange(0.5, maxval, 1), n/num_chg*100, linewidth=2, color=colors['CHG'], label='CHG')
    n, bins = np.histogram(cov_chh, bins=np.linspace(0, maxval, maxval+1))
    n = np.cumsum(n[::-1])[::-1]
    ax.plot(np.arange(0.5, maxval, 1), n/num_chh*100, linewidth=2, color=colors['CHH'], label='CHH')
    ax.set_xlim(0, maxval)
    ax.set_ylim(0, 100)
    ax.tick_params(direction='out', top='off', right='off', length=5, width=2, labelsize='large')
    ax.set_xlabel('Coverage (x)', size='large', weight='bold')
    ax.set_ylabel('Percentage (%)', size='large', weight='bold')
    ax.legend(loc='upper right', prop={'size': 'small'})
    plt.tight_layout()
    return ax

def main():
    parser = get_parser()
    args = parser.parse_args()
    num_cg, num_chg, num_chh = get_ctxnum(args.fasta)
    cov_cg, cov_chg, cov_chh = read_ctxcov(args.cgmap)
    ctxcov_ax = plot_ctxcov(num_cg, num_chg, num_chh, cov_cg, cov_chg, cov_chh)
    fig = ctxcov_ax.get_figure()
    prefix = os.path.splitext(os.path.basename(args.cgmap))[0]
    fig.savefig('{}.cov.png'.format(prefix), dpi=300)
    plt.close(fig)

if __name__ == '__main__':
    main()

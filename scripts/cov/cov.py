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

def main():
    parser = get_parser()
    args = parser.parse_args()
    num_cg, num_chg, num_chh = get_ctxnum(args.fasta)
    cgmap = pd.read_table(args.cgmap, header=None, usecols=[3, 7], names=['ctx', 'cov'])
    cg = cgmap[cgmap['ctx'] == 'CG']['cov']
    chg = cgmap[cgmap['ctx'] == 'CHG']['cov']
    chh = cgmap[cgmap['ctx'] == 'CHH']['cov']
    plt.switch_backend('Agg')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_position(('outward', 8))
    ax.spines['left'].set_position(('outward', 8))
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    #n, bins, patches = ax.hist((cg, chg, chh), bins=np.linspace(0, 40, 9), normed=True, cumulative=1, color=[(38/255, 173/255, 84/255), (44/255, 180/255, 234/255), (249/255, 42/255, 54/255)], edgecolor='w', label=['CG', 'CHG', 'CHH'])
    n, bins = np.histogram(cg, bins=np.linspace(0, 40, 41))
    n = np.cumsum(n[::-1])[::-1]
    ax.plot(np.arange(0.5, 40, 1), n/num_cg*100, linewidth=2, color=(38/255, 173/255, 84/255), label='CG')
    n, bins = np.histogram(chg, bins=np.linspace(0, 40, 41))
    n = np.cumsum(n[::-1])[::-1]
    ax.plot(np.arange(0.5, 40, 1), n/num_chg*100, linewidth=2, color=(44/255, 180/255, 234/255), label='CHG')
    n, bins = np.histogram(chh, bins=np.linspace(0, 40, 41))
    n = np.cumsum(n[::-1])[::-1]
    ax.plot(np.arange(0.5, 40, 1), n/num_chh*100, linewidth=2, color=(249/255, 42/255, 54/255), label='CHH')
    ax.set_ylim(0, 100)
    ax.tick_params(direction='out', top='off', right='off', length=5, width=2, labelsize='large')
    ax.set_xlabel('Coverage (x)', size='large', weight='bold')
    ax.set_ylabel('Percentage (%)', size='large', weight='bold')
    ax.legend(loc='upper right', prop={'size': 'small'})
    plt.tight_layout()
    root = os.path.splitext(os.path.basename(args.cgmap))[0]
    plt.savefig('{}.cov.png'.format(root), dpi=300)
    return n, bins

if __name__ == '__main__':
    n, bins = main()

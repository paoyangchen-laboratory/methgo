#!/usr/bin/env python
from __future__ import division
import os
import math
import copy
import random
import argparse
import pysam
import numpy as np
import scipy.stats
import pandas as pd
import matplotlib.pyplot as plt

def extract_chrnum(chr):
    if chr.lower().startswith('chr'):
        try:
            return int(chr[3:])
        except ValueError:
            return chr[3:]
    try:
        return int(chr)
    except ValueError:
        return chr[3:]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-w', '--winsize', type=int, default=200000)
    parser.add_argument('-p', '--pvalue', type=float, default=0.05)
    parser.add_argument('-s', '--succession', type=int, default=3)
    parser.add_argument('refindex')
    parser.add_argument('bamfile')
    args = parser.parse_args()

    chrlen = {}
    with open(args.refindex) as infile:
        for line in infile:
            line = line.strip().split('\t')
            chrlen[line[0]] = int(line[1])

    gcov = {}
    for chr in chrlen:
        gcov[chr] = [float('nan')]*chrlen[chr]

    bam = pysam.Samfile(args.bamfile, 'rb')
    chr = None
    for col in bam.pileup():
        pos = col.pos
        cov = col.n
        chr = bam.getrname(col.tid)
        gcov[chr][pos] = cov

    pos = 0
    chrs = sorted(chrlen.keys(), key=extract_chrnum)
    win_x = []
    win_chr = []
    win_cov = []
    for chr in chrs:
        start = 0
        while (start + args.winsize) <= chrlen[chr]:
            win_x.append(pos+(args.winsize/2))
            win_chr.append((chr, start, start+args.winsize))
            tmp = np.array(gcov[chr][start:start+args.winsize])
            tmp = tmp[~np.isnan(tmp)]
            #win_cov.append(np.mean(tmp))
            win_cov.append(tmp.sum())
            start += args.winsize
            pos += args.winsize

    shuffle_cov = []
    sample_size = 100
    for _ in xrange(sample_size):
        temp = copy.copy(win_cov)
        random.shuffle(temp)
        shuffle_cov.append(temp)

    chr_index = pd.MultiIndex.from_tuples(win_chr, names=['chr', 'start', 'end'])
    shuffle_table = pd.DataFrame(zip(*shuffle_cov), index=chr_index)
    data = pd.DataFrame(win_cov, index=chr_index, columns=['cov'], dtype=int)
    data['mean'] = shuffle_table.mean(axis=1)
    data['std'] = shuffle_table.std(axis=1)
    data['zscore'] = data.apply(lambda s: (s['cov']-s['mean'])/(s['std']/math.sqrt(sample_size)), axis=1)
    data['pval'] = data['zscore'].apply(scipy.stats.norm.sf)
    data = data[['cov', 'pval']]
    data = data[data['pval'] < args.pvalue]
    index = data.index.tolist()
    cnv_cand = []
    temp = []
    prev = index[0]
    for next in index[1:]:
        if prev[0] == next[0]:
            if (next[1] - prev[1]) == args.winsize:
                if len(temp) > 0:
                    if temp[-1] != prev:
                        temp.append(prev)
                else:
                    temp.append(prev)
                temp.append(next)
            else:
                if len(temp) >= args.succession:
                    cnv_cand.extend(temp)
                temp = []
        prev = next
    data = data.ix[pd.MultiIndex.from_tuples(cnv_cand, names=['chr', 'start', 'end'])]
    root = os.path.splitext(os.path.basename(args.bamfile))[0]
    data.to_csv('{}.cnv.txt'.format(root), sep='\t')

    plt.switch_backend('Agg')
    vlines = [0]
    for i, chr in enumerate(chrs):
        vlines.append(vlines[i] + len(gcov[chr]))
    fig = plt.figure(figsize=(16, 4.5))
    ax = fig.add_subplot(111)
    ax.plot(win_x, win_cov, linewidth=1.5, color=(233/255, 69/255, 44/255))
    ax.set_xlim(0, vlines[-1])
    for pos in vlines[1:-1]:
        ax.axvline(x=pos, linestyle='--', linewidth=1.5, color='gray')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.spines['left'].set_position(('outward', 10))
    for label in ax.xaxis.get_ticklabels():
        label.set_fontweight('bold')
    for label in ax.yaxis.get_ticklabels():
        label.set_fontweight('bold')
    ax.tick_params(direction='out', length=6, width=2, labelsize='large', top='off', right='off', bottom='off')
    ax.set_xticks([(vlines[i] + vlines[i+1])/2 for i in xrange(len(vlines) - 1)])
    ax.set_xticklabels([extract_chrnum(chr) for chr in chrs], fontsize='large', fontweight='bold')
    #ax.set_xlabel('Chromosome', fontsize='xx-large', fontweight='bold')
    ax.set_ylabel('Coverage (base)', fontsize='large', fontweight='bold')
    fig.tight_layout()
    plt.savefig('{}.cnv.png'.format(root), dpi=300)
    plt.close(fig)

if __name__ == '__main__':
    main()

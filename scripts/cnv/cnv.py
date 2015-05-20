#!/usr/bin/env python
from __future__ import division
import os
import argparse
import pysam
import numpy as np
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-w', '--winsize', default=4000000)
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
    chrs = map(str, range(1, 23)) + ['X', 'Y']
    win_x = []
    win_cov = []
    for chr in chrs:
        start = 0
        while (start + args.winsize) <= chrlen[chr]:
            win_x.append(pos+(args.winsize/2))
            tmp = np.array(gcov[chr][start:start+args.winsize])
            tmp = tmp[~np.isnan(tmp)]
            win_cov.append(np.mean(tmp))
            start += args.winsize
            pos += args.winsize

    root = os.path.splitext(os.path.basename(args.bamfile))[0]
    vlines = [0]
    for i, chr in enumerate(chrs):
        vlines.append(vlines[i] + len(gcov[chr]))
    fig = plt.figure(figsize=(16, 4.5))
    ax = fig.add_subplot(111)
    ax.plot(win_x, win_cov, linewidth=1.5, color=(233/255, 69/255, 44/255))
    ax.set_xlim(0, vlines[-1])
    ax.set_ylim(0, 60)
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
    ax.set_xticklabels(chrs)
    ax.set_xlabel('Chromosome', fontsize='xx-large', fontweight='bold')
    ax.set_ylabel('Coverage', fontsize='xx-large', fontweight='bold')
    fig.tight_layout()
    plt.savefig('{}.cov.png'.format(root), dpi=300)
    plt.close(fig)

if __name__ == '__main__':
    main()

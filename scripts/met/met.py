#!/usr/bin/env python
from __future__ import division
import os
import re
import argparse
from itertools import izip, compress
from collections import defaultdict
import numpy as np
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt

def get_parser():
    """
    Create a parser and add arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--depth', type=int, default=4, help='minimum read depth')
    parser.add_argument('-p', '--pmtsize', type=int, default=1000, help='promoter size')
    parser.add_argument('gtf', help='GTF file')
    parser.add_argument('fasta', help='reference genome FASTA file')
    parser.add_argument('cgmap', help='CGmap file')
    return parser

def const_gtftree(gtffile):
    """
    Read a GTF file and convert it to a nested dictionary
    """
    gtftree = defaultdict(lambda: defaultdict(list))
    with open(gtffile) as infile:
        for line in infile:
            if not line.startswith('#'):
                gene_id = None
                transcript_id = None
                line = line.strip().split('\t')
                chr = line[0]
                feature = line[2]
                start = int(line[3]) - 1
                end = int(line[4])
                strand = line[6]
                attributes = line[8].split(';')
                if feature == 'exon':
                    for atb in attributes:
                        if 'gene_id' in atb:
                            gene_id = atb.strip().split()[1][1:-1]
                        elif 'transcript_id' in atb:
                            transcript_id = atb.strip().split()[1][1:-1]
                    if gene_id and transcript_id:
                        gtftree[chr][(gene_id, strand)].append((start, end))
    return gtftree

def const_ctxstr(reffile):
    """
    Construct methylation context strings from a reference genome FASTA file
    """
    with open(reffile) as infile:
        fasta = SeqIO.to_dict(SeqIO.parse(infile, 'fasta'))
        for chr in fasta:
            fasta[chr] = str(fasta[chr].seq).upper()
    ctxstr = {}
    for chr in fasta:
        ctxstr[chr] = ['-']*len(fasta[chr])
        cg = [match.start() for match in re.finditer(r'(?=(CG))', fasta[chr])]
        for pos in cg:
            ctxstr[chr][pos] = 'X'
        chg = [match.start() for match in re.finditer(r'(?=(C[ACT]G))', fasta[chr])]
        for pos in chg:
            ctxstr[chr][pos] = 'Y'
        chh = [match.start() for match in re.finditer(r'(?=(C[ACT][ACT]))', fasta[chr])]
        for pos in chh:
            ctxstr[chr][pos] = 'Z'
        rcg = [match.start()-1 for match in re.finditer(r'(?<=(CG))', fasta[chr])]
        for pos in rcg:
            ctxstr[chr][pos] = 'x'
        rchg = [match.start()-1 for match in re.finditer(r'(?<=(C[AGT]G))', fasta[chr])]
        for pos in rchg:
            ctxstr[chr][pos] = 'y'
        rchh = [match.start()-1 for match in re.finditer(r'(?<=([AGT][AGT]G))', fasta[chr])]
        for pos in rchh:
            ctxstr[chr][pos] = 'z'
    for chr in ctxstr:
        ctxstr[chr] = ''.join(ctxstr[chr])
    return ctxstr

def const_cgmap(ctxstr, cgmapfile, readdepth=4):
    """
    Construct lists of methylation levels from a CGmap file for rapid access
    """
    cgmap = {}
    with open(cgmapfile) as infile:
        for chr in ctxstr.keys():
            cgmap[chr] = ['-' for _ in xrange(len(ctxstr[chr]))]
        for line in infile:
            line = line.strip().split()
            chr = line[0]
            pos = int(line[2]) - 1 # Transfer to 0-based
            context = line[3]
            level = float(line[5])
            depth = int(line[7])
            if context in ['CG', 'CHG', 'CHH'] and depth >= readdepth:
                cgmap[chr][pos] = level
    return cgmap

def calc_bulk(ctxstr, cgmap):
    """
    Compute the global methylation level in CG/CHG/CHH
    """
    inv_ctxs = {'X': 'CG', 'Y': 'CHG', 'Z': 'CHH'}
    bulk = defaultdict(list)
    for chr in set(ctxstr) & set(cgmap):
        for tag, mlevel in izip(ctxstr[chr], cgmap[chr]):
            tag = tag.upper()
            if tag in inv_ctxs and mlevel != '-':
                bulk[inv_ctxs[tag]].append(mlevel)
    return bulk

def calc_mlevel(ctxstr, cgmap, gtftree, pmtsize=1000):
    """
    Compute the mean methylation level of promoter/gene/exon/intron/IGN in each gene
    """
    inv_ctxs = {'X': 'CG', 'Y': 'CHG', 'Z': 'CHH'}
    ign = defaultdict(list)
    mtable = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
    counter = defaultdict(lambda: defaultdict(int))
    for chr in set(ctxstr) & set(cgmap) & set(gtftree):
        mask = [1]*len(cgmap[chr])
        for (gene_id, strand) in gtftree[chr]:
            feature_mlevels = defaultdict(lambda: defaultdict(list))
            gstart = min(gtftree[chr][(gene_id, strand)])[0]
            gend = max(gtftree[chr][(gene_id, strand)])[1]
            mask[gstart:gend] = [0]*(gend - gstart)
            if strand == '+':
                for (pos, (tag, mlevel)) in enumerate(izip(ctxstr[chr][gstart-pmtsize:gstart], cgmap[chr][gstart-pmtsize:gstart])):
                    tag = tag.upper()
                    if tag in inv_ctxs and mlevel != '-':
                        feature_mlevels[inv_ctxs[tag]]['pmt'].append(mlevel)
            elif strand == '-':
                for (pos, (tag, mlevel)) in enumerate(izip(ctxstr[chr][gend:gend+pmtsize], cgmap[chr][gend:gend+pmtsize])):
                    tag = tag.upper()
                    if tag in inv_ctxs and mlevel != '-':
                        feature_mlevels[inv_ctxs[tag]]['pmt'].append(mlevel)
            for (pos, (tag, mlevel)) in enumerate(izip(ctxstr[chr][gstart:gend], cgmap[chr][gstart:gend])):
                tag = tag.upper()
                inexon = False
                if tag in inv_ctxs and mlevel != '-':
                    feature_mlevels[inv_ctxs[tag]]['gene'].append(mlevel)
                    for exon in gtftree[chr][(gene_id, strand)]:
                        if exon[0] <= pos+gstart < exon[1]:
                            feature_mlevels[inv_ctxs[tag]]['exon'].append(mlevel)
                            inexon = True
                            break
                    if not inexon:
                        feature_mlevels[inv_ctxs[tag]]['intron'].append(mlevel)
            for ctx in ['CG', 'CHG', 'CHH']:
                for feature in feature_mlevels[ctx]:
                    counter[ctx][feature] += len(feature_mlevels[ctx][feature])
                    mtable[ctx][gene_id][feature] = np.mean(feature_mlevels[ctx][feature])
        for (pos, (tag, mlevel)) in enumerate(izip(ctxstr[chr], cgmap[chr])):
            tag = tag.upper()
            if (tag in inv_ctxs) and (mask[pos] == 1) and (mlevel != '-'):
                ign[inv_ctxs[tag]].append(mlevel)
    for ctx in ign:
        print '{}: {}'.format(ctx, len(ign[ctx]))
        for feature in counter[ctx]:
            print '{} {}: {}'.format(ctx, feature, counter[ctx][feature])
        ign[ctx] = np.mean(ign[ctx])
    cg_table = pd.DataFrame(mtable['CG']).T
    cg_table = cg_table[['pmt', 'gene', 'exon', 'intron']]
    chg_table = pd.DataFrame(mtable['CHG']).T
    chg_table = chg_table[['pmt', 'gene', 'exon', 'intron']]
    chh_table = pd.DataFrame(mtable['CHH']).T
    chh_table = chh_table[['pmt', 'gene', 'exon', 'intron']]
    return ign, cg_table, chg_table, chh_table

def plot_bar(dataframe, bulk, ctx):
    colors = { 'CG': ( 38/255, 173/255,  84/255),
              'CHG': ( 44/255, 180/255, 234/255),
              'CHH': (249/255,  42/255,  54/255)}
    dataframe = dataframe*100
    plt.switch_backend('Agg')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax = dataframe.plot(ax=ax, kind='bar', grid=False, rot=0, color=colors[ctx], ylim=(0, 100))
    ax.set_ylabel('Methylation Level (%)', fontsize='xx-large', fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    #ax.spines['bottom'].set_position(('outward', 5))
    #ax.spines['left'].set_position(('outward', 5))
    ax.tick_params(direction='out', length=6, width=2, labelsize='xx-large', top='off', right='off')
    for label in ax.xaxis.get_ticklabels():
        label.set_fontweight('bold')
    for label in ax.yaxis.get_ticklabels():
        label.set_fontweight('bold')
    ax.set_title(ctx, fontsize='xx-large', weight='bold')
    #ax.axhline(y=np.mean(bulk[ctx])*100, linewidth=2, linestyle='--', color='k')
    fig.tight_layout()
    return ax

def plot_feature_mlevel(bulk, ign, cg_table, chg_table, chh_table):
    cg = cg_table.mean()
    cg = cg.set_value('genome', np.mean(bulk['CG']))
    cg = cg.set_value('IGN', ign['CG'])
    cg = cg[['genome', 'pmt', 'gene', 'exon', 'intron', 'IGN']]
    cg_ax = plot_bar(cg, bulk, 'CG')
    chg = chg_table.mean()
    chg = chg.set_value('genome', np.mean(bulk['CHG']))
    chg = chg.set_value('IGN', ign['CHG'])
    chg = chg[['genome', 'pmt', 'gene', 'exon', 'intron', 'IGN']]
    chg_ax = plot_bar(chg, bulk, 'CHG')
    chh = chh_table.mean()
    chh = chh.set_value('genome', np.mean(bulk['CHH']))
    chh = chh.set_value('IGN', ign['CHH'])
    chh = chh[['genome', 'pmt', 'gene', 'exon', 'intron', 'IGN']]
    chh_ax = plot_bar(chh, bulk, 'CHH')
    return cg_ax, chg_ax, chh_ax

def plot_bulkmean(bulk):
    bulk_mean = {}
    for ctx in ['CG', 'CHG', 'CHH']:
        bulk_mean[ctx] = np.mean(bulk[ctx])
    bulk_mean = pd.Series(bulk_mean)*100
    colors = { 'CG': ( 38/255, 173/255,  84/255),
              'CHG': ( 44/255, 180/255, 234/255),
              'CHH': (249/255,  42/255,  54/255)}
    plt.switch_backend('Agg')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax = bulk_mean.plot(ax=ax, kind='bar', grid=False, rot=0, color=[colors[ctx] for ctx in ['CG', 'CHG', 'CHH']], ylim=(0, 100))
    ax.set_ylabel('Methylation Level (%)', fontsize='xx-large', fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.tick_params(direction='out', length=6, width=2, labelsize='xx-large', top='off', right='off')
    for label in ax.xaxis.get_ticklabels():
        label.set_fontweight('bold')
    for label in ax.yaxis.get_ticklabels():
        label.set_fontweight('bold')
    fig.tight_layout()
    return ax

def plot_bulkhist(bulk):
    colors = { 'CG': ( 38/255, 173/255,  84/255),
              'CHG': ( 44/255, 180/255, 234/255),
              'CHH': (249/255,  42/255,  54/255)}
    plt.switch_backend('Agg')
    fig = plt.figure(figsize=(8, 3))
    axes = {}
    for i, ctx in enumerate(['CG', 'CHG', 'CHH']):
        if i == 0:
            axes[ctx] = fig.add_axes((0.15, 0.25, 0.25, 0.65))
            #axes[ctx] = fig.add_subplot(131)
            axes[ctx].hist(bulk[ctx], weights=np.repeat(1.0/len(bulk[ctx]), len(bulk[ctx])), color=colors[ctx])
            axes[ctx].spines['top'].set_visible(False)
            axes[ctx].spines['right'].set_visible(False)
            axes[ctx].spines['bottom'].set_linewidth(2)
            axes[ctx].spines['left'].set_linewidth(2)
            axes[ctx].spines['left'].set_position(('outward', 10))
            plt.setp(axes[ctx].get_xticklabels(), visible=False)
            axes[ctx].tick_params(axis='y', direction='out', right='off', length=6, width=2, labelsize='xx-large')
            axes[ctx].tick_params(axis='x', top='off', bottom='off')
            for label in axes[ctx].yaxis.get_ticklabels():
                label.set_fontweight('bold')
            axes[ctx].set_ylabel('Fraction', fontsize='xx-large', fontweight='bold')
        else:
            axes[ctx] = fig.add_axes((0.15 + (0.25 + 0.025) * i, 0.25, 0.25, 0.65))
            #axes[ctx] = fig.add_subplot(1, 3, i+1)
            axes[ctx].hist(bulk[ctx], weights=np.repeat(1.0/len(bulk[ctx]), len(bulk[ctx])), color=colors[ctx])
            axes[ctx].spines['top'].set_visible(False)
            axes[ctx].spines['left'].set_visible(False)
            axes[ctx].spines['right'].set_visible(False)
            axes[ctx].spines['bottom'].set_linewidth(2)
            axes[ctx].spines['left'].set_linewidth(2)
            plt.setp(axes[ctx].get_xticklabels(), visible=False)
            plt.setp(axes[ctx].get_yticklabels(), visible=False)
            axes[ctx].tick_params(top='off', bottom='off', left='off', right='off')
        axes[ctx].set_ylim(0, 1)
        axes[ctx].set_yticks(np.arange(0, 1.2, 0.2))
        axes[ctx].set_xlim(-0.025, 1.025)
        axes[ctx].set_xlabel(ctx, fontsize='xx-large', fontweight='bold')
        fig.suptitle('Methylation Level (0 -> 100%)', x=0.55, y=0.1, fontsize='xx-large', fontweight='bold')
    return fig

def calc_genomewide(ctxstr, cgmap, winsize=200000):
    inv_ctxs = {'X': 'CG', 'Y': 'CHG', 'Z': 'CHH'}
    win_mlevel = defaultdict(list)
    win_x = []
    pos = 0
    if 'chr' in ctxstr.keys()[0].lower():
        chrs = sorted(ctxstr.keys(), key=lambda s: s[3:])
    else:
        chrs = sorted(ctxstr.keys())
    #chrs = map(str, range(1, 23)) + ['X', 'Y']
    for chr in chrs:
        start = 0
        while (start + winsize) <= len(ctxstr[chr]):
            win_x.append(pos+(winsize/2))
            tmp = defaultdict(list)
            for tag, mlevel in izip(ctxstr[chr][start:start+winsize], cgmap[chr][start:start+winsize]):
                tag = tag.upper()
                if tag in inv_ctxs and mlevel != '-':
                    tmp[inv_ctxs[tag]].append(mlevel)
            for ctx in ['CG', 'CHG', 'CHH']:
                win_mlevel[ctx].append(np.mean(tmp[ctx])*100)
            start += winsize
            pos += winsize
    return win_x, win_mlevel

def plot_genomewide(ctxstr, gpos, gmlevel):
    colors = { 'CG': ( 38/255, 173/255,  84/255),
              'CHG': ( 44/255, 180/255, 234/255),
              'CHH': (249/255,  42/255,  54/255)}
    if 'chr' in ctxstr.keys()[0].lower():
        chrs = sorted(ctxstr.keys(), key=lambda s: s[3:])
    else:
        chrs = sorted(ctxstr.keys())
    #chrs = map(str, range(1, 23)) + ['X', 'Y']
    vlines = [0]
    for i, chr in enumerate(chrs):
        vlines.append(vlines[i] + len(ctxstr[chr]))
    plt.switch_backend('Agg')
    fig = plt.figure(figsize=(16, 4.5))
    ax = fig.add_subplot(111)
    ax.plot(gpos, gmlevel['CG'], color=colors['CG'], linewidth=1.5, label='CG')
    ax.plot(gpos, gmlevel['CHG'], color=colors['CHG'], linewidth=1.5, label='CHG')
    ax.plot(gpos, gmlevel['CHH'], color=colors['CHH'], linewidth=1.5, label='CHH')
    ax.set_ylim(0, 100)
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
    ax.set_xticklabels(chrs)
    ax.set_xlabel('Chromosome', fontsize='xx-large', fontweight='bold')
    ax.set_ylabel('Methylation Level (%)', fontsize='xx-large', fontweight='bold')
    ax.legend(loc='upper right', fontsize='large', frameon=False)
    fig.tight_layout()
    return ax

def main():
    parser = get_parser()
    args = parser.parse_args()
    root = os.path.splitext(os.path.basename(args.cgmap))[0]
    ctxstr = const_ctxstr(args.fasta)
    cgmap = const_cgmap(ctxstr, args.cgmap, args.depth)
    gtftree = const_gtftree(args.gtf)
    bulk = calc_bulk(ctxstr, cgmap)
    plt.switch_backend('Agg')
    bulk_ax = plot_bulkmean(bulk)
    fig = bulk_ax.get_figure()
    fig.savefig('{}.bulk.mean.png'.format(root), dpi=300)
    plt.close(fig)
    bulk_fig = plot_bulkhist(bulk)
    bulk_fig.savefig('{}.bulk.hist.png'.format(root), dpi=300)
    plt.close(fig)
    ign, cg_table, chg_table, chh_table = calc_mlevel(ctxstr, cgmap, gtftree, args.pmtsize)
    cg_table.to_csv('{}.feature.CG.txt'.format(root), sep='\t', float_format='%.3f')
    chg_table.to_csv('{}.feature.CHG.txt'.format(root), sep='\t', float_format='%.3f')
    chh_table.to_csv('{}.feature.CHH.txt'.format(root), sep='\t', float_format='%.3f')
    cg_ax, chg_ax, chh_ax = plot_feature_mlevel(bulk, ign, cg_table, chg_table, chh_table)
    fig = cg_ax.get_figure()
    fig.savefig('{}.feature.CG.png'.format(root), dpi=300)
    plt.close(fig)
    fig = chg_ax.get_figure()
    fig.savefig('{}.feature.CHG.png'.format(root), dpi=300)
    plt.close(fig)
    fig = chh_ax.get_figure()
    fig.savefig('{}.feature.CHH.png'.format(root), dpi=300)
    plt.close(fig)
    gpos, gmlevel = calc_genomewide(ctxstr, cgmap)
    gax = plot_genomewide(ctxstr, gpos, gmlevel)
    fig = gax.get_figure()
    fig.savefig('{}.genomewide.png'.format(root), dpi=300)
    plt.close(fig)

if __name__ == '__main__':
    main()

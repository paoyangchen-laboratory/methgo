#!/usr/bin/env python
from __future__ import division
import os
import re
import string
import argparse
from collections import defaultdict
import pysam
import pandas as pd
from pyfasta import Fasta

def alleleCount(baseList,refNuc):
    pos = defaultdict(int)
    neg = defaultdict(int)
    if refNuc in ['A', 'T']:
        for (base, isReverse) in baseList:
            if isReverse:    # negative strand
                neg[base] += 1
            else:    # positive strand
                pos[base] += 1
    elif refNuc == 'C': # only negative strand
        for (base, isReverse) in baseList:
            if isReverse:
                neg[base] += 1
    elif refNuc == 'G': # only positive strand
        for (base, isReverse) in baseList:
            if not isReverse:
                pos[base] += 1
    aCount = pos['A'] + neg['A']
    tCount = pos['T'] + neg['T']
    cCount = pos['C'] + neg['C']
    gCount = pos['G'] + neg['G']
    total = aCount + tCount + cCount + gCount
    posCov = sum([pos[base] for base in ['A', 'T', 'C', 'G']])
    negCov = sum([neg[base] for base in ['A', 'T', 'C', 'G']])
    return aCount, tCount, cCount, gCount, total, posCov, negCov

def snpDetermine(aCount, tCount, cCount, gCount, total, refNuc, majorAlleleFreq, buffer):
    freqA = aCount/total
    freqT = tCount/total
    freqC = cCount/total
    freqG = gCount/total
    #freqCounter = Counter(A=freqA, T=freqT, C=freqC, G=freqG)
    freqList = sorted([(freqA, 'A'), (freqT, 'T'), (freqG, 'G'), (freqC, 'C')], reverse=True)
    freqDict = {'A':freqA, 'T':freqT, 'C':freqC, 'G':freqG}
    primaryAllele = ''
    secondaryAllele = ''
    #(primary, secondary) = freqCounter.most_common(2)
    primary, secondary = freqList[:2]
    #homozygous case
    if primary[0] >= majorAlleleFreq and primary[1] != refNuc:
        #print max(freqList)
        primaryAllele = primary[1]
        secondaryAllele = 'NA'
        snpType = 'homo'
    #heterozygous case
    elif (0.5-buffer <= primary[0] <= 0.5+buffer) and (0.5-buffer <= secondary[0] <= 0.5+buffer):
        snpType = 'het'
        primaryAllele = primary[1]
        secondaryAllele = secondary[1]
    else:
        snpType = None
    return (freqList,freqDict,primaryAllele,secondaryAllele,snpType)
	
def snpOutput(bamfile, genomeFile, coverage=5,majorAlleleFreq=0.9,buffer=0.1):  
    bam = pysam.Samfile(bamfile, 'rb')
    root = os.path.splitext(os.path.basename(bamfile))[0]
    homFile = open('{}.homozygous.txt'.format(root), 'w')
    hetFile = open('{}.heterozygous.txt'.format(root), 'w')
    openFile = Fasta(genomeFile)
    for col in bam.pileup():
        chr = bam.getrname(col.tid)
        pos = col.pos
        cov = col.n
        refNuc = openFile[chr][pos]
        baseList = []
        for pileupRead in col.pileups:
            if not pileupRead.is_del:
                isReverse = pileupRead.alignment.is_reverse
                base = pileupRead.alignment.seq[pileupRead.qpos]
                baseList += [(base, isReverse)]
        (aCount,tCount,cCount,gCount,total,posCov,negCov) = alleleCount(baseList,refNuc)
        if total >= coverage:
            (freqList,freqDict,allele1,allele2,snpType) = snpDetermine(aCount,tCount,cCount,gCount,total,refNuc,majorAlleleFreq,buffer)
            eachLine = (('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n')%(chr,refNuc,pos,allele1,allele2,snpType,total,posCov,negCov,('\t'.join(map(str,(freqDict[base] for base in ['A','T','C','G']))))))
            if snpType == 'het':
			    hetFile.write(eachLine)
            elif snpType == 'homo':
				homFile.write(eachLine)	
    homFile.close()
    hetFile.close()

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('bamfile', help='input BAM file')
    parser.add_argument('-c', '--coverage', type=int, default=5, help='coverage or minimum number of reads desired')
    parser.add_argument('-m', '--majorAlleleFreq',type=float, default=0.9, help='frequency to be considered homozygous allele')
    parser.add_argument('-b', '--buffer',type=float,default=0.1, help='buffer on either side of 0.5 to be considered heterozygous allele')
    parser.add_argument('-g', '--genomeFile', help='input FASTA file')
    return parser

def main():
    parser = get_parser()
    args = parser.parse_args()
    snpOutput(args.bamfile,args.genomeFile,args.coverage,args.majorAlleleFreq,args.buffer)

if __name__ == '__main__':
    main()

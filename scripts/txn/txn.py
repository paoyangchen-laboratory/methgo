#!/usr/bin/env python
from __future__ import division
import pybedtools
import argparse
import os
import matplotlib.pyplot as plt
import numpy as np

N_TF_RANGE = 1500 * 2
STR_ALL = 'ALL'
N_RNAVG = 20

def convertDir(strDir):
    if strDir == 'C':
        return '+'
    return '-'

def splitCGmap(strCGmap,nMinCov):
    dCntxtToFout = {}
    dCntxtToCount = {}
    lOutFiles = []
    for strLine in open(strCGmap):
        lcols = strLine.rstrip().split('\t')
        strChr = lcols[0]
        strDir = convertDir(lcols[1])
        nStart = int(lcols[2])
        nEnd = nStart+1
        strCntxt = lcols[3]
        if strCntxt == '--':
            continue
        fLevel = float(lcols[5])*100
        if int(lcols[7]) < nMinCov:
            continue
        if dCntxtToFout.has_key(strCntxt) == False:
            strOutFile = strCGmap.rstrip('CGmap')+strCntxt+'.bed'
            lOutFiles.append(strOutFile)
            dCntxtToFout[strCntxt] = open(strOutFile,'w')
            dCntxtToCount[strCntxt] = 0
        dCntxtToFout[strCntxt].write('%s\t%d\t%d\t%s_%d\t%f\t%s\n'%(strChr,nStart,nEnd,strCntxt,dCntxtToCount[strCntxt]+1,fLevel,strDir))
        dCntxtToCount[strCntxt]+=1
    for strCntxt in dCntxtToFout.keys():
        dCntxtToFout[strCntxt].close()
    return lOutFiles

def getTFXY(strCntxt,strTF,dCntxtToMeth,dCntxtToCount):
    lX = []
    lY = []
    lXAvg = []
    lYAvg = []

    for i in range(N_TF_RANGE):
        nCount = dCntxtToCount[strCntxt][strTF][i]
        fSum = dCntxtToMeth[strCntxt][strTF][i]
        if nCount == 0:
            continue
        lX.append(i)
        lY.append(fSum/nCount)

    for i in np.arange(0,N_TF_RANGE,N_RNAVG):
        nCount = sum(dCntxtToCount[strCntxt][strTF][i:i+N_RNAVG])
        fSum = sum(dCntxtToMeth[strCntxt][strTF][i:i+N_RNAVG])
        if nCount == 0:
            continue
        lXAvg.append(i+(N_RNAVG/2))
        lYAvg.append(fSum/nCount)
    return lX,lY,lXAvg,lYAvg

def getCGVals(strCntxt,strCGBed,strAnnotBedPath,dCntxtToMeth,dCntxtToCount,lTFs):
    obsbed = pybedtools.BedTool(strCGBed)

    for strFile in os.listdir(strAnnotBedPath):
        lf_tfmethsum = []
        lf_tfmethcount = []
        for i in range(N_TF_RANGE):
            lf_tfmethcount.append(0)
            lf_tfmethsum.append(0.0)
        if strFile.endswith('.bed') == False:
            continue
        head,tail = os.path.split(strFile)
        if (tail.split('.')[-2] in lTFs) == False:
            continue
        strTF = strFile.split('.')[-2]
        tfbed = pybedtools.BedTool(os.path.join(strAnnotBedPath,strFile))
        tf_intersect = obsbed.intersect(tfbed,wa=True,wb=True)
        for strLine in tf_intersect:
            lcols = str(strLine).split('\t')
            nCGPos = int(lcols[1])
            nTFStart = int(lcols[7])
            nIndex = nCGPos - nTFStart
            fMeth = float(lcols[4])
            lf_tfmethcount[nIndex]+=1
            lf_tfmethsum[nIndex]+=fMeth
        if strTF in lTFs:
            dCntxtToMeth[strCntxt][strTF] = lf_tfmethsum
            dCntxtToCount[strCntxt][strTF] = lf_tfmethcount
    return

def getGlobalMeth(strBedFile):
    fVal = 0.0
    nCount = 0
    for strLine in open(strBedFile):
        fVal += float(strLine.rstrip().split('\t')[4])
        nCount+=1
    return fVal/nCount

def getCntxtToCountMeth(strBedFile,strPathToTXN,lTXN):
    dCntxtToCount = {}
    dCntxtToMeth = {}
    dCntxtGlobal = {}
    strDir, strFile = os.path.split(strBedFile)
    strCntxt = strFile.split('.')[-2]
    if dCntxtToCount.has_key(strCntxt) == False:
        dCntxtToCount[strCntxt] = {}
        dCntxtToMeth[strCntxt] = {}
    dCntxtGlobal[strCntxt] = getGlobalMeth(strBedFile)
    getCGVals(strCntxt,strBedFile,strPathToTXN,dCntxtToMeth,dCntxtToCount,lTXN)
    return dCntxtToCount,dCntxtToMeth,dCntxtGlobal

def merge(dCntxtToCount,dCntxtToMeth,dCntxtGlobal,dDCntxtToCount,dDCntxtToMeth,dDCntxtGlobal):
    for strCntxt in dCntxtToCount.keys():
        for strTF in dCntxtToCount[strCntxt].keys():
            dCntxtToCount[strCntxt][strTF+'.Inf'] = dDCntxtToCount[strCntxt][strTF]
            dCntxtToMeth[strCntxt][strTF+'.Inf'] = dDCntxtToMeth[strCntxt][strTF]

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--txnfiles', type=str, help='path to txn files')
    parser.add_argument('-l', '--txns', type=str, help='list of txn labels')
    parser.add_argument('-c', '--cgmap', type=str, help='path to cgmap file')
    parser.add_argument('-m', '--mincov', type=int, help='minimum coverage', default=4)
    return parser

def procMethBed(strMethBed,ltxns,strPathToTXN):
    strColors = 'rgbcmykw'
    colors = [(200/255, 0/255, 0/255), (53/255, 99/255, 188/255)]
    dCntxtToCount,dCntxtToMeth,dCntxtGlobal = getCntxtToCountMeth(strMethBed,strPathToTXN,ltxns)
    for strCntxt in dCntxtToCount.keys():
        linelist = []
        plt.switch_backend('Agg')
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        plt.axhline(y=dCntxtGlobal[strCntxt],color='gray',ls='dashed',linewidth=2)
        for i in range(len(ltxns)):
            strTxn = ltxns[i]
            lX,lY,lXAvg,lYAvg = getTFXY(strCntxt,strTxn,dCntxtToMeth,dCntxtToCount)
            if len(lX) > 0:
                ax.scatter(lX, lY, s=10, c=colors[i], alpha=0.1, linewidth=0)
                p0,=ax.plot(lXAvg,lYAvg, color=colors[i],label=strTxn, linewidth=2, alpha=1)
                linelist.append(p0)
        ax.legend(frameon=False)

        ax.annotate("1500bp",xy=(0.0, -1.06), xycoords='axes fraction',xytext=(.20, -1.075), textcoords='axes fraction',arrowprops=dict(arrowstyle="->",connectionstyle="arc3"),)
        ax.annotate("",xy=(0.5, -1.06), xycoords='axes fraction',xytext=(.30, -1.06), textcoords='axes fraction',arrowprops=dict(arrowstyle="->",connectionstyle="arc3"),)
        ax.annotate("1500bp",xy=(1, -1.06), xycoords='axes fraction',xytext=(.70, -1.075), textcoords='axes fraction',arrowprops=dict(arrowstyle="->",connectionstyle="arc3"),)
        ax.annotate("",xy=(0.5, -1.06), xycoords='axes fraction',xytext=(.70, -1.06), textcoords='axes fraction',arrowprops=dict(arrowstyle="->",connectionstyle="arc3"),)
        ax.xaxis.set_ticks([1500])
        ax.xaxis.set_ticklabels([''])
        ax.tick_params(direction='out', length=6, width=2, labelsize='xx-large', top='off', right='off')
        ax.set_title(strCntxt, fontsize='xx-large', weight='bold')
        plt.xlabel('Binding Site', fontsize='xx-large', weight='bold')
        plt.ylabel('Methylation Level (%)', fontsize='xx-large', weight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_linewidth(2)
        ax.spines['left'].set_linewidth(2)
        ax.spines['left'].set_position(('outward', 10))
        ax.tick_params(direction='out', length=6, width=2, labelsize='xx-large', top='off', right='off')
        for label in ax.xaxis.get_ticklabels():
            label.set_fontweight('bold')
        for label in ax.yaxis.get_ticklabels():
            label.set_fontweight('bold')
        ax.set_xlim(0,3000)
        ax.set_ylim(bottom=0)
        plt.tight_layout()
        fig.savefig(strMethBed.rstrip('bed')+'txn.png', dpi=300)

def main():
    parser = get_parser()
    args = parser.parse_args()

    fRandAvg = 0

    #split CGmap files into bed files by context
    lOutFiles = splitCGmap(args.cgmap,args.mincov)

    #process files
    for strMethBedFile in lOutFiles:
        procMethBed(strMethBedFile,args.txns.split(','),args.txnfiles)

if __name__ == '__main__':
    main()

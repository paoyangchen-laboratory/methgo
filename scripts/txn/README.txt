TXN
The script will generate meta plots at user selected transcription factor targets for CG, CHG, and CHH contexts.

Output Formats:
.CG.txn.pdf
.CHG.txn.pdf
.CHH.txn.pdf

Output Description:
The average methylation level of the sample file is displayed as a dashed grey line
The middle of the x-axis is the center of the transcription binding site 
The start of the x-axis displays the methylation 1500bp upstream of the transcription factor binding center.
The end of the x-axis displays the methylation level 1500bp downstream of the transcription factor binding center.

Prerequisite: 
pybedtools, argparse, matplotlib, numpy

Usage:
$python runModule.py txn --txnfiles <path to bed files of transcription factor binding centers padded by 1500bp> --txns <list of transcription factor binding sites to plot> --cgmap <BSSeeker CGmap methylation file>

Example:
$python runModule.py txn --txnfiles ./hg19_txn --txns RAND,ATF1 --cgmap test.CGmap

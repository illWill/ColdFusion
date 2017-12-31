import imp
import sys
import os
import datetime
import time
import subprocess

import readTools
import blatTools
imp.reload(readTools)
imp.reload(blatTools)

def alignAndCall(fileName,sourceDir,destinationDir):
    print("processing {}".format(fileName))
    blatPath = '/Users/wjeck/blat/blat'
    genomePath = '/Users/wjeck/Genomes/hg19/hg19.2bit'
    fastaFile = destinationDir + "/" + fileName + ".fa"
    reads = readTools.loadNanoporeFastq(sourceDir + "/" + fileName)
    fasta = open(fastaFile,"w")
    fasta.write(reads.asFasta())
    fasta.close()
    blatFile = destinationDir + "/" + fileName + ".blat"
    runBlat = '{} -stepSize=5 -repMatch=2253 -minScore=0 -minIdentity=0 "{}" "{}" "{}"'.format(blatPath,genomePath, fastaFile, blatFile)
    print(runBlat)
    subprocess.call(runBlat,shell=True)

startTime = datetime.datetime.now()
processed = []
a = []
sourceDir = '/Users/wjeck/Dropbox (Partners HealthCare)/MGH Nanopore/ComboRUN/fastq/pass/'
destinationDir = '/Users/wjeck/Dropbox (Partners HealthCare)/MGH Nanopore/ComboRUN/fastq/pass/blats'

if not os.path.exists(destinationDir):
    os.makedirs(destinationDir)

slept = 0
while slept < 60:
    files = os.listdir(sourceDir)
    unprocessed = [f for f in files if f not in processed and f.endswith('.fastq')]
    print(unprocessed)
    if len(unprocessed) < 2: 
        time.sleep(10)
        slept += 10
    else:
        slept = 0
        unprocessedTimes = [os.path.getmtime(sourceDir + "/" + f) for f in unprocessed]
        a = list(zip(unprocessedTimes,unprocessed))
        a.sort()
        alignAndCall(a[0][1] , sourceDir, destinationDir)
        processed.append(a[0][1])

unprocessed = [f for f in files if f not in processed and f.endswith('.fastq')]
alignAndCall(unprocessed[0] , sourceDir, destinationDir)
processed.append(a[0][1])
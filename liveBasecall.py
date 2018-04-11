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

startCode = time.time()

def callLive(subdirectory,sourceRoot,targetDir):
    sourcePath = sourceRoot + "/" + subdirectory
    albacorePy = 'read_fast5_basecaller.py --flowcell FLO-MIN107 --kit SQK-LSK108 -o fastq -q 0 --worker_threads 8 --save_path ' + targetDir + " --input " + sourcePath 
    albacoreBsubOptions = "-q big-multi -n 4 -R 'rusage[mem=16000]'"
    albacoreJobcode = startCode + ":" + subdirectory + ":albacore"
    albacoreRunline = 'bsub -J "{}" {} "{}"'.format(albacoreJobcode, albacoreBsubOptions, albacorePy)
    
    deMultiplexer = 'python3 demultiplex.py ' + targetDir  #TODO: Demultiplex directories rather than files
    deMultiplexerJobCode = startCode + ":" + subdirectory + ":demultiplexer"
    deMultiplexerRunline = 'bsub -J "{}" -w "{}" "{}"'.format(deMultiplexerJobCode,albacoreJobcode,deMultiplexer)
    
    #bwaCall = '/PHShome/wrj2/bwa/bwa mem /PHShome/wrj2/Genomes/hg19/coldfuse/coldFusion3.fasta' #needs redirect
    print(albacoreRunline)
    print(deMultiplexerRunline)
    #subprocess.call(albacoreRunline, shell=True)
    #subprocess.call(deMultiplexerRunline, shell=True)

startCode = str(time.time())
processed = []
a = []
sourceDir = sys.argv[1]
destinationDir = sys.argv[2]

if not os.path.exists(destinationDir):
    os.makedirs(destinationDir)

slept = 0
while slept < 60:
    files = os.listdir(sourceDir)
    unprocessedCandidates = [f for f in files if os.path.isdir(sourceDir + "/" + f) and f.isdigit()]
    unprocessed = [f for f in unprocessedCandidates if f not in processed]
    while len(unprocessed) >= 2: 
        slept = 0
        print(unprocessed)
        unprocessedTimes = [os.path.getmtime(sourceDir + "/" + f) for f in unprocessed]
        a = list(zip(unprocessedTimes,unprocessed))
        a.sort()
        targetDir = destinationDir + "/" + a[0][1]
        if not os.path.exists(targetDir):
            os.makedirs(targetDir)
        callLive(a[0][1] , sourceDir, targetDir)
        processed.append(a[0][1])
        unprocessed = [f for f in unprocessedCandidates if f not in processed]
    time.sleep(10)
    slept += 10

unprocessed = [f for f in files if f not in processed and os.path.isdir(sourceDir + "/" + f) and f.isdigit()]
for f in unprocessed:
    callLive(unprocessed[0] , sourceDir, destinationDir)
    processed.append(a[0][1])

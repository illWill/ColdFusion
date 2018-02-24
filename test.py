import imp
import sys
import os
import datetime
import re
import operator
import itertools
from collections import Counter

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

import readTools
import blatTools
import alignment
import cpp_function_wrappers
imp.reload(readTools)
imp.reload(blatTools)
imp.reload(alignment)
imp.reload(cpp_function_wrappers)

def printRed(text):
    RED   = "\033[1;31m"  
    BLUE  = "\033[1;34m"
    CYAN  = "\033[1;36m"
    GREEN = "\033[0;32m"
    RESET = "\033[0;0m"
    BOLD    = "\033[;1m"
    REVERSE = "\033[;7m"
    return RED + text + RESET

def printGreen(text):
    RED   = "\033[1;31m"  
    BLUE  = "\033[1;34m"
    CYAN  = "\033[1;36m"
    GREEN = "\033[0;32m"
    RESET = "\033[0;0m"
    BOLD    = "\033[;1m"
    REVERSE = "\033[;7m"
    return GREEN + text + RESET

def bestAdapterAlign(read,adapterDict):
    '''Provided a read sequence and a dictionary of adapters, it identifies the 
    best adapter alignment and its distance from the other adapter alignments'''
    align = (0.0, 0.0, 0, 0)
    gap = 0
    bestName = ''
    for name,adapt in adapterDict.items():
        align2 = cpp_function_wrappers.adapter_alignment_simple(read, adapt)
        if align2[0] > align[0]:
            bestName = name
            gap = align2[0] - align[0]
            align = align2
        elif align2[0] == align[0]:
            gap = 0
    return (align,gap,bestName)

def readSplitAlignment(seq,align):
    '''Splits up a read into overall_identity, prefix, adapter, and suffix, based on the provided alignment'''
    return (align[0], seq[0:align[2]], seq[align[2]:align[3]], seq[align[3]:])

def readSplitAdapter(seq,adapter):
    '''Identifies the best alignment for an adapter, and considers both and orientations'''
    alignf = cpp_function_wrappers.adapter_alignment_simple(seq, adapter)
    alignr = cpp_function_wrappers.adapter_alignment_simple(seq, alignment.revcomp(adapter))
    align = max(alignf,alignr)
    return readSplitAlignment(seq,align)

def printSegmentedRead(readText):
    '''Takes a segmented read with tags and prints it out to highlight adapters'''
    RED   = "\033[1;31m"
    RESET = "\033[0;0m"
    readText = re.sub(r'(<[^/])',RED + r'\1',readText)
    readText = re.sub(r'(</[^>]+>)', r'\1' + RESET,readText)
    print(readText + "\n")

def parseReadFile(fn,adapterFile):
    prefix = ".".join(fn.split(".")[:-1])
    reads = readTools.loadNanoporeFastq(fn)
    
    adapterFasta = readTools.loadFasta(adapterFile)
    adapters = {read.name: read.sequence for read in adapterFasta.readList}
    adapters.update ( {'~' + read.name : alignment.revcomp(read.sequence) for read in adapterFasta.readList} )
    
    adapterLead = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
    p5 = 'AATGATACGGCGACCACCGAGATCTACAC'
    p7 = 'CAAGCAGAAGACGGCATACGAGAT'
    
    p7withBarcode = 'CAAGCAGAAGACGGCATACGAGAT TCGCCTTA             GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT' #not used, but for reference
    p5withBarcode = 'AATGATACGGCGACCACCGAGATCTACAC TAGATCGCNNWNNWNN ACACTCTTTCCCTACACGACGCTCTTCCGATCT' #not used, but for reference
    
    p7Lead = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
    p5Lead = 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'
    
    adapters.update({'p5':p5, 'p7':p7, 'p5Lead':p5Lead,'~p5':alignment.revcomp(p5), '~p7':alignment.revcomp(p7), '~p5Lead':alignment.revcomp(p5Lead)})
    
    def recursiveReadSplit(read):
      '''Takes in read sequence and recursively identifies adapter matches and adds HTML style tags'''
      # p5align = readSplitAdapter(read, p5)
      # if(p5align[0] > 75):
      #     return recursiveReadSplit(p5align[1]) + printGreen('<p5>' + p5align[2] + '</p5>') + recursiveReadSplit(p5align[3]) 
      # p7align = readSplitAdapter(read, p7)
      # if(p7align[0] > 75):
      #     return recursiveReadSplit(p7align[1]) + printGreen('<p7>' + p7align[2] + '</p7>') + recursiveReadSplit(p7align[3])
      adaptersAlign = bestAdapterAlign(read,adapters)
      adaptersSplit = readSplitAlignment(read,adaptersAlign[0])
      if(adaptersSplit[0] > 75):
          prefix = recursiveReadSplit(adaptersSplit[1]) 
          match = '<' + adaptersAlign[2] + '>' + adaptersSplit[2] + '</' + adaptersAlign[2] + '>'
          suffix = recursiveReadSplit(adaptersSplit[3]) 
          return prefix + match + suffix
      p7LeadAlign = readSplitAdapter(read, p7Lead)
      if(p7LeadAlign[0] > 75):
          return recursiveReadSplit(p7LeadAlign[1]) + '<p7Lead>' + p7LeadAlign[2] + '</p7Lead>' + recursiveReadSplit(p7LeadAlign[3])
      return read
    
    allParsed = []
    coded = {}
    total = 0
    for read in reads.readList:
        read_split = recursiveReadSplit(read.sequence)
        allParsed.append(read_split)
        #printSegmentedRead( read_split )
        find_forward = re.finditer(r'</(?P<barcode>BC\d\d)>(?P<sequence>[ACTG]+)<~p5Lead>',read_split)
        find_forward2 = re.finditer(r'</(?P<barcode>BC\d\d)>(?P<sequence>[ACTG]+)$',read_split)
        find_double = re.finditer(r'</(?P<barcode>BC\d\d)>(?P<sequence>[ACTG]+)<~(?P=barcode)>',read_split)
        find_reverse = re.finditer(r'</p5Lead>(?P<sequence>[ACTG]+)<~(?P<barcode>BC\d\d)>',read_split)
        find_reverse2 = re.finditer(r'^(?P<sequence>[ACTG]+)<~(?P<barcode>BC\d\d)>',read_split)
        for (i,match) in enumerate(itertools.chain(find_forward, find_forward2, find_double, find_reverse, find_reverse2)):
            total+=1
            if total%1000 == 0:
                print(total)
            if match.group('barcode') in coded:
                coded[match.group('barcode')].append( (read.name + '_FRAGMENT' + str(i), match.group('sequence')))
            else:
                coded[match.group('barcode')] = [ (read.name + '_FRAGMENT' + str(i), match.group('sequence')) ]
    for a in coded.keys():
        f = open(prefix + '.' + a + '.fasta','w')
        for (name,sequence) in coded[a]:
            junk = f.write('>' + name + "\n" + sequence + "\n")
        f.close()

#x1 = [sorted(alignment) for alignment in aligns]
#maxed  = [a[11][0] for a in x1]
#gapped = [a[11][0] - a[10][0] for a in x1]
#n, bins, patches = plt.hist(gapped, 50, normed=1, facecolor='green', alpha=0.75)

'''System for processing a coldfusion mapped sam file... preliminary'''
def findBestFusion(samFileName):
    try:
        sam = open(samFileName,'r')
    except IOError:
        return({})
    readHits = {}
    for line in sam:
        if line[0] == '@':
            continue
        fields = line.split("\t")
        if fields[0] in readHits:
            if not fields[2] in readHits[fields[0]]:
                readHits[fields[0]].append(fields[2])
        else:
            readHits[fields[0]] = [fields[2]]
    sam.close()
    hits = [tuple(sorted(readHits[r])) for r in readHits.keys() if len(readHits[r]) >= 2]
    hit_dict = Counter(hits)
    import operator
    hit_sorted = sorted(hit_dict.items(), key=lambda x: -x[1])
    best_hit = hit_sorted[0]
    #print(best_hit)
    #reads_on_best_hit = [r for r in readHits.keys() if len(readHits[r]) >= 2 and best_hit[0] == (readHits[r][0], readHits[r][1])]
    return(hit_sorted)

def narrowRecursiveSplitForDoubledFasta():   ####NOT CURRENTLY USED, I don't know what to do with it. 
    coded2 = {}
    total = 0
    for read in reads.readList:
        total+=1
        if total%1000 == 0:
            print(total)
        read_split = recursiveReadSplit(read.sequence)
        #printSegmentedRead( read_split )
        find_double = re.finditer(r'</(?P<barcode>BC\d\d)>(?P<sequence>[ACTG]+)<~(?P=barcode)>',read_split)
        for (i,match) in enumerate(find_double):
            if match.group('barcode') in coded2:
                coded2[match.group('barcode')].append( (read.name + '_FRAGMENT' + str(i), match.group('sequence')))
            else:
                coded2[match.group('barcode')] = [ (read.name + '_FRAGMENT' + str(i), match.group('sequence')) ]
    for a in coded2.keys():
        f = open(a + '.doubled.fasta','w')
        for (name,sequence) in coded[a]:
            junk = f.write('>' + name + "\n" + sequence + "\n")
        f.close()


def groupQuery(blatlet):
    qgroups = []
    for hit in blatlet.records:
        if len(qgroups) == 0:
            qgroups.append( [hit.qStart,hit.qEnd,blat([hit])] )
            next
        for group in qgroups:
            if hit.qStart > group[1] or hit.qEnd < group[0]: #nonoverlapping in query
                qgroups.append( [hit.qStart,hit.qEnd,blat([hit])] )
            else: 
                group[0] = min(group[0],hit.qStart)
                group[1] = max(group[1],hit.qEnd)
                group[2].records.append(hit)
    return qgroups

if len(sys.argv) > 1:
    fileName = sys.argv[1]
    adapterFile= r'/PHShome/wrj2/Scripts/ColdFusion/NanoporeBarcodes.fa'
    parseReadFile(fileName,adapterFile)
else:
    BClist = ['01','02','03','04','05','06','07','08','09','10','11','12']
    fnList = ['RAW.BC' + i + '.sam' for i in BClist]
    for i in fnList:
        print(i)
        hits = findBestFusion(i)
        for hit in hits:
            if hit[1] > 2:
                print(hit)



class BlatRecord:
    def __init__(self, s):
        array = s.split()
        self.match = int(array[0])
        self.mismatch = int(array[1])
        self.repmatch = int(array[2])
        self.ns = int(array[3])
        self.qGapCount = int(array[4])
        self.qGapBases = int(array[5])
        self.tGapCount = int(array[6])
        self.tGapBases = int(array[7])
        self.strand = array[8]
        self.qName = array[9]
        self.qSize = int(array[10])
        self.qStart = int(array[11])
        self.qEnd = int(array[12])
        self.qRegion = Region(self.qName,self.qStart,self.qEnd)
        self.tName = array[13]
        self.tSize = int(array[14])
        self.tStart = int(array[15])
        self.tEnd = int(array[16])
        self.tRegion = Region(self.tName,self.tStart,self.tEnd)
        self.blockCount = int(array[17])
        self.blockSizes = array[18].split(",")
        self.qStarts = array[19].split(",")
        self.tStarts = array[20].split(",")
        return
    
    def __str__(self):
        return "\t".join([str(self.match),
        str(self.mismatch),
        str(self.repmatch),
        str(self.ns),
        str(self.qGapCount),
        str(self.qGapBases),
        str(self.tGapCount),
        str(self.tGapBases),
        self.strand,
        self.qName,
        str(self.qSize),
        str(self.qStart),
        str(self.qEnd),
        self.tName,
        str(self.tSize),
        str(self.tStart),
        str(self.tEnd),
        str(self.blockCount),
        ",".join(self.blockSizes),
        ",".join(self.qStarts),
        ",".join(self.tStarts)])

class Region:
    def __init__(self, name, start, end):
        if start > end:
            raise ValueError("Start location > end location")
        self.name = name
        self.start = start
        self.end = end
    
    def __len__(self):
        return self.end - self.start
    
    def overlaps(self,region):
        return self.name==region.name and self.start < region.end and self.end > region.start
    
    def __contains__(self,region):
        return self.name==region.name and self.start <= region.start and region.end <= self.end

class Blat:
    
    def __init__(self, rec):
        self.records = rec
        self.headerText = '''psLayout version 3

match	mis- 	rep. 	N's	Q gap	Q gap	T gap	T gap	strand	Q        	Q   	Q    	Q  	T        	T   	T    	T  	block	blockSizes 	qStarts	 tStarts
     	match	match	   	count	bases	count	bases	      	name     	size	start	end	name     	size	start	end	count
---------------------------------------------------------------------------------------------------------------------------------------------------------------'''
        return
    
    def __str__(self):
        return "\n".join([str(r) for r in self.records])
    
    def __len__(self):
        return len(self.records)
    
    def __bool__(self):
        return bool(self.records)
    
    def qerryNameHits(self,queryName):
        return Blat([r for r in self.records if r.qName==readName])
    
    def targetNameHits(self,targetName):
        return Blat([r for r in self.records if r.tName==targetName])
        
    def queryList(self):
        return [r.qName for r in self.records]
    
    def targetList(self):
        return [r.tName for r in self.records]
    
    def genQueryGroups(self):
        i = 0
        while i < len(self.records):
            queryHits = []
            qName = self.records[i].qName
            while i < len(self.records) and self.records[i].qName == qName:
                queryHits.append(self.records[i])
                i += 1
            yield Blat(queryHits)
    
    def cull(self,scoreCutoff):
        i = 0
        while i < len(self.records):
            r = self.records[i]
            if r.match < scoreCutoff:
                self.records.pop(i)
            else:
                i += 1
    
    def targetRegionContained(self,region):
        return [contained for contained in self.records if contained.tRegion in region]
    
    def queryRegionContained(self,region):
        return [contained for contained in self.records if contained.qRegion in region]
    
    def targetRegionHits(self):
        return [hit for hit in self.records if hit.tRegion.overlaps(region)]
    
    def queryRegionHits(self):
        return [hit for hit in self.records if hit.qRegion.overlaps(region)]

def loadBlat(fileName):
    file = open(fileName,'r')
    text = file.read()
    rec = text.split("\n")[5:]
    records = []
    for r in rec:
        if r == "":
            continue
        records.append(BlatRecord(r))
    return(Blat(records))

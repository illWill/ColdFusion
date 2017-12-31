import re

class ReadGroup:
    def __init__(self,l):
        self.readList = l
    
    def averageLength(self):
        return sum([len(r.sequence) for r in self.readList])/len(self.readList)
    
    def findRead(self,readName): 
        for r in self.readList:
            if r.name == readName:
                return r
        return None
    
    def asFasta(self):
        s = ""
        for r in self.readList:
            s = s + r.asFasta()
        return s

class Read:
    def __init__(self,readName,readSequence,readQuality):
        self.name = readName
        self.sequence = re.sub('\s','', readSequence,count=0)
        self.quality = re.sub('\s','',readQuality,count=0)
    
    def __str__(self):
        return self.name + "\n" + self.sequence + "\n"
    
    def asFasta(self):
        return ">{}\n{}\n".format(self.name,self.sequence)

class NanoporeRead(Read):
    def __init__(self,readName,readSequence,readQuality):
        '''Supports reads with metadata in the readName line'''
        dat = readName.split(" ")
        self.meta = dict([r.split("=") for r in dat[1:]])
        Read.__init__(self, dat[0],readSequence,readQuality)
    
    def __str__(self):
        return self.name + " " + " ".join([key + "=" + self.meta[key] for key in self.meta]) + "\n" + self.sequence

class NanoporeReadGroup(ReadGroup):
    '''Stub class at the moment'''

def loadFasta(fileName):
    mfile = open(fileName,'r')
    text = mfile.read()
    fast = text.split("\n")
    records = []
    i = 0
    while i < len(fast):
        readName = fast[i][1:]
        i += 1
        readSequence = ""
        while i < len(fast) and (len(fast[i]) == 0 or fast[i][0] != '>')    :   
            readSequence += fast[i]
            i += 1
        records.append(Read(readName,readSequence,''))
    return(ReadGroup(records))

def loadNanoporeFasta(fileName):
    mfile = open(fileName,'r')
    text = mfile.read()
    fast = text.split("\n")
    records = []
    i = 0
    while i < len(fast):
        readName = fast[i][1:]
        i += 1
        readSequence = ""
        while i < len(fast) and (len(fast[i]) == 0 or fast[i][0] != '>')    :   
            readSequence += fast[i]
            i += 1
        records.append(NanoporeRead(readName,readSequence,''))
    return(ReadGroup(records))

def loadNanoporeFastq(fileName):
    mfile = open(fileName,'r')
    text = mfile.read()
    fast = text.split("\n")
    records = []
    i = 0
    while i < len(fast) and len(fast[i]) > 0:
        readName = fast[i][1:]
        readSequence = fast[i+1]
        readQuality = fast[i+3]
        records.append(NanoporeRead(readName, readSequence, readQuality))
        i += 4
    return NanoporeReadGroup(records)


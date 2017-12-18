class Primer (object):
    def __init__(self, kind, sequnce):
        self.kind = kind
        self.sequnce = sequnce

    def aCounter(self):
        aCount = 0
        for index in self.sequnce:
            if index == 'A':
                aCount += 1
        return aCount

    def tCounter(self):
        tCount = 0
        for index in self.sequnce:
            if index == 'T':
                tCount += 1
        return tCount

    def gCounter(self):
        gCount = 0
        for index in self.sequnce:
            if index == 'G':
                gCount += 1
        return gCount

    def cCounter(self):
        cCount = 0
        for index in self.sequnce:
            if index == 'C':
                cCount += 1
        return cCount

    def precentGC(self):
        prcentGC = int(((self.cCounter()+self.gCounter())/len(self.sequnce))*100)
        return  prcentGC

    def primerTm(self):
        tm = int((64.9 +41 * (self.gCounter() + self.cCounter() - 16.4))/(len(self.sequnce)))
        return tm
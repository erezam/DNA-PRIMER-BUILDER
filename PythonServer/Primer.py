from __future__ import division


class Primer (object):
    def __init__(self, kind, sequnce , start_index):
        self.kind = kind
        self.sequnce = sequnce
        self.length = len(sequnce)
        self.score = 0
        self.start_index = start_index

    def a_counter(self):
        aCount = 0
        for index in self.sequnce:
            if index == 'A':
                aCount += 1
        return aCount

    def t_counter(self):
        tCount = 0
        for index in self.sequnce:
            if index == 'T':
                tCount += 1
        return tCount

    def g_counter(self):
        gCount = 0
        for index in self.sequnce:
            if index == 'G':
                gCount += 1
        return gCount

    def c_counter(self):
        cCount = 0
        for index in self.sequnce:
            if index == 'C':
                cCount += 1
        return cCount

    def precent_gc(self):
        prcentGC = int(((self.c_counter()+self.g_counter())/len(self.sequnce))*100)
        return prcentGC

    def primer_tm(self):
        tm = int((64.9 + 41 * (self.g_counter() + self.c_counter() - 16.4)/(len(self.sequnce))))
        return tm

    def add_score(self, add):
        self.score += add

    def printPrimer(self):
        print "Kind: %s , Seq: %s , Tm : %s , GC : %s , Score : %s , Start index: %s" % \
              (self.kind, self.sequnce,self.primer_tm(), self.precent_gc(), self.score, self.start_index)
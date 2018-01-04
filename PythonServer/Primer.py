from __future__ import division

import jason as jason

config = jason.load(open("config.json"))

class Primer (object):
    def __init__(self, id, kind, sequnce , start_index, pair_id=0):
        self.id = id
        self.kind = kind
        self.sequnce = sequnce
        self.length = len(sequnce)
        self.score = 0
        self.start_index = start_index
        self.pair_id = pair_id

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

# ==============================================================================================

    def precent_gc(self):
        prcentGC = int(((self.c_counter()+self.g_counter())/len(self.sequnce))*100)
        return prcentGC

# ==============================================================================================

    def primer_tm(self):
        tm = int((64.9 + 41 * (self.g_counter() + self.c_counter() - 16.4)/(len(self.sequnce))))
        return tm

# ====================================calculate the score of the primer, based on multiplication of the score of each parameter ===========

    def get_primer_score(self):
        tmp_score = self.get_tm_score() * self.get_gc_score() * self.get_length_score()
        self.score = tmp_score

# ================================calculate score for each parameter============================================
    def get_tm_score(self):
        Avg = config["Tm"]["Avg"]
        min = config["Tm"]["Min"]
        max = config["Tm"]["Max"]
        std = config["Tm"]["STDDEV"]
        if self.primer_tm() < min or self.primer_tm() > max:
                return 0
        elif self.primer_tm() == Avg:
            return self.primer_tm()
        elif self.primer_tm() + std == Avg or self.primer_tm() - std == Avg:
            return float(self.primer_tm())*0.95
        elif self.primer_tm()+ 1.5*float(std) == Avg or self.primer_tm()- 1.5*float(std) == Avg:
            return float(self.primer_tm())*0.9
        elif self.primer_tm() + 2 * std == Avg or self.primer_tm() - 2 * std == Avg:
            return float(self.primer_tm()) * 0.85
        else:
            return float(self.primer_tm()) * 0.8

    def get_gc_score(self):
        Avg = config["GC Percent"]["Avg"]
        min = config["GC Percent"]["Min"]
        max = config["GC Percent"]["Max"]
        std = config["GC Percent"]["STDDEV"]
        if self.precent_gc() < min or self.precent_gc() > max:
            return 0
        elif self.precent_gc() == Avg:
            return self.precent_gc()
        elif self.precent_gc() + std == Avg or self.precent_gc() - std == Avg:
            return float(self.primer_tm()) * 0.95
        elif self.precent_gc() + 1.5 * float(std) == Avg or self.precent_gc() - 1.5 * float(std) == Avg:
            return float(self.precent_gc()) * 0.9
        elif self.precent_gc() + 2 * std == Avg or self.precent_gc() - 2 * std == Avg:
            return float(self.precent_gc()) * 0.85
        else:
            return float(self.precent_gc()) * 0.8

    def get_length_score(self):
        Avg = config["Length"]["Avg"]
        min = config["Length"]["Min"]
        max = config["Length"]["Max"]
        std = config["Length"]["STDDEV"]
        if self.length < min or self.length > max:
            return 0
        elif self.length == Avg:
            return self.primer_tm()
        elif self.length + std == Avg or self.length - std == Avg:
            return float(self.length) * 0.95
        elif self.length + 1.5 * float(std) == Avg or self.length - 1.5 * float(std) == Avg:
            return float(self.length) * 0.9
        elif self.length + 2 * std == Avg or self.length - 2 * std == Avg:
            return float(self.length) * 0.85
        else:
            return float(self.length) * 0.8

            # ================================print primer============================================
    def printPrimer(self):
        print "Id:%s, Pair id:%s, Kind: %s , Seq: %s , Tm : %s , GC : %s , Score : %s , Start index: %s" % \
              (self.id, self.pair_id, self.kind, self.sequnce,self.primer_tm(), self.precent_gc(), self.score, self.start_index)
from __future__ import division
import json

config = json.load(open("config.json"))

class Primer (object):
    def __init__(self, id, kind, sequnce , start_index, pair_id=0):
        self.id = id
        self.kind = kind
        self.sequnce = sequnce
        self.length = len(sequnce)
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

# ========== calculate the score of the primer, based on multiplication of the score of each parameter ===========

    def get_primer_score(self):
        return float(self.get_tm_score()) * float(self.get_gc_score()) * float(self.get_length_score())


# ================================ calculate score for each parameter============================================

# =========================== get Tm score ======================================================================
    def get_tm_score(self):
        Avg = int(config["Tm"]["Avg"])
        min = int(config["Tm"]["Min"])
        max = int(config["Tm"]["Max"])
        std = float(config["Tm"]["STDDEV"])
        if self.primer_tm() < min or self.primer_tm() > max:
            return 0
        elif self.primer_tm() + 0.5 * float(std) >= Avg >= self.primer_tm() - 0.5 * float(std):
            return 1
        elif self.primer_tm() + 1 * float(std) >= Avg >= self.primer_tm() - 1 * float(std):
            return 0.95
        elif self.primer_tm() + 1.5 * float(std) >= Avg >= self.primer_tm() - 1.5 * float(std):
            return 0.9
        elif self.primer_tm() + 2 * float(std) >= Avg >= self.primer_tm() - 2 * float(std):
            return 0.85
        else:
            return 0.8

# =========================== get Gc score ======================================================================
    def get_gc_score(self):
        Avg = int(config["GC Percent"]["Avg"])
        min = int(config["GC Percent"]["Min"])
        max = int(config["GC Percent"]["Max"])
        std = float(config["GC Percent"]["STDDEV"])
        if self.precent_gc() < min or self.precent_gc() > max:
            return 0
        elif self.precent_gc() + 0.5 * float(std) >= Avg >= self.precent_gc() - 0.5 * float(std):
            return 1
        elif self.precent_gc() + 1 * float(std) >= Avg >= self.precent_gc() - 1 * float(std):
            return 0.95
        elif self.precent_gc() + 1.5 * float(std) >= Avg >= self.precent_gc() - 1.5 * float(std):
            return 0.9
        elif self.precent_gc() + 2 * float(std) >= Avg >= self.precent_gc() - 2 * float(std):
            return 0.85
        else:
            return 0.8

# =========================== get len score ======================================================================
    def get_length_score(self):
        Avg = int(config["Length"]["Avg"])
        min = int(config["Length"]["Min"])
        max = int(config["Length"]["Max"])
        std = float(config["Length"]["STDDEV"])
        if self.length < min or self.length > max:
            return 0
        elif self.length + 0.5 * float(std) >= Avg >= self.length - 0.5 * float(std):
            return 1
        elif self.length + 1 * float(std) >= Avg >= self.length - 1 * float(std):
            return 0.95
        elif self.length + 1.5 * float(std) >= Avg >= self.length - 1.5 * float(std):
            return 0.9
        elif self.length + 2 * float(std) >= Avg >= self.length - 2 * float(std):
            return 0.85
        else:
            return 0.8

# ================================print primer============================================
    def printPrimer(self):
        print "Id:%s, Pair id:%s, Kind: %s , Seq: %s , Tm : %s , GC : %s , Start index: %s" % \
              (self.id, self.pair_id, self.kind, self.sequnce,self.primer_tm(), self.precent_gc(), self.start_index)
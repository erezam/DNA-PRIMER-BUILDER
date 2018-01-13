import json

from PythonServer.Primer import Primer

config = json.load(open("config.json"))

class Primer_set (object):
    def __init__(self, forward_primer, reverse_primer):
        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer

    def tm_dif(self):
        tm_difference = abs(self.forward_primer.primer_tm() - self.reverse_primer.primer_tm())
        return tm_difference

    def set_print(self):
        print "Forward id: %s - Reverse id:%s" %(self.forward_primer.id , self.reverse_primer.id)

    def get_amplicon_length(self):
        return self.reverse_primer.start_index - self.forward_primer.start_index

    def write_to_file(self, file_name):
        write_file = open(file_name, "a")
        write_file.write(
            "----------------------------------------------------------------------------------------------\n" +
            "FORWARD :\n" +
            "Sequence: %s length: %s Tm: %s GC: %s start index: %s \n" % (
            self.forward_primer.sequence, self.forward_primer.length, self.forward_primer.primer_tm(),
            self.forward_primer.precent_gc(), self.forward_primer.start_index) +
            "REVERSE:\n" +
            "Sequence: %s length: %s Tm: %s GC: %s start index: %s \n" % (
            self.reverse_primer.sequence, self.reverse_primer.length, self.reverse_primer.primer_tm(),
            self.reverse_primer.precent_gc(), self.reverse_primer.start_index) +
            "tm difference between F & R: %s \n" % (self.tm_dif()) +
            "Amplicon length: %s \n" % (self.get_amplicon_length()) +
            "Score: %s \n" % (self.get_primers_set_score()))
        write_file.close()

    def get_primers_set_score(self):
        return int((self.forward_primer.get_primer_score() *
                    self.reverse_primer.get_primer_score() * self.get_amplicon_length_score())*100)


    def get_amplicon_length_score(self):
        Avg = int(config["Amplicon Length"]["Avg"])
        min = int(config["Amplicon Length"]["Min"])
        max = int(config["Amplicon Length"]["Max"])
        std = float(config["Amplicon Length"]["STDDEV"])
        if self.get_amplicon_length() < min or self.get_amplicon_length() > max:
            return 0
        elif self.get_amplicon_length() + 0.5 * float(std) >= Avg >= self.get_amplicon_length() - 0.5 * float(std):
            return 1
        elif self.get_amplicon_length() + 1 * float(std) >= Avg >= self.get_amplicon_length() - 1 * float(std):
            return 0.95
        elif self.get_amplicon_length() + 1.5 * float(std) >= Avg >= self.get_amplicon_length() - 1.5 * float(std):
            return 0.9
        elif self.get_amplicon_length() + 2 * float(std) >= Avg >= self.get_amplicon_length() - 2 * float(std):
            return 0.85
        else:
            return 0.8


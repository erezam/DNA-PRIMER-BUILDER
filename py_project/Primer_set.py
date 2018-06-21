import json
import ScorePredictor

config = json.load(open("config.json"))

class Primer_set (object):
    def __init__(self, forward_primer, reverse_primer, id):
        self.id = id
        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer

    def tm_dif(self):
        tm_difference = abs(self.forward_primer.primer_tm() - self.reverse_primer.primer_tm())
        return tm_difference

    def set_print(self):
        print("Forward id: %s - Reverse id:%s" %(self.forward_primer.id , self.reverse_primer.id))

    def get_amplicon_length(self):
        return (self.reverse_primer.start_index+self.reverse_primer.length) - self.forward_primer.start_index

    def write_to_file(self, file_name):
        write_file = open(file_name, "a")
        write_file.write(
            "----------------------------------------------------------------------------------\n" +
            "set ID: {}\n".format(self.id) +
            "FORWARD:\n" +
            "\tid:%s, Sequence:%s, length: %s\n\tTm:%.2f, GC:%s, Palindrome:%s, start index:%s\n" % (
            self.forward_primer.id, self.forward_primer.sequence, self.forward_primer.length, self.forward_primer.primer_tm(),
            self.forward_primer.precent_gc(), self.forward_primer.palindrome_length, self.forward_primer.start_index) +
            "REVERSE:\n" +
            "\tid:%s, Sequence:%s, length: %s\n\tTm:%.2f, GC:%s, Palindrome:%s, start index:%s\n" % (
                self.reverse_primer.id, self.reverse_primer.sequence, self.reverse_primer.length, self.reverse_primer.primer_tm(),
            self.reverse_primer.precent_gc(), self.forward_primer.palindrome_length, self.reverse_primer.start_index) +
            "tm difference between F & R: %.2f \n" % (self.tm_dif()) +
            "Amplicon length: %s \n" % (self.get_amplicon_length()) +
            "Score: %s \n" % (self.get_primers_set_score()))
        write_file.close()

# =========================== get primers set score ======================================================================
    def get_primers_set_score(self):
        #forward_scr = self.forward_primer.get_primer_score()
        #reverse_scr = self.reverse_primer.get_primer_score()
        #amplicon_scr = 100*self.get_amplicon_length_score()
        #return int(0.45*forward_scr + 0.45*reverse_scr + 0.1*amplicon_scr)
        score = ScorePredictor.logreg.predict_proba([ScorePredictor.normal_set(self, config)])
        return int(score[0][1]*100)

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


from PythonServer.Primer import Primer

class Primer_set (object):
    def __init__(self, forward_primer, reverse_primer):
        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer
    def tm_dif(self):
        tm_difference = abs(self.forward_primer.primer_tm() - self.reverse_primer.primer_tm())
        return tm_difference

    def set_print(self):
        print "Forward id: %s - Reverse id:%s" %(self.forward_primer.id , self.reverse_primer.id)

    def write_to_file(self,file_name):
        write_file = open(file_name, "a")
        write_file.write("----------------------------------------------------------------------------------------------\n" +
                         "FORWARD :\n" +
                         "Sequence: %s length: %s start index: %s \n"%(self.forward_primer.sequnce, self.forward_primer.length, self.forward_primer.start_index) +
                         "REVERSE:\n" +
                         "Sequence: %s length: %s start index: %s \n"%(self.reverse_primer.sequnce, self.reverse_primer.length, self.reverse_primer.start_index) +
                         "tm difference between F & R: %s \n" % (self.tm_dif()))

        write_file.close()

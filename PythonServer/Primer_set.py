from PythonServer.Primer import Primer

class Primer_set (object):
    def __init__(self, forward_primer, reverse_primer):
        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer

    def tm_dif(self):
        a = abs(self.forward_primer.primer_tm() - self.reverse_primer.primer_tm())
        return a

    def set_print(self):
        print "Forward id: %s - Reverse id:%s" %(self.forward_primer.id , self.reverse_primer.id)
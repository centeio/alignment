from .matrices import *

matrix_dict =   { 1: BLOSUM62, 2: DNAfull, 3: PAM250}


class Parameters:
    def __init__(self, gapext=-0.5, gapopen=-10, matrix=1, stype='protein'):
        self.gapext = gapext
        self.gapopen = gapopen
        self.matrix = matrix_dict[matrix]
        self.stype = stype

    def score(self, a, b):
        assert len(a) == len(b) == 1 #assures that is just one letter
        return self.matrix[a][b]

    def __str__(self):
        return "matrix = {}\nmatrix = {}\nstype = {}\ngapopen = {}\ngapext = {}".format(
            self.matrix, self.stype, self.gapopen, self.gapext
            )


class Alignment():
    def __init__(self, seq1, seq2):
        self.seq1 = seq1
        self.seq2 = seq2
        self.matches = "not calculated yet"
        self.mismatches = "not calculated yet"
        self.gaps = "not calculated yet"
    
    def calculate_mat_mis_gaps(self):
        self.gaps = 0
        self.matches = 0
        self.mismatches = 0
        
        for i in range(0,len(self.seq1)):
            if(self.seq1[i] == self.seq2[i]):
                self.matches = self.matches + 1
            elif(self.seq1[i] != '-' and self.seq2[i] != '-'):
                self.mismatches = self.mismatches + 1
            else:
                self.gaps = self.gaps + 1
                
    
    def __str__(self):
        return "[ALIGNMENT] \n\
        SEQ1: {} \n\
        SEQ2: {} \n\
        #gaps: {} \t #matches: {} \t #mismatches: {}".format(self.seq1,self.seq2,self.gaps,self.matches,self.mismatches)
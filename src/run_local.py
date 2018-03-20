#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from util.matrices import *
from util.alignment import Alignment
from util import io

#INF = sys.maxsize
INF = 1000
DEBUG = True

####################
#### PARAMETERS ####
####################

matrix_dict =   { 'BLOSUM62': BLOSUM62, 'DNAFULL': DNAfull, 'PAM250': PAM250}

class Parameters:
    def __init__(self, gap=-0.5, gapopen=-10, matrix='BLOSUM62', stype='protein'):
        self.gap = gap
        self.gapopen = gapopen
        self.matrix = matrix_dict[matrix]
        self.stype = stype

    def score(self, a, b):
        assert len(a) == len(b) == 1 #assures that is just one letter
        return self.matrix[a][b]

    def __str__(self):
        return "matrix = {}\nmatrix = {}\nstype = {}\ngap = {}\ngapext = {}".format(
            self.matrix, self.stype, self.gap, self.gapext
            )

##########################
#### GLOBAL ALIGNMENT ####
##########################

def global_align(seq1, seq2, Parameters=Parameters()):
    M = create_matrix(len(seq1), len(seq2))
    Ix = create_matrix(len(seq1), len(seq2))
    Iy = create_matrix(len(seq1), len(seq2))
    P = create_matrix(len(seq1), len(seq2))
    
    #initializes matrices
    for i in range(0, len(seq1)+1):
        for j in range(0, len(seq2)+1):
            if i==0 and j==0:
                Iy[i][j] = Parameters.gapopen + (Parameters.gap*j)
                Ix[i][j] = Parameters.gapopen + (Parameters.gap*i)
            elif i==0:
                Iy[i][j] = Parameters.gapopen + (Parameters.gap*j)
                Ix[i][j] = -INF
                M[i][j] = - INF
            elif j==0:
                Ix[i][j] = Parameters.gapopen + (Parameters.gap*i)
                Iy[i][j] = -INF
                M[i][j] = - INF

    M[0][0] = 0



    #Affine Score
    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):

                score = Parameters.score(seq1[i-1],seq2[j-1])
                diag = M[i-1][j-1] + score
                esq = Ix[i-1][j-1] + score
                cim = Iy[i-1][j-1] + score

                M[i][j] = max (
                    esq,
                    diag,
                    cim 
                )

                Ix[i][j] = max (
                    M[i-1][j] + Parameters.gapopen + Parameters.gap,
                    Ix[i-1][j] + Parameters.gap
                )

                Iy[i][j] = max (
                    M[i][j-1] + Parameters.gapopen + Parameters.gap,
                    Iy[i][j-1] + Parameters.gap
                )

    if (DEBUG):
        print("\nIx")
        print_matrix("_"+seq1,"_"+seq2,Ix)
        print("\nIy")
        print_matrix("_"+seq1,"_"+seq2,Iy)
        print("\nP")
        print_matrix("_"+seq1,"_"+seq2,P)
        print("\nM")
        print_matrix("_"+seq1,"_"+seq2,M)


    return M, Ix, Iy

def traceback(M, Ix, Iy ,seq1,seq2,Parameters):

    # As duas sequencias
    alignedseq1 = ""
    alignedseq2 = ""

    # Os ultimos indices
    i = len(seq1)
    j = len(seq2)

    #Matriz
    matrix = "M"
    
    print(seq1)
    print(seq2)

    while ((i is not 0) or (j is not 0)):
        print(i,j, matrix)
        if matrix == "M":
            score = Parameters.score(seq1[i-1],seq2[j-1])

            diag = M[i-1][j-1] + score
            esq = Ix[i-1][j-1] + score
            cim = Iy[i-1][j-1] + score

            if M[i][j] == esq:
                matrix = "Ix"
                alignedseq2 = alignedseq2 + "-"
                alignedseq1 = alignedseq1 + seq1[i-1]
                
            elif M[i][j] == diag:
                matrix = "M"
                alignedseq2 = alignedseq2 + seq2[j-1]
                alignedseq1 = alignedseq1 + seq1[i-1]

            elif M[i][j] == cim:
                matrix = "Iy"
                alignedseq1 = alignedseq1 + '-'
                alignedseq2 = alignedseq2 + seq2[j-1]

            i = i-1
            j = j-1 
        
        elif matrix == "Ix":
            if Ix[i][j] == M[i-1][j] + Parameters.gapopen + Parameters.gap:
                matrix = "M"
                alignedseq2 = alignedseq2 + seq2[j-1]
                alignedseq1 = alignedseq1 + seq1[i-1]
                
            elif Ix[i][j] == Ix[i-1][j] + Parameters.gap:
                matrix = "Ix"
                alignedseq2 = alignedseq2 + "-"
                alignedseq1 = alignedseq1 + seq1[i-1]
                
            i = i-1
        
        elif matrix == "Iy":
            if Iy[i][j] == M[i][j-1] + Parameters.gapopen + Parameters.gap:
                matrix = "M"
                alignedseq2 = alignedseq2 + seq2[j-1]
                alignedseq1 = alignedseq1 + seq1[i-1]
                
            elif Iy[i][j] == Iy[i][j-1] + Parameters.gap:
                matrix = "Iy"
                alignedseq1 = alignedseq1 + '-'
                alignedseq2 = alignedseq2 + seq2[j-1]
        
            j = j-1

        print(alignedseq1)
        print(alignedseq2)
        print(".")






    #Revertendo a String        
    alignedseq1 = alignedseq1[::-1]
    alignedseq2 = alignedseq2[::-1]

    return alignedseq1, alignedseq2




if __name__ == '__main__':
    #Definicoes dos parametros
    par = Parameters(gapopen=-10,gap=-0.5,matrix='BLOSUM62',stype='protein')

    #Sequencias
    seq1=io.read_fasta(io.read_file("../inputs/dummy1.fasta"))
    seq2 =io.read_fasta(io.read_file("../inputs/dummy2.fasta"))
    
    #Matriz de apontadores
    M, Ix, Iy = global_align(seq1, seq2, par)

    #Sequencias alinhadas
    alignedseq1, alignedseq2 = traceback(M, Ix, Iy, seq1, seq2, par)

    result = Alignment(alignedseq1,alignedseq2)
    result.calculate_mat_mis_gaps()
    print(str(result))
    io.write_file("../outputs/local_output.txt",str(result))
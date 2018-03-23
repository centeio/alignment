#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from util.matrices import *
from util.alignment import Alignment
from util import io

INF = sys.maxsize
#INF = 1000
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
#### LOCAL ALIGNMENT ####
##########################

def local_align(seq1, seq2, Parameters=Parameters()):
    M = create_matrix(len(seq1), len(seq2))
    Ix = create_matrix(len(seq1), len(seq2))
    Iy = create_matrix(len(seq1), len(seq2))
  
    #initializes matrices
    for i in range(0, len(seq1)+1):
        for j in range(0, len(seq2)+1):
            if i==0 and j==0:
                Iy[i][j] = Parameters.gapopen + (Parameters.gap*j)
                Ix[i][j] = Parameters.gapopen + (Parameters.gap*i)
                M[i][j] = 0
            elif i==0:
                Iy[i][j] = Parameters.gapopen + (Parameters.gap*j)
                Ix[i][j] = - INF
                M[i][j] = - INF
            elif j==0:
                Ix[i][j] = Parameters.gapopen + (Parameters.gap*i)
                Iy[i][j] = - INF
                M[i][j] = - INF




    #Max score
    maxmatrix = "M"
    max_score = 0
    li = 0
    lj = 0

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
                    cim,
                    0 
                )

                if M[i][j] >= max_score:
                    max_score = M[i][j]
                    li = i
                    lj = j
                    maxmatrix = "M"

                Ix[i][j] = max (
                    M[i-1][j] + Parameters.gapopen + Parameters.gap,
                    Ix[i-1][j] + Parameters.gap
                )

                if Ix[i][j] >= max_score:
                    max_score = Ix[i][j]
                    li = i
                    lj = j
                    maxmatrix = "Ix"

                Iy[i][j] = max (
                    M[i][j-1] + Parameters.gapopen + Parameters.gap,
                    Iy[i][j-1] + Parameters.gap
                )

                if Iy[i][j] >= max_score:
                    max_score = Iy[i][j]
                    li = i
                    lj = j
                    maxmatrix = "Iy"

        
    if (DEBUG):
        print("\nIx")
        print_matrix("_"+seq1,"_"+seq2,Ix)
        print("\nIy")
        print_matrix("_"+seq1,"_"+seq2,Iy)
        print("\nM")
        print_matrix("_"+seq1,"_"+seq2,M)


    return M, Ix, Iy, maxmatrix, li, lj

def traceback(M, Ix, Iy , seq1, seq2, Parameters, maxmatrix, li, lj):

    # As duas sequencias
    alignedseq1 = ""
    alignedseq2 = ""

    # Os indices do max score
    i = li
    j = lj


    #Matriz
    matrix = maxmatrix

    while True:
        print(matrix, i, j)
        if matrix == "M":
            alignedseq2 = alignedseq2 + seq2[j-1]
            alignedseq1 = alignedseq1 + seq1[i-1]

            if(M[i][j] == 0):
                break;

            score = Parameters.score(seq1[i-1],seq2[j-1])

            diag = M[i-1][j-1] + score
            esq = Ix[i-1][j-1] + score
            cim = Iy[i-1][j-1] + score

            if M[i][j] == esq:
                matrix = "Ix"
                
            elif M[i][j] == diag:
                matrix = "M"

            elif M[i][j] == cim:
                matrix = "Iy"
            
            i = i-1
            j = j-1 
        
        elif matrix == "Ix":
            alignedseq2 = alignedseq2 + "-"
            alignedseq1 = alignedseq1 + seq1[i-1]

            if(Ix[i][j] == 0):
                break;

            if Ix[i][j] == M[i-1][j] + Parameters.gapopen + Parameters.gap:
                matrix = "M"
                
            elif Ix[i][j] == Ix[i-1][j] + Parameters.gap:
                matrix = "Ix"

                
            i = i-1
        
        elif matrix == "Iy":
            alignedseq1 = alignedseq1 + '-'
            alignedseq2 = alignedseq2 + seq2[j-1]

            if(Iy[i][j] == 0):
                break;

            if Iy[i][j] == M[i][j-1] + Parameters.gapopen + Parameters.gap:
                matrix = "M"
                
            elif Iy[i][j] == Iy[i][j-1] + Parameters.gap:
                matrix = "Iy"
        
            j = j-1



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
    M, Ix, Iy, maxmatrix, li, lj = local_align(seq1, seq2, par)

    #Sequencias alinhadas
    alignedseq1, alignedseq2 = traceback(M, Ix, Iy, seq1, seq2, par, maxmatrix, li, lj)

    result = Alignment(alignedseq1,alignedseq2)
    result.calculate_mat_mis_gaps()
    print(str(result))
    io.write_file("../outputs/locally_local_output.txt",str(result))
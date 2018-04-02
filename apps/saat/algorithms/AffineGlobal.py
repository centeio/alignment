#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from .util.alignment import Alignment, Parameters
from .util import io
from .util.matrices import *
 
INF = sys.maxsize
#INF = 1000
DEBUG = False

##########################
#### GLOBAL ALIGNMENT ####
##########################


class AffineGlobal():
    def __init__(self,multiples=False):
        self.multiples = multiples

    def affine_global_align(self,seq1, seq2, Parameters=Parameters()):
        M = create_matrix(len(seq1), len(seq2))
        Ix = create_matrix(len(seq1), len(seq2))
        Iy = create_matrix(len(seq1), len(seq2))
        P = create_matrix(len(seq1), len(seq2))
        
        #initializes matrices
        for i in range(0, len(seq1)+1):
            for j in range(0, len(seq2)+1):
                if i==0:
                    Iy[i][j] = Parameters.gapopen + (Parameters.gapext*j)
                    Ix[i][j] = -INF
                    M[i][j] = - INF
                if j==0:
                    Ix[i][j] = Parameters.gapopen + (Parameters.gapext*i)
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
                        M[i-1][j] + Parameters.gapopen + Parameters.gapext,
                        Ix[i-1][j] + Parameters.gapext
                    )

                    Iy[i][j] = max (
                        M[i][j-1] + Parameters.gapopen + Parameters.gapext,
                        Iy[i][j-1] + Parameters.gapext
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
        

    def traceback(self, M, Ix, Iy, seq1, seq2, Parameters):

        # As duas sequencias
        alignedseq1 = ""
        alignedseq2 = ""

        # Os ultimos indices
        i = len(seq1)
        j = len(seq2)

        #Matriz
        if(M[i][j] >= Ix[i][j] and M[i][j] >= Iy[i][j]):
            matrix = "M"
        elif(Ix[i][j] > Iy[i][j]):
            matrix = "Ix"
        else:
            matrix = "Iy"


        while ((i is not 0) or (j is not 0)):
            if matrix == "M":
                alignedseq2 = alignedseq2 + seq2[j-1]
                alignedseq1 = alignedseq1 + seq1[i-1]

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

                if Ix[i][j] == M[i-1][j] + Parameters.gapopen + Parameters.gapext:
                    matrix = "M"
                    
                elif Ix[i][j] == Ix[i-1][j] + Parameters.gapext:
                    matrix = "Ix"

                    
                i = i-1
            
            elif matrix == "Iy":
                alignedseq1 = alignedseq1 + '-'
                alignedseq2 = alignedseq2 + seq2[j-1]

                if Iy[i][j] == M[i][j-1] + Parameters.gapopen + Parameters.gapext:
                    matrix = "M"
                    
                elif Iy[i][j] == Iy[i][j-1] + Parameters.gapext:
                    matrix = "Iy"
            
                j = j-1



        #Revertendo a String        
        alignedseq1 = alignedseq1[::-1]
        alignedseq2 = alignedseq2[::-1]

        return alignedseq1, alignedseq2

    def run(self, seq1, seq2, par):
        seq1 = io.read_fasta(seq1)
        seq2 = io.read_fasta(seq2)

        M, Ix, Iy = self.affine_global_align(seq1, seq2, par)

        alignedseq1, alignedseq2 = self.traceback(M, Ix, Iy, seq1, seq2, par)
        result = Alignment(alignedseq1,alignedseq2)
        result.calculate_mat_mis_gaps()

        return result


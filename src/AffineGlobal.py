#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from util.matrices import *
from util.alignment import Alignment
from util.alignment import Parameters
from util import io

INF = sys.maxsize
#INF = 1000
DEBUG = False


#########################
#### GLOBAL ALIGNMENT ####
#########################

def global_align(seq1, seq2, Parameters=Parameters()):
    M = create_matrix(len(seq1), len(seq2))
    Ix = create_matrix(len(seq1), len(seq2))
    Iy = create_matrix(len(seq1), len(seq2))
  
    #initializes matrices
    for i in range(0, len(seq1)+1):
        for j in range(0, len(seq2)+1):
            if i==0 and j==0:
                M[i][j] = 0
            elif i==0:
                Iy[i][j] = 0
                Ix[i][j] = - INF
                M[i][j] = - INF
            elif j==0:
                Ix[i][j] = 0
                Iy[i][j] = -INF
                M[i][j] = - INF


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
        print("\nM")
        print_matrix("_"+seq1,"_"+seq2,M)


    return M, Ix, Iy

def traceback_left(M, Ix, Iy ,seq1,seq2,Parameters):

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

            if Ix[i][j] == M[i-1][j] + Parameters.gapopen + Parameters.gap:
                matrix = "M"
                
            elif Ix[i][j] == Ix[i-1][j] + Parameters.gap:
                matrix = "Ix"

                
            i = i-1
        
        elif matrix == "Iy":
            alignedseq1 = alignedseq1 + '-'
            alignedseq2 = alignedseq2 + seq2[j-1]

            if Iy[i][j] == M[i][j-1] + Parameters.gapopen + Parameters.gap:
                matrix = "M"
                
            elif Iy[i][j] == Iy[i][j-1] + Parameters.gap:
                matrix = "Iy"
        
            j = j-1



    #Revertendo a String        
    alignedseq1 = alignedseq1[::-1]
    alignedseq2 = alignedseq2[::-1]

    return alignedseq1, alignedseq2

def traceback_up(M, Ix, Iy ,seq1,seq2,Parameters):

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

            if M[i][j] == cim:
                matrix = "Iy"

            elif M[i][j] == diag:
                matrix = "M"

            elif M[i][j] == esq:
                matrix = "Ix"

            i = i-1
            j = j-1 
        
        elif matrix == "Ix":
            alignedseq2 = alignedseq2 + "-"
            alignedseq1 = alignedseq1 + seq1[i-1]
                
            if Ix[i][j] == Ix[i-1][j] + Parameters.gap:
                matrix = "Ix"

            elif Ix[i][j] == M[i-1][j] + Parameters.gapopen + Parameters.gap:
                matrix = "M"
                
            i = i-1
        
        elif matrix == "Iy":
            alignedseq1 = alignedseq1 + '-'
            alignedseq2 = alignedseq2 + seq2[j-1]
                
            if Iy[i][j] == Iy[i][j-1] + Parameters.gap:
                matrix = "Iy"

            elif Iy[i][j] == M[i][j-1] + Parameters.gapopen + Parameters.gap:
                matrix = "M"

        
            j = j-1



    #Revertendo a String        
    alignedseq1 = alignedseq1[::-1]
    alignedseq2 = alignedseq2[::-1]

    return alignedseq1, alignedseq2




if __name__ == '__main__':
    #Definicoes dos parametros
    par = Parameters(gapopen=-10,gap=-0.5,matrix='BLOSUM62',stype='protein')

    #Sequencias
    seq1=io.read_fasta(io.read_file("../inputs/default1.fasta"))
    seq2=io.read_fasta(io.read_file("../inputs/default2.fasta"))
    
    #Matriz de apontadores
    M, Ix, Iy = global_align(seq1, seq2, par)

    #Sequencias alinhadas
    leftalignedseq1, leftalignedseq2 = traceback_left(M, Ix, Iy, seq1, seq2, par)
    upalignedseq1, upalignedseq2 = traceback_up(M, Ix, Iy, seq1, seq2, par)

    result = Alignment(leftalignedseq1,leftalignedseq2, "LEFT")
    result.calculate_mat_mis_gaps()
    io.write_file("../outputs/locally_global_affine_output.txt",str(result))

    result = Alignment(upalignedseq1,upalignedseq2, "UP")
    result.calculate_mat_mis_gaps()
    io.append_file("../outputs/locally_global_affine_output.txt",str(result))
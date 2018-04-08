#!/usr/bin/env python
# -*- coding: utf-8 -*-

from util.alignment import Alignment
from util.alignment import Parameters
from util import io
from util.matrices import create_matrix

##########################
#### GLOBAL ALIGNMENT ####
##########################

def global_align(seq1, seq2, Parameters=Parameters()):
    M = create_matrix(len(seq1), len(seq2))

    for i in range(0, len(seq1)+1):
        for j in range(0, len(seq2)+1):
            if i == 0 and j == 0:
                M[i][j] = 0
            elif i == 0:
                M[i][j] = M[i][j - 1] + Parameters.gap
            elif j == 0:
                M[i][j] = M[i - 1][j] + Parameters.gap
            else:
                vj = M[i][j-1] + Parameters.gap
                vi = M[i-1][j] + Parameters.gap
                vd = M[i-1][j-1] + Parameters.score(seq1[i-1], seq2[j-1])

                M[i][j] = max(vj, vi, vd)


    return M

def traceback(M,seq1,seq2,Parameters):
    """
    Esta funcao encontra o melhor aliamente seguindo o algoritmo a seguir:

        0- verificar se a celula corrente esta na borda
        1- obter os n vizinhos
        2- obter o valor da celula corrente
        3- Usar a formula
               a) M[i][j-1] + Parameters.gap == atual ?
               b) M[i-1][j-1] + Parameters.score(seq1[i-1], seq2[j-1]) == atual ?
               c) M[i-1][j] + Parameters.gap == atual ?
        4 - Criar strings alignedseq1 e alignedseq2
        5 - caso a) (esquerda)
             alignedseq2 + "_"
             alignedseq1 + seq1[j-1]   #"_"+seq1
             j = j-1
        
             caso b) (diagonal)
              alignedseq1 + seq1[j-1]
              alignedseq2 + seq2[i-1]
              j = j-1
              i = i-1
        
             caso c) (cima)
              alignedseq1 + "_"
              alignedseq2 + seq2[i-1]
              i = i-1
        6 - reverter a string

    OBS: J É LINHA e I É COLUNA
    """

    # As duas sequencias
    alignedseq1 = ""
    alignedseq2 = ""

    # Os ultimos indices
    i = len(seq1)
    j = len(seq2)
    
    while ((i is not 0) or (j is not 0)):
        #Primeira linha
        if j==0:
            alignedseq2 = alignedseq2 + "-"
            alignedseq1 = alignedseq1 + seq1[i-1]
            i = i-1
            continue
        
        #Primeira coluna
        if i==0:
            alignedseq1 = alignedseq1 + "-"
            alignedseq2 = alignedseq2 + seq2[j-1]
            j = j-1
            continue
       
        #Esquerda
        if (M[i][j] == M[i-1][j] + Parameters.gap):
            alignedseq2 = alignedseq2 + "-"
            alignedseq1 = alignedseq1 + seq1[i-1]
            i = i-1
            continue
            
        #Diagonal
        if (M[i][j] == M[i-1][j-1] + Parameters.score(seq1[i-1], seq2[j-1])):
            alignedseq2 = alignedseq2 + seq2[j-1]
            alignedseq1 = alignedseq1 + seq1[i-1]
            j = j-1
            i = i-1
            continue

        #Cima
        if (M[i][j] == M[i][j-1] + Parameters.gap):
            alignedseq1 = alignedseq1 + '-'
            alignedseq2 = alignedseq2 + seq2[j-1]
            j=j-1
            continue
        

    #Revertendo a String        
    alignedseq1 = alignedseq1[::-1]
    alignedseq2 = alignedseq2[::-1]

    return alignedseq1, alignedseq2




if __name__ == '__main__':
    #Definicoes dos parametros
    par = Parameters(gap=-10,matrix='BLOSUM62',stype='protein')

    
    seq1=io.read_fasta(io.read_file("../inputs/NM_002688.fasta"))
    seq2=io.read_fasta(io.read_file("../inputs/XM_001166286.fasta"))
    
    matrix = global_align(seq1, seq2, par)

    alignedseq1, alignedseq2 = traceback(matrix, seq1, seq2, par)

    result = Alignment(alignedseq1,alignedseq2, "DEFAULT")
    result.calculate_mat_mis_gaps()

    io.write_file("../outputs/locally_global_linear_output.txt",str(result))
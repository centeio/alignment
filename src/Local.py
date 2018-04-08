#!/usr/bin/env python
# -*- coding: utf-8 -*-


from util.alignment import Alignment
from util.alignment import Parameters
from util import io
from util.matrices import create_matrix
########################
#### LOCAL ALIGNMENT ####
########################

def local_align(seq1, seq2, Parameters=Parameters()):
    M = create_matrix(len(seq1), len(seq2))

    score = 0
    li = 0
    lj = 0

    # fill in A in the right order
    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):

            # the local alignment recurrance rule:
            M[i][j] = max(
               M[i][j-1] + Parameters.gap,
               M[i-1][j] + Parameters.gap,
               M[i-1][j-1] + Parameters.score(seq1[i-1], seq2[j-1]),
               0
            )

            # track the cell with the largest score
            if M[i][j] >= score:
                score = M[i][j]
                li = i
                lj = j


    return M, li, lj



def traceback(M,seq1,seq2,i,j,Parameters):
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

    while M[i][j] != 0:


        #Primeira linha ou esquerda
        if j==0 or (M[i][j] == M[i-1][j] + Parameters.gap):
            #if M[i-1][j] is not 0:
            alignedseq2 = alignedseq2 + "-"
            alignedseq1 = alignedseq1 + seq1[i-1]
            i = i-1
            continue
        
        #Primeira coluna
        if i==0:
            #if M[i][j-1] is not 0:
            alignedseq1 = alignedseq1 + "-"
            alignedseq2 = alignedseq2 + seq2[j-1]
            j = j-1
            continue

    
        #Diagonal
        if M[i][j] == M[i-1][j-1] + Parameters.score(seq1[i-1], seq2[j-1]):
            #if M[i-1][j-1] is not 0:
            alignedseq2 = alignedseq2 + seq2[j-1]
            alignedseq1 = alignedseq1 + seq1[i-1]
            j = j-1
            i = i-1
            continue

        #Cima
        if M[i][j] == M[i][j-1] + Parameters.gap:
            #if M[i][j-1] is not 0:
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
    par = Parameters(gap=-10,matrix='BLOSUM62',stype='dna')

    
    seq1=io.read_fasta(io.read_file("../inputs/default1.fasta"))
    seq2=io.read_fasta(io.read_file("../inputs/default2.fasta"))
    
    matrix, i, j = local_align(seq1, seq2, par)

    alignedseq1, alignedseq2 = traceback(matrix, seq1, seq2, i, j, par)

    result = Alignment(alignedseq1,alignedseq2, "DEFAULT")
    result.calculate_mat_mis_gaps()

    io.write_file("../outputs/locally_local_linear_output.txt",str(result))
#!/usr/bin/env python
# -*- coding: utf-8 -*-

########################
#### UTIL FUNCTIONS ####
########################

def create_matrix(n,m):
    """
    Equivalent of
    int[][]M = [n][m];
    for (int i; i<n; i++){
        for (int j; j<m; j++){
            M[i][j] = 0;
        }
    }
    """
    return [[0]*(m+1) for i in range(n+1)]

####################
#### PARAMETERS ####
####################

matrix_dict =   { 'BLOSUM62': BLOSUM62, 'DNAFULL': DNAFULL, 'PAM': PAM}

class Parameters:
    def __init__(self, gap=-2, matrix='BLOSUM62', stype='protein'):
        self.gap = gap
        self.matrix = matrix_dict[matrix]
        self.stype = stype

    def score(self, a, b):
        assert len(a) == len(b) == 1 #assures that is just one letter
        return matrix[a][b]

    def __str__(self):
        return "matrix = {}\nmatrix = {}\nstype = {}\ngap = {}\n".format(
            self.matrix, self.stype, self.gap
            )

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

    print_matrix("_"+seq1,"_"+seq2,M)

    return M

def align(M,seq1,seq2,Parameters):
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
    
    while i!=0 and j!=0:

        #Primeira linha ou esquerda
        if j==0 or (M[i][j] == M[i-1][j] + Parameters.gap):
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

def print_matrix(seq1, seq2, M):
    """
    This function prints the matrix along with sequences
    Example:
              A     A     A
        C     0     0     0
        C     0     0     0
        C     0     0     0
    """
    # print the empty space
    # print the top row
    print(" ", end=" ")
    for c in seq1:
        print(c, end=" ")
    print("")

    # print the remaining
    for j in range(len(M[0])):
        print(seq2[j], end=" ")
        for i in range(len(M)):
            print(M[i][j], end=" ")
        print("")



if __name__ == '__main__':
    #Definicoes dos parametros
    par = Parameters(gap=-2,matrix='BLOSUM62',stype='protein')

    #Sequencias
    #>tr|B7Z8R9|B7Z8R9_HUMAN cDNA FLJ57474, highly similar to Homo sapiens plasticity
    # related gene 3 (PRG-3), transcript variant 1, mRNA OS=Homo sapiens PE=2 SV=1
    seq1 = "MAGTVLLAYYFECTDTFQVHIQGFFCQDGDLMKPYPGTEEESFITPLVLYCVLAATPTAI\
    IFIGEISMYFIKSTRESLIAQEKTILTGECCYLNPLLRRIIRFTGVFAFGLFATDIFVNA\
    GQVVTGHLTPYFLTVCKPNYTSADCQAHHQFINNGNICTGDLEVIEKARRSFPSKHAALS\
    IYSALYATMYITSTIKTKSSRLAKPVLCLGTLCTAFLTGLNRVSEYRNHCSDVIAGFILG\
    TAVALFLGMCVVHNFKGTQGSPSKPKPEDPRGVPLMAFPRIESPLETLSAQNHSASMTEVT"


    #>sp|Q9Y2Y8|PRG3_HUMAN Proteoglycan 3 OS=Homo sapiens GN=PRG3 PE=1 SV=2
    seq2 = "MQCLLLLPFLLLGTVSALHLENDAPHLESLETQADLGQDLDSSKEQERDLALTEEVIQAE\
    GEEVKASACQDNFEDEEAMESDPAALDKDFQCPREEDIVEVQGSPRCKICRYLLVRTPKT\
    FAEAQNVCSRCYGGNLVSIHDFNFNYRIQCCTSTVNQAQVWIGGNLRGWFLWKRFCWTDG\
    SHWNFAYWSPGQPGNGQGSCVALCTKGGYWRRAQCDKQLPFVCSF"

    #Matriz
    matrix = global_align(seq1, seq2, par)
    #Sequencias alinhadas
    alignedseq1, alignedseq2 = align(matrix, seq1, seq2, par)

    print(alignedseq1)
    print(alignedseq2)
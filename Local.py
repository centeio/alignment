"""
Python code for local alignment

Usage:
    local_align("AAA","AAC",Parameters(gap=-2,mismatch=-1,match=1))
"""

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
    print ("{}".format(" "),end=" ")
    # print the top row
    for c in seq1:
        print (c,end=" ")
    
    print("",end="\n")

    # print the remaining
    for j in range(len(M[0])):
        print (seq2[j],end=" ")
        for i in range(len(M)):
            print (M[i][j], end=" ")
        print("",end="\n")

####################
#### PARAMETERS ####
####################

class Parameters:
    def __init__(self, gap=-2,mismatch=-1,match=1):
        self.gap = gap
        self.mismatch = mismatch
        self.match= match

    def score(self, a, b):
        assert len(a) == len(b) == 1 #assures that is just one letter

        if a==b:
            return self.match
        else:
            return self.mismatch

    def __str__(self):
        return "match = {}\nmismatch = {}\ngap = {}\n".format(
                self.match, self.mismatch, self.gap
        )

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

    print ("Parameters\n{}".format(str(Parameters)))
    print ("Score = {}".format(score))
    print ("Location in the matrix = {},{}".format(li,lj))
    # just putting an * in the location
    #M[li][lj] = "*"+str(M[li][lj])
    print ("Matrix =")
    print_matrix("_"+seq1, "_"+seq2, M)

    return M, li, lj



def align(M,seq1,seq2,i,j,Parameters):
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


    par = Parameters()

    #>tr|B7Z8R9|B7Z8R9_HUMAN cDNA FLJ57474, highly similar to Homo sapiens plasticity
    # related gene 3 (PRG-3), transcript variant 1, mRNA OS=Homo sapiens PE=2 SV=1
    a = "MAGTVLLAYYFECTDTFQVHIQGFFCQDGDLMKPYPGTEEESFITPLVLYCVLAATPTAI\
    IFIGEISMYFIKSTRESLIAQEKTILTGECCYLNPLLRRIIRFTGVFAFGLFATDIFVNA\
    GQVVTGHLTPYFLTVCKPNYTSADCQAHHQFINNGNICTGDLEVIEKARRSFPSKHAALS\
    IYSALYATMYITSTIKTKSSRLAKPVLCLGTLCTAFLTGLNRVSEYRNHCSDVIAGFILG\
    TAVALFLGMCVVHNFKGTQGSPSKPKPEDPRGVPLMAFPRIESPLETLSAQNHSASMTEVT"
    
    
    #>sp|Q9Y2Y8|PRG3_HUMAN Proteoglycan 3 OS=Homo sapiens GN=PRG3 PE=1 SV=2
    b = "MQCLLLLPFLLLGTVSALHLENDAPHLESLETQADLGQDLDSSKEQERDLALTEEVIQAE\
    GEEVKASACQDNFEDEEAMESDPAALDKDFQCPREEDIVEVQGSPRCKICRYLLVRTPKT\
    FAEAQNVCSRCYGGNLVSIHDFNFNYRIQCCTSTVNQAQVWIGGNLRGWFLWKRFCWTDG\
    SHWNFAYWSPGQPGNGQGSCVALCTKGGYWRRAQCDKQLPFVCSF"

    matrix, i, j = local_align(a,b,par)
    alignedseq1, alignedseq2 = align(matrix,a,b,i,j,par)
    print(alignedseq1)
    print(alignedseq2)
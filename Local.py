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
    return [[0]*m for i in xrange(n)]

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
    print "{:^5}".format(" "),
    # print the top row
    for c in seq1:
        print "{:^5}".format(c),
    print

    # print the remaining
    for j in xrange(len(M[0])):
        print "{:>3}".format(seq2[j]),
        for i in xrange(len(M)):
            print "{:>5}".format(M[i][j]),
        print

####################
#### PARAMETERS ####
####################

class Parameters:
    def __init__(self, gap=-2,mismatch=-1,match=1):
        self.gap = gap
        self.mismatch = mismatch
        self.match=1

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
    location = (0,0)

    # fill in A in the right order
    for i in xrange(1, len(seq1)):
        for j in xrange(1, len(seq2)):

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
                location = (i,j)

    print "Parameters\n{}".format(Parameters)
    print "Score = {}".format(score)
    print "Location in the matrix = {}".format(location)
    # just putting an * in the location
    M[location[0]][location[1]] = "*"+str(M[location[0]][location[1]])
    print "Matrix ="
    print_matrix(seq1, seq2, M)


    """
    TODO: get the sequences
    alignedseq1, alignedseq2 = reverse_search(seq1, seq2,location,M)
    """

    
    return score


if __name__ == '__main__':

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
    
    local_align(a,b,Parameters())
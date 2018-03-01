
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


global_align("AACA", "AG", Parameters())

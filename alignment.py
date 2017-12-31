#This software is a free software. Thus, it is licensed under GNU General Public License.
#Python implementation to Smith-Waterman Algorithm for Homework 1 of Bioinformatics class.
#Forrest Bao, Sept. 26 <http://fsbao.net> <forrest.bao aT gmail.com>

# zeros() was origianlly from NumPy.
# This version is implemented by alevchuk 2011-04-10
def zeros(shape):
    retval = []
    for x in range(shape[0]):
        retval.append([])
        for y in range(shape[1]):
            retval[-1].append(0)
    return retval

match_award       = 8
mismatch_penalty  = -4
gap_penalty       = -5
extension_penalty = -2

def revcomp(seq):
    comp = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
    return(''.join([comp[q] for q in seq[::-1]]))

def match_score(alpha, beta):
    if alpha == beta or alpha == 'N' or beta == 'N':
        return match_award
    else:
        return mismatch_penalty

def water(seq1, seq2):
    m, n = len(seq1), len(seq2)  # length of two sequences
    
    # Generate DP table and traceback path pointer matrix
    score = zeros((m+1, n+1))      # the DP table
    pointer = zeros((m+1, n+1))    # to store the traceback path
    for i in range(1, m+1):
        pointer[i][0] = 1
    for j in range(1, n+1):
        pointer[0][j] = 2
    
    # Calculate DP table and mark pointers
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            score_diagonal = score[i-1][j-1] + match_score(seq1[i-1], seq2[j-1])
            if (i == m):
                score_left = score[i][j-1]
            elif (pointer[i][j-1] == 2):
                score_left = score[i][j-1] + extension_penalty
            else:
                score_left = score[i][j-1] + gap_penalty
            if (j == n):
                score_up = score[i-1][j]
            elif (pointer[i-1][j] == 1):
                score_up = score[i-1][j] + extension_penalty
            else:
                score_up = score[i-1][j] + gap_penalty
            
            score[i][j] = max(0,score_left, score_up, score_diagonal)
            if score[i][j] == 0:
                pointer[i][j] = 0 # 0 means end of the path
            if score[i][j] == score_up:
                pointer[i][j] = 1 # 1 means trace up
            if score[i][j] == score_left:
                pointer[i][j] = 2 # 2 means trace left
            if score[i][j] == score_diagonal:
                pointer[i][j] = 3 # 3 means trace diagonal
    
    align1, align2 = '', ''    # initial sequences
    
    i,j = m , n
    
    #traceback, follow pointers (note that this constructs the reverse string)
    while pointer[i][j] != 0:
        if pointer[i][j] == 3:
            align1 += seq1[i-1]
            align2 += seq2[j-1]
            i -= 1
            j -= 1
        elif pointer[i][j] == 2:
            align1 += '-'
            align2 += seq2[j-1]
            j -= 1
        elif pointer[i][j] == 1:
            align1 += seq1[i-1]
            align2 += '-'
            i -= 1
    
    return((score[m][n],align1[::-1],align2[::-1]))
    

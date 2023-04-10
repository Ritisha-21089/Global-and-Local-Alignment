import numpy as np

# Initializing the two sequences as arrays:
s1 = "GATGCGCAG"
s2 = "GGCAGTA"

# the scoring scheme is as follows: (same as question2)
Match = +2
Mismatch = -1
Gap = -3

# initializing the matrix :
n = len(s1)
m = len(s2)

# As per the Smith-Waterman algorithm, the initialization for the matrix is to assign 0 to first row & col
matrix = np.zeros((len(s1) + 1, len(s2) + 1))

for i in range(n):
    for j in range(m):
        # match_score : check if the letter matches and assign match/mismatch (replace with 0 if value is negative)
        if s1[i] == s2[j]:
            match_score = max(0, Match)
        else:
            match_score = max(0, Mismatch)

        # choosing the max value for each cell (ensuring value !=0), thereby generating the matrix :
        matrix[i + 1][j + 1] = max(max(0, matrix[i][j] + match_score), max(0, matrix[i][j + 1] + Gap),
                                   max(0, matrix[i + 1][j] + Gap))

# Now, the matrix has been filled, we need to perform the traceback after identifying the score (which will be the max
# value in the matrix):
max_score = 0

# Find the cell with the maximum score and its position (represented by [a][b])
for i in range(len(s1) + 1):
    for j in range(len(s2) + 1):
        if matrix[i][j] > max_score:
            max_score = matrix[i][j]
            a = i
            b = j

aligned_seqs = []
# There are more than one possibility of aligning, and we will be tracing all of those below :

def traceback(a, b, aligned_s1, aligned_s2):
    if matrix[a][b] == 0:
        aligned_seqs.append((aligned_s1, aligned_s2))
        return
    if a > 0 and b > 0 and matrix[a][b] == matrix[a - 1][b - 1] + (Match if s1[a - 1] == s2[b - 1] else Mismatch):
        traceback(a - 1, b - 1, s1[a - 1] + aligned_s1, s2[b - 1] + aligned_s2)
    if a > 0 and matrix[a][b] == max(0, matrix[a - 1][b] + Gap):
        traceback(a - 1, b, s1[a - 1] + aligned_s1, '-' + aligned_s2)
    if b > 0 and matrix[a][b] == max(0, matrix[a][b - 1] + Gap):
        traceback(a, b - 1, '-' + aligned_s1, s2[b - 1] + aligned_s2)


traceback(a, b, '', '')

print()
print('LOCAL ALLIGNMENT MATRIX : ')
print()
matrix_print = matrix.transpose()
for i in range(m + 1):
    print(matrix_print[i])
print()
print("SCORE : ", max_score)
print()
print("Alignment(s) for the score generated above (optimal alignments): ")
for i, seq in enumerate(aligned_seqs):
    print(f'Alignment {i + 1}:')
    print(seq[0])
    print(seq[1])
    print('\n')

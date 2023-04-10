import numpy as np

s1 = "GATGCGCAG"
s2 = "GGCAGTA"

s2='0'+s2
s1='0'+s1

scheme=[2, -3, -1]
#match mismatch gap

m= len(s1)
n= len(s2)

mat = np.zeros((n , m ))
# As per the Needleman-Wunsch algorithm, first row and column initialised with gap penalties
for i in range(1, n):
    mat[i][0] = mat[i-1][0]+ scheme[2]
for j in range(1, m):
    mat[0][j] = mat[0][j-1]+ scheme[2]

score=[0,0,0] #diag left right
for i in range(1, n):
    for j in range(1, m):
        if s2[i] == s1[j]:
           score[0]= scheme[0] + mat[i-1][j-1]
        else:
           score[0]= scheme[1] + mat[i-1][j-1]
            
        score[1]= mat[i][j-1] + scheme[2]
        score[2]= mat[i-1][j] + scheme[2]
        
        mat[i][j] = max(score)

#We will now trace back the possible paths taken during the alignment to identify all optimal alignments.

aligned_sequences = []

# Recursive function to trace all possible paths
def alignments(i, j, aligned_s2, aligned_s1):
    # Add the aligned sequences to the list
    if i == 0 and j == 0:
        aligned_sequences.append((aligned_s2, aligned_s1))
        return
    
    # Trace the diagonal path
    if i > 0 and j > 0 and mat[i][j] == mat[i - 1][j - 1] + (scheme[1], scheme[0])[s1[j] == s2[i]]:
        alignments(i - 1, j - 1, s2[i] + aligned_s2, s1[j] + aligned_s1)

    # Trace the gap in sequence 1 path
    if i > 0 and mat[i][j] == mat[i - 1][j] + scheme[2]:
        alignments(i - 1, j, s2[i] + aligned_s2, '-' + aligned_s1)

    # Trace the gap in sequence 2 path
    if j > 0 and mat[i][j] == mat[i][j - 1] + scheme[2]:
        alignments(i, j - 1, '-' + aligned_s2, s1[j] +  aligned_s1)

alignments(n-1, m-1, '', '')

#------------------------------------------------------------------------------------------------------------------
print()
print('GLOBAL ALLIGNMENT MATRIX: ')
print()
for e in list(s1):
    print(' ',e,'', end="")
print()

s3= (list(s2))
for i in range (len(s3)):
    print(s3[i], mat[i])

print()
print("SCORE: ", mat[n-1][m-1])
print()


print("Optimal Alignment(s): ")
count=0
print()
for seq in (aligned_sequences):
    print('Alignment ',count + 1,':')
    print(seq[1])
    print(seq[0])
    print('\n')
    count+=1

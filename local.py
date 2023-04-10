import numpy as np

s1 = "GATGCGCAG"
s2 = "GGCAGTA"

s2='0'+s2
s1='0'+s1

scheme=[2, -1, -3]
#match mismatch gap

m= len(s1)
n= len(s2)

mat = np.zeros((n , m ))
# As per the Needleman-Wunsch algorithm, first row and column initialised with gap penalties

score=[0,0,0,0] #diag left right
for i in range(1, n):
    for j in range(1, m):
        if s2[i] == s1[j]:
            score[0] = mat[i-1][j-1] + scheme[0]
        else:
            score[0] = mat[i-1][j-1] + scheme[1]

        score[1]=max(0, mat[i][j-1]+ scheme[2])
        score[2]=max(0, mat[i-1][j]+ scheme[2])
        
        mat[i][j] = max(score)

#We will now trace back the possible paths taken during the alignment to identify all optimal alignments.

aligned_sequences = []

# Recursive function to trace all possible paths
def alignments(i, j, aligned_s2, aligned_s1):
    if i == 0 and j == 0:
        aligned_sequences.append((aligned_s2, aligned_s1)) # Add the aligned sequences to the list
        return
    
    # Trace the diagonal path
    if i > 0 and j > 0 and mat[i][j] == mat[i - 1][j - 1] + (scheme[1], scheme[0])[s1[j] == s2[i]]: 
        alignments(i - 1, j - 1, s2[i] + aligned_s2, s1[j] + aligned_s1) 

    # Trace the gap in sequence 1 path
    if i > 0 and mat[i][j] == max(0, mat[i - 1][j] + scheme[2]):
        alignments(i - 1, j, s2[i] + aligned_s2, '-' + aligned_s1)

    # Trace the gap in sequence 2 path
    if j > 0 and mat[i][j] == max(0, mat[i][j - 1] + scheme[2]):
        alignments(i, j - 1, '-' + aligned_s2, s1[j] +  aligned_s1)

alignments(n-1, m-1, '', '')

#------------------------------------------------------------------------------------------------------------------
print()
print('LOCAL ALLIGNMENT MATRIX: ')
print()
for e in list(s1):
    print(' ',e, end="")
print()

s3= (list(s2))
max_score=0
for i in range (len(s3)):
    print(s3[i], mat[i])
    max_score= max( max(mat[i]), max_score) # Calculate the max score in the matrix

print()
print("SCORE: ", max_score)
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

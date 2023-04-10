# Global-and-Local-Alignment

The main difference between the Local and Global Alignment is the way they are scored and aligned. 
In Global Alignment, the entire length of both sequences is aligned to find the optimal alignment, whereas in Local Alignment, only the best matching sub-sequences are aligned.

In the given codes, the Local Alignment algorithm is implemented using the Needleman-Wunsch algorithm, whereas the Global Alignment algorithm is implemented using a modified version of the Needleman-Wunsch algorithm.

The major differences in the implementation are:

1.) Initialization of the first row and column with gap penalties: In Local Alignment, the first row and column are initialized with 0, whereas in Global Alignment, they are initialized with gap penalties.

2.) Scoring matrix values: In Local Alignment, negative scores are set to 0. In contrast, in Global Alignment, negative scores are not set to 0 but are included in the scoring matrix.

3.) Traceback paths: In Local Alignment, traceback starts from the cell with the maximum score in the matrix, whereas in Global Alignment, it starts from the last cell in the matrix.

4.) The scoring scheme: Both algorithms use a scoring scheme that assigns positive scores to matches and negative scores to mismatches and gaps, but the values for these scores differ in both algorithms.

Overall, the Local Alignment algorithm is used when searching for similarities between smaller sub-sequences, whereas the Global Alignment algorithm is used when comparing entire sequences.

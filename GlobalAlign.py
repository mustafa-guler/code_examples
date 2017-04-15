#Uses the standard Global Alignment dynmaic programming solution to align
#two protein sequences. Uses the BLOSUM62 scoring function and an indel
#penalty of -5.
#Does not taken into account longer gaps.
#
#How to Run:
#python GlobalAlign.py <input_file1> <input_file2>
#
#each input file should be in fasta format
#These input files were taken from uniprot database

import sys

#reads in the first sequence
with open(sys.argv[1], "r") as f:
    f.readline()
    seq1 = ""
    for line in f:
        seq1 += line.strip()

#reads in the second sequence
with open(sys.argv[2], "r") as f:
    f.readline()
    seq2 = ""
    for line in f:
        seq2 += line.strip()

#blosum stores the values in the BLOSUM62 matrix
#letters stores the single letter amino acid codes in the order they appear
blosum = []
letters = []
with open("BLOSUM62.txt", "r") as f:
    letters = f.readline().strip().split()
    for line in f:
        blosum.append(map(int,line.strip().split()))

#stores the indel penalty
indel = -5

#creates the alignment matrix and backtrack matrix
alignment = [[0] * (len(seq2) + 1) for i in range(len(seq1)+ 1)]
backtrack = [[0] * (len(seq2)+1) for i in range(len(seq1)+1)]

#intializes the first row and column to the indel penalties
#this would be equivalent to aligning either string to an empty string
for i in range(len(seq1)):
    alignment[i][0] = i * indel;
    backtrack[i][0] = "dwn"

for j in range(len(seq2)):
    alignment[0][j] = j * indel
    backtrack[0][j] = "rgt"


for i in range(1, len(seq1)+1):
    for j in range(1, len(seq2)+1):
        #finds wherew in the BLOSUM62 matrix the current lettters are
        ind1 = letters.index(seq1[i-1])
        ind2 = letters.index(seq2[j-1])

        #applies the recurrence relationship
        alignment[i][j] = max(alignment[i-1][j]+indel, alignment[i][j-1]+indel,
                              alignment[i-1][j-1]+blosum[ind1][ind2])

        #fills in the backtracking matrix
        if alignment[i][j] == alignment[i-1][j]+indel:
            backtrack[i][j] = "dwn"
        elif alignment[i][j] == alignment[i][j-1]+indel:
            backtrack[i][j] = "rgt"
        else:
            backtrack[i][j] = "diag"

#prints the score
print alignment[len(seq1)][len(seq2)]

#used to store the alignments themselves
align1 = []
align2 = []

#indices used to go backwards in bactrack matrix
i = len(seq1)
j - len(seq2)

while i > 0 or j > 0:
    #decrease both if prev is diagonal
    if backtrack[i][j] == "diag":
        align1.append(seq1[i-1])
        align2.append(seq2[j-1])
        i -= 1
        j -= 1
    #decrease only column if prev is right
    elif backtrack[i][j] == "rgt":
        align1.append("-")
        align2.append(seq2[j-1])
        j -= 1
    #decrease only row if prev is down
    else:
        align1.append(seq1[i-1])
        align2.append("-")
        i -= 1

#need to reverse these before printing
print "".join(align1[::-1])
print "".join(align2[::-1])

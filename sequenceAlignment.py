"""

Global sequence alignment based on dynamic programming. 
    python sequenceAlignment.py -i1 (seq1 file) -i2 (seq2 file) -s (scoring file) -g (gap penalty) 

    @Chuankai Zhao, czhao37@illinois.edu
    
"""

import argparse
import numpy as np

#Read the sequence fasta file. 
def readseq(seqname):
    seqfile = open(seqname,"r")
    lines   = [ line.rstrip() for line in seqfile ]
    lines   = lines[1:]
    seq     = ""
    for line in lines:
        seq = seq + line
    return seq

#Read the substitution matrix.
def readscores(score, gap):
    scorefile = open(score,"r")
    lines     = [ line.rstrip() for line in scorefile]
    lines     = lines[1:]
    scores    = []
    gap       = float(gap)
    for line in lines:
        line  = line.split("\t")
        line.append(gap)
        scores.append(line[1:])
    scores.append([gap, gap, gap, gap, 0])
    for i in range(5):
        for j in range(5):
            scores[i][j] = float(scores[i][j])
    scores     = np.array(scores)
    return scores

#Align the sequences.
def seqalign(seq1, seq2, scores, dict):
    M = len(seq1)
    N = len(seq2)
    max_scores = np.zeros((M+1,N+1))
    last_poses = []
    for row in range(M+1):
        last_poses_row = []
        for column in range(N+1):
            if row == 0:
                if column == 0:
                    max_scores[row][column] = 0.
                    last_pos = [row,column]
                    last_poses_row.append(last_pos)
                if column > 0: 
                    max_scores[row][column] = max_scores[row][column-1] + scores[dict["-"]][dict[seq2[column-1]]]
                    last_pos = [row,column - 1]
                    last_poses_row.append(last_pos)
            if row > 0:
                if column == 0: 
                    max_scores[row][column] = max_scores[row-1][column] + scores[dict[seq1[row-1]]][dict["-"]]
                    last_pos = [row, column]
                    last_poses_row.append(last_pos)
                if column > 0:
                    last_pos  = [row-1, column-1]
                    max       = max_scores[row-1][column-1] + scores[dict[seq1[row-1]]][dict[seq2[column-1]]]
                    max_right = max_scores[row-1][column]   + scores[dict[seq1[row-1]]][dict["-"]]
                    if max_right > max:
                        max = max_right
                        last_pos = [row-1, column]
                    max_left  = max_scores[row][column-1] + scores[dict["-"]][dict[seq2[column-2]]]
                    if max_left  > max:
                        max = max_left
                        last_pos = [row, column-1]
                    max_scores[row][column] = max
                    last_poses_row.append(last_pos)
        last_poses.append(last_poses_row)

    trace_path = [[M,N]]
    last_pos = last_poses[M][N]

    while True:
        if last_pos[0] == 0 and last_pos[1] == 0:  
            break
        if last_pos[0] > 0 or last_pos[1] > 0:
            trace_path.append(last_pos)
            last_pos = last_poses[last_pos[0]][last_pos[1]]
   
    trace_path_new = []
    for i in reversed(trace_path):
        trace_path_new.append(i)
    trace_path = trace_path_new
    Num = len(trace_path)
    
    coords = [0,0]
    Aligned_Seq1 = ""
    Aligned_Seq2 = ""
    for i in range(Num):
        element = trace_path[i]
        if element[0] > coords[0]:
            Aligned_Seq1 = Aligned_Seq1 + seq1[element[0]-1] 
        else:
            Aligned_Seq1 = Aligned_Seq1 + "-"
        if element[1] > coords[1]:
            Aligned_Seq2 = Aligned_Seq2 + seq2[element[1]-1]
        else:
            Aligned_Seq2 = Aligned_Seq2 + "-"
        coords  = element
  
    print(max_scores) 
    max_score = max_scores[M][N]

    return Aligned_Seq1, Aligned_Seq2, max_score

#Define the interactive command line. 
def parse_cmdln():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i1','--input1',dest='seq1',
                         help='File of sequence 1.', required=True)
    parser.add_argument('-i2','--input2',dest='seq2',
                         help='File of sequence 2.', required=True)
    parser.add_argument('-s', '--score', dest='score',
                         help='File of scoring matrix.', required=True)
    parser.add_argument('-g', '--gap', dest='gap',
                         help='Gap penalty for sequence alignment.', required=True)
    args   = parser.parse_args()
    return args

#Main program. 
if __name__ == '__main__':
    options = parse_cmdln()
    seq1    = readseq(options.seq1)
    seq2    = readseq(options.seq2)
    dict    = { "A":0 , "C":1 , "G":2 , "T":3  , "-":4 }
    scores  = readscores(options.score, options.gap)
    aligned_seq1, aligned_seq2, max_score = seqalign(seq1, seq2, scores, dict)
    file    = open("aligned.txt","w")
    file.write("The optimal alignment between given sequences has score %r.\n" % max_score)
    file.write(aligned_seq1)
    file.write("\n")
    file.write(aligned_seq2)
    file.write("\n")
    file.close()

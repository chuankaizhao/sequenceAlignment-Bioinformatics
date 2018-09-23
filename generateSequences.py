"""

Generating random sequences for testing the global sequence alignment algorithm. 
    python generateSequences.py 
    @Chuankai Zhao, czhao37@illinois.edu

"""

import random

#Create a random sequence.
def creat_seq(L):
    seq = ""
    for i in range(L):
        random.seed()
        num = random.random()
        if  0. <= num < 0.25: 
            seq = seq + "A"
        if  0.25 <= num < 0.5:
            seq = seq + "C"
        if  0.5  <= num < 0.75:
            seq = seq + "G"
        if  0.75 <= num <= 1.0: 
            seq = seq + "T"
    return seq

#Randomly mutate the sequences. 
def mutate_seq(seq,L):
    N = int(L/10)
    for i in range(N):
        random.seed()
        pos   = random.randrange(0,len(seq))
        var_1 = random.random()
        if  0. <= var_1 < 0.5:
            list  = ['A','C','G','T']
            list.remove(seq[pos])
            var_2 = random.random()
            if 0. <= var_2 < 1./3:
                seq = seq[0:pos] + list[0] + seq[pos+1:]
            if 1./3 <= var_2 < 2./3:
                seq = seq[0:pos] + list[1] + seq[pos+1:]
            if 2./3 <= var_2 < 1.:
                seq = seq[0:pos] + list[2] + seq[pos+1:]
        if 0.5 <= var_1 <= 1.:
            if pos == 0:
                seq = seq[1:]
            if pos != 0:
                seq = seq[0:pos] + seq[pos+1:]
    return seq

#Main program. 
if __name__ == "__main__":
    L    = 100
    seq1 = creat_seq(L)
    seq2 = mutate_seq(seq1,L)
    seq3 = mutate_seq(seq1,L)
    f1   = open("seq2.fasta", "w")
    f1.write(">seq2\n")
    f1.write(seq2)
    f1.close()
    f2   = open("seq3.fasta", "w")
    f2.write(">seq3\n")
    f2.write(seq3)
    f2.close()

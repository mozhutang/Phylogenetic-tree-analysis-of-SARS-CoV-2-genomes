import copy
import multiprocessing
import random
import sys

import numpy
from Bio import SeqIO


def main():

    all_species = []
    for seq_record in SeqIO.parse("cov.fasta", "fasta"):
        dna_seq=str(seq_record.seq)
        name=seq_record.description.split()[0]
        all_species.append([name,dna_seq])

    Tree(all_species)

class Tree:
    def __init__( self, dna_seq ):
        
        self.name = []
        self.seq = []
        for i in range(len(dna_seq)):
            self.name.append(dna_seq[i][0])
            self.seq.append(dna_seq[i][1])

        self.tree = generate_tree( self.seq )
        T = 10
        clades = {}
        clade_inti( self.tree, clades )

        for i in range ( 0, T) :
            sample_seq = sampling(self.seq)
            sample_tree = generate_tree(sample_seq)

            print("Bootstrap Tree ", i)
            print(sample_tree)

            clade_update(sample_tree, clades)

        confidences = {}
        for i, j in clades.items():
            confidences[ i ] = j / T
        print("confidences:", confidences)

        print("Phylogenetic tree:", self.tree)
        for i in range(len(self.name)):
            print(i, "=", self.name[i])
        
def generate_tree( seq ):
    d_matrix = numpy.zeros((len(seq), len(seq)), dtype=float)
    def store(result):
        d = result[0]
        i = result[1]
        j = result[2]
        d_matrix[i, j] = d

    pool = multiprocessing.Pool()

    for i in range( 0, len( seq ) ):
        for j in range( i+1, len( seq ) ):
            pool.apply_async(calculate_distance, args =  (seq[i], seq[j], i, j), callback = store)
    
    pool.close()
    pool.join()
    return calculate_tree( d_matrix )

def calculate_distance( seq1, seq2, i, j):
    return [dis( seq1, seq2 ), i, j]

def sampling( seq ):
    sample = ["" for i in seq]

    for i in range(len(seq[0])):
        index = random.randint(0, len(seq[0]) - 1)
        for j in range(len(seq)):
            sample[j] += seq[j][index]

    return sample

def dis( seq1, seq2): 
    d=0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            d+=1
    return d
#TODO: maybe more complicated method for distance


def clade_inti( tree, clade ):
    if type( tree[0] ) == tuple:
        clade_inti( tree[0], clade )
    if type( tree[1] ) == tuple:
        clade_inti( tree[1], clade )
    clade[tree] = 0

def clade_update( tree, clade ):
    if type( tree[0] ) == tuple:
        clade_update( tree[0], clade )
    if type( tree[1] ) == tuple:
        clade_update( tree[1], clade )

    if tree in clade:
        clade[ tree ] += 1

def calculate_tree( dis ):
    dis = dis.tolist()

    t = []
    map = []
    for i in range( 0, len( dis ) ):
        t.append(i)
        map.append( [i] )
    copy_dis = copy.deepcopy(dis)
    
    while len(copy_dis)>2:
        min = copy_dis[0][1]
        min_index = [0, 1]

        for i in range( 0, len( copy_dis ) ):
            for j in range( i+1, len( copy_dis[i] ) ):
                if copy_dis[i][j] < min:
                    min_index[0] = i
                    min_index[1] = j
                    min = copy_dis[i][j]
        X = min_index[0]
        Y = min_index[1]

        t[X] = ( t[X], t[Y] )
        t.pop(Y)

        map[X] += map[Y]
        map.pop(Y)
        copy_dis.pop( Y )
        for i in range( 0, len( copy_dis ) ):
            copy_dis[i].pop( Y )

        for i in range(0, X):
            copy_dis[i][X] = cal_avg( map[i], map[X], dis )
        for j in range( X+1, len( copy_dis ) ):
            copy_dis[X][j] = cal_avg( map[X], map[j], dis )
    return ( t[0], t[1] )

def cal_avg( seq1, seq2, dis ):
    sum = 0
    count = 0
    for i in seq1:
        for j in seq2:
            sum += dis[min( i, j )][max( i, j )]
            count += 1
    return sum / count


if __name__ == '__main__':
    main()

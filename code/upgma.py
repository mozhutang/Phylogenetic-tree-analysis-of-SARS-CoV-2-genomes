import copy
import numpy
import multiprocessing

def calculate_upgma( distance_matrix ):
    print(distance_matrix)

    distance_matrix = distance_matrix.tolist()

    # will hold the newick expression of the tree
    seq_tree = []

    # will map the indexes of the new matrix to the list 
    # of sequences that they represent
    seq_map = []

    # start them off with just 0 through len of matrix
    for i in range( 0, len( distance_matrix ) ):
        seq_tree.append(i)
        seq_map.append( [i] )

    new_matrix = copy.deepcopy( distance_matrix )
    
    while len( new_matrix ) > 2:

        min_index = find_lowest_distance( new_matrix )

        X = min_index[0]
        Y = min_index[1]

        seq_tree[X] = ( seq_tree[X], seq_tree[Y] )
        seq_tree.pop(Y)

        seq_map[X] += seq_map[Y]
        seq_map.pop(Y)

        # remove the row and column that will not be used
        # to store the average distance once both sequences
        # X and Y are combined
        drop_seq( new_matrix, Y )

        for i in range(0, X):
            new_matrix[i][X] = \
            find_avg_distance( seq_map[i], seq_map[X], distance_matrix )

        for j in range( X+1, len( new_matrix ) ):
            new_matrix[X][j] = \
            find_avg_distance( seq_map[X], seq_map[j], distance_matrix )

    print( seq_tree[0], seq_tree[1])

    return ( seq_tree[0], seq_tree[1] )

def find_lowest_distance( distance_matrix ):
    min_val = distance_matrix[0][1]
    min_index = [0, 1]

    for i in range( 0, len( distance_matrix ) ):
        for j in range( i+1, len( distance_matrix[i] ) ):
            if distance_matrix[i][j] < min_val:
                min_index[0] = i
                min_index[1] = j
                min_val = distance_matrix[i][j]

    return min_index

def drop_seq( matrix, coordinate ):
    # remove the row at the coordinate
    matrix.pop( coordinate )

    # remove the column at the coordinate
    for i in range( 0, len( matrix ) ):
        matrix[i].pop( coordinate )

def find_avg_distance( seq_list_1, seq_list_2, distance_matrix ):
    sum = 0
    count = 0
    for seq_1_index in seq_list_1:
        for seq_2_index in seq_list_2:
            # since our matrix is symetric and we only fill the 
            # values in the upper right half, this min max thing 
            # is to ensure we use the half of the matrix that 
            # isn't all of the zeros
            x = min( seq_1_index, seq_2_index )
            y = max( seq_1_index, seq_2_index )

            sum += distance_matrix[x][y]
            count += 1

    return sum / count

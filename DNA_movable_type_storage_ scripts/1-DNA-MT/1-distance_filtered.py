#This script filters out DNA-MTs with an editing distance of at least 6 from a pool of 20,000 potential sequences
import sys
import distance
import time
import random

def comp_all(sequence, mw_list,n=10):
    res = False
    for s in mw_list:
        sim = distance.levenshtein(s, sequence)  ###
        if (sim <n):
            res = True
            break
    return res

if __name__ == '__main__':
    fdis_seq=open("distance_only_6.res","w")
    fdis_seq6 = open("seq_distance_only_6.res", "w")
    fdis_time=open("distance_only-6.time","w")
    ok_count = 0
    seq_path = sys.argv[1]
    fsequences = open(seq_path, "r")
    sequences = {}
    sequences_list=[]
    for seq in fsequences.readlines():
        sequences[seq.strip('\n').split(' ')[0]] = seq.strip('\n').split(' ')[1]
        sequences_list.append(seq.strip('\n').split(' ')[0])

    print("read file",len(sequences))

    sequences_length=len(sequences)

    seed_length=10
    seed_sequence={}
    while(len(seed_sequence)<=seed_length):
        ran_int=random.randint(0,sequences_length-1)
        seed_sequence[sequences_list[ran_int]]=0


    for seed in seed_sequence:
        distance_seq_list={}
        max_seq=0
        start=time.perf_counter()
        for s in sequences:
            if (len(distance_seq_list) == 0):
                distance_seq_list[seed] = sequences[seed]
            elif (comp_all(s, distance_seq_list.keys(), n=6) == False):
                distance_seq_list[s] = sequences[s]
        end = time.perf_counter()
        print(seed,end - start, file=fdis_time)
        seed_sequence[seed]=len(distance_seq_list)
        if(max_seq<seed_sequence[seed]):
            max_seq=seed_sequence[seed]

    for key in seed_sequence:
        print(key,seed_sequence[key],file=fdis_seq)
    print(max_seq,file=fdis_seq)

    for key in distance_seq_list.keys():
        print(key,distance_seq_list[key],file=fdis_seq6)










        
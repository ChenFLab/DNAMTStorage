#This script filtered 20,000 potential DNA-MTs
import random
import math
import time
import os
import distance
fseed20=open("mw_produce200000.res","w")
ftime=open("mw_produce200000.time","w")


def DNA_complement(sequence):  ###
    complement = {'C': 'G', 'G': 'C', 'T': 'A', 'A': 'T'}
    res = ''
    for chr in sequence:
        res = res + complement[chr]
    return res


def DNA_reverse(sequence):  ##
    return sequence[::-1]


#
def repeatable_permutation(n=10, s=[0, 0, 0, 0, 0, 0, 0, 0, 0,0], cur=0, k=4,
                           seeds=[]):  #
    dna_map = ['A', 'T', 'G', 'C']
    if (cur == n):
        seq = ''.join(dna_map[elem] for elem in s)
        if (seq.find("AAAAA") == -1 and seq.find("TTTTT") == -1 and seq.find("CCCC") == -1 and seq.find("GGGG") == -1):
            if (seq_quality(seq,down=0.2,up=0.8)):
                if (has_rev_compl(seq, 4) == False):
                    if (has_rev_2chains_connect(seq, 3) == False):
                        if (has_rev_2chains_gap(seq, 5) == False):
                            bases = 0
                            if (seq.find('A') != -1):
                                bases = bases + 1
                            if (seq.find('T') != -1):
                                bases = bases + 1
                            if (seq.find('C') != -1):
                                bases = bases + 1
                            if (seq.find('G') != -1):
                                bases = bases + 1
                            if (bases >= 3):
                                seeds.append(seq)
    else:
        for i in range(0, k):
            s[cur] = i
            repeatable_permutation(n=n, s=s, cur=cur + 1, k=4)
    return seeds

def seq_quality(s,down=0.4,up=0.6):  ###
    res=False
    ca = s.count('A')
    ct = s.count('T')
    cg = s.count('G')
    cc = s.count('C')
    score = round(float(float((cc + cg))/ float(ca + ct + cc + cg)),2)
    if(score-down>0.000001 and score-up<0.000001):
        res=True
    return res
#######
def has_subseq(seq, conser, RECsite):
    res = False
    count = 0  #
    for ss in conser:
        ss = ss.upper()
        ss_count=seq.count(ss)
        count = ss_count + count
    if (count > 0 and RECsite == False):
        res = True
    if (count > 1 and RECsite == True):
        res = True
    return res

def has_rev_compl(sequence, n=4):
    res = False
    l = len(sequence)
    i = 0
    if (l >= 2 * n):
        while (i < l - n):
            s1 = sequence[i:i + n]
            s2 = sequence[i + n:]
            rs2 = DNA_reverse(DNA_complement(s2))
            if (s1 in rs2):
                res = True
            i = i + 1
    return res

#########################################################5'-3' 3'-5'
def has_rev_2chains_gap(s, vn=8):
    res = False
    s_rev = s[::-1]
    i = 0
    l = len(s)
    while (i < l and res == False):
        count=0
        sr1 = s_rev[0:l - i]
        k = 0
        while (k < len(sr1)):
            if (sr1[k] == DNA_complement(s[i + k])):
                count = count + 1
            k = k + 1
        if (count >= vn):
            res = True
            break
        i = i + 1

    j = 0
    while (j < l and res == False):
        count = 0
        sr2 = s_rev[j:l]
        k = 0
        while (k < len(sr2)):
            if (sr2[k] == DNA_complement(s[k])):
                count = count + 1
            k = k + 1
        if (count >= vn):
            res = True
            break
        j = j + 1
    return res

def has_rev_2chains_connect(s, n=3):
    res = False
    s_rev = DNA_complement(s)[::-1]
    i = 0
    l = len(s)
    while (i <= l - n):
        if (s[i:i + n] in s_rev):
            #print(s[i:i + n],s_rev)
            res = True
            break
        i = i + 1
    return res

def comp_all(sequence, mw_list,n=10):
    res = False
    for s in mw_list:
        sim = distance.levenshtein(s, sequence)  ###
        if (sim <n):
            res = True
            break
    return res

########################################################
def seq_random(n, conser,left, right, number,seeds,movable_dict,RECsite=False):  #
    length = len(seeds)
    sptimes = math.ceil(n / len(seeds[0]))
    number_seq = []
    while (True):
        if (len(number_seq) >= number):
            return number_seq
            break
        s = ''
        for i in range(0, sptimes):
            rt = random.randint(0, length - 1)
            s = s + seeds[rt]
        s=s[0:n]
        left = left.upper()
        right = right.upper()
        ss = left + s + right


        if (ss.find('AAAAA') == -1 and ss.find('TTTTT') == -1 and ss.find('CCCC') == -1 and ss.find('GGGG') == -1):
            #print("1 polymer yes")  #
            if (has_subseq(seq=ss, conser=conser, RECsite=RECsite) == False ):
                #print("2 subseq:yes")  #
                if (seq_quality(s,down=0.4,up=0.6)):  #
                    #print("3 gc content yes")
                    if (has_rev_compl(ss, 4) == False):
                        #print("4 hairpin:yes")  ##hairpin4
                        if (RECsite == True and ss.find(conser[0]) != -1):
                            self_dimer_s = s + right
                        elif (RECsite == True and ss.find(conser[1]) != -1):
                            self_dimer_s = left + s
                        else:
                            self_dimer_s = ss
                        if (has_rev_2chains_connect(self_dimer_s, 3) == False):
                            #print("5 self-dimer-connect yes")
                            if (has_rev_2chains_gap(self_dimer_s, 8) == False):
                                #print("6 self-dimer-gap yes")
                                if("AGTACT" not in ss):##Scal
                                    if(s not in movable_dict.keys()):
                                        number_seq.append(s)

    return number_seq


if __name__ == '__main__':
    if (os.path.exists("seeds10.res") == False):
        start = time.perf_counter()
        fseed10 = open("seeds10.res", "w")
        seeds = repeatable_permutation()
        for seed in seeds:
            print(seed, file=fseed10)
        fseed10.close()
        end = time.perf_counter()
        print("seed10", end - start)
    start = time.perf_counter()
    seeds = []
    fseed = open("seeds10.res", "r")
    for seed in fseed.readlines():
        seeds.append(seed.strip('\n'))
    end = time.perf_counter()
    print("read seeds10 file", end - start)
    movable_dict={}
    conser = (
        "AAGCTT",
        "GAATTC"
    )
    re_conser = (
        "AAGCTT",
        "GAATTC"
    )
    mw1_count=0
    mw2_count=0
    mw2_count = 0
    mw3_count = 0
    mw4_count = 0
    mw5_count = 0
    mw6_count = 0
    mw7_count = 0
    mw8_count = 0
    mw9_count = 0
    mw10_count = 0
    while(len(movable_dict)<=200000):

        start=time.perf_counter()
        mw1_list = seq_random(20, conser=conser,  left='aagctt', right='ctcc', number=64,movable_dict=movable_dict,seeds=seeds, RECsite=True)
        for mw1 in mw1_list:
            if(mw1 not in movable_dict.keys()):
                mw1_count=mw1_count+1
                movable_dict[mw1]="mw1-"+str(mw1_count)
        end = time.perf_counter()
        print("mw1-64", end - start,file=ftime)

        start = time.perf_counter()
        mw2_list = seq_random(20, conser=conser,  left='ctcc', right='acta', number=64,movable_dict=movable_dict, seeds=seeds)
        for mw2 in mw2_list:
            if(mw2 not in movable_dict.keys()):
                mw2_count=mw2_count+1
                movable_dict[mw2]="mw2-"+str(mw2_count)
        end = time.perf_counter()
        print("mw2-64", end - start,file=ftime)

        start = time.perf_counter()
        mw3_list = seq_random(20, conser=conser,  left='acta', right='ggta', number=64, movable_dict=movable_dict,seeds=seeds)
        for mw3 in mw3_list:
            if(mw3 not in movable_dict.keys()):
                mw3_count=mw3_count+1
                movable_dict[mw3]="mw3-"+str(mw3_count)
        end = time.perf_counter()
        print("mw3-64", end - start,file=ftime)

        start = time.perf_counter()
        mw4_list = seq_random(20, conser=conser,  left="ggta", right="atag", number=64,movable_dict=movable_dict,seeds=seeds)
        for mw4 in mw4_list:
            if(mw4 not in movable_dict.keys()):
                mw4_count=mw4_count+1
                movable_dict[mw4]="mw4-"+str(mw4_count)
        end = time.perf_counter()
        print("mw4-64", end - start,file=ftime)

        start = time.perf_counter()
        mw5_list = seq_random(20, conser=conser,  left="atag", right="aagt", number=64,movable_dict=movable_dict, seeds=seeds)
        for mw5 in mw5_list:
            if(mw5 not in movable_dict.keys()):
                mw5_count=mw5_count+1
                movable_dict[mw5]="mw5-"+str(mw5_count)
        end = time.perf_counter()
        print("mw5-64", end - start,file=ftime)

        start = time.perf_counter()
        mw6_list = seq_random(20, conser=conser,  left="aagt", right="ctga", number=64,movable_dict=movable_dict,seeds=seeds)
        for mw6 in mw6_list:
            if(mw6 not in movable_dict.keys()):
                mw6_count=mw6_count+1
                movable_dict[mw6]="mw6-"+str(mw6_count)
        end = time.perf_counter()
        print("mw6-64", end - start,file=ftime)

        start = time.perf_counter()
        mw7_list = seq_random(20, conser=conser, left='ctga', right='caac', number=64,movable_dict=movable_dict,seeds=seeds)
        for mw7 in mw7_list:
            if(mw7 not in movable_dict.keys()):
                mw7_count=mw7_count+1
                movable_dict[mw7]="mw7-"+str(mw7_count)
        end = time.perf_counter()
        print("mw7-64", end - start,file=ftime)

        start = time.perf_counter()
        mw8_list = seq_random(20, conser=conser, left='caac', right='agga', number=64, movable_dict=movable_dict,
                              seeds=seeds)
        for mw8 in mw8_list:
            if (mw8 not in movable_dict.keys()):
                mw8_count = mw8_count + 1
                movable_dict[mw8] = "mw8-" + str(mw8_count)
        end = time.perf_counter()
        print("mw8-64", end - start, file=ftime)

        start = time.perf_counter()
        mw9_list = seq_random(20, conser=conser, left='agga', right='tgtg', number=64, movable_dict=movable_dict,
                              seeds=seeds)
        for mw9 in mw9_list:
            if (mw9 not in movable_dict.keys()):
                mw9_count = mw9_count + 1
                movable_dict[mw9] = "mw9-" + str(mw9_count)
        end = time.perf_counter()
        print("mw9-64", end - start, file=ftime)

        start = time.perf_counter()
        mw10_list = seq_random(20, conser=conser, left='tgtg', right='gaattc', number=64, movable_dict=movable_dict,
                              seeds=seeds, RECsite=True)
        for mw10 in mw10_list:
            if (mw10 not in movable_dict.keys()):
                mw10_count = mw10_count + 1
                movable_dict[mw10] = "mw10-" + str(mw10_count)
        end = time.perf_counter()
        print("mw10-64", end - start, file=ftime)



    for mw in movable_dict.keys():
        print(mw ,movable_dict[mw],file=fseed20)
    fseed20.close()












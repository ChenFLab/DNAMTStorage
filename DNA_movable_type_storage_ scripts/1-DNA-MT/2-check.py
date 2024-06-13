#This script check all DNA-MTs
import time
import sys
import distance


def DNA_complement(sequence):  ###
    complement = {'C': 'G', 'G': 'C', 'T': 'A', 'A': 'T'}
    res = ''
    for chr in sequence:
        res = res + complement[chr]
    return res


def DNA_reverse(sequence):  ##
    return sequence[::-1]


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
def check_seq(s, conser,left, right, RECsite=False):  #
    res=False
    left = left.upper()
    right = right.upper()
    ss = left + s + right
    if (ss.find('AAAAA') == -1 and ss.find('TTTTT') == -1 and ss.find('CCCC') == -1 and ss.find('GGGG') == -1):
        if (has_subseq(seq=ss, conser=conser, RECsite=RECsite) == False ):
            if (seq_quality(s,down=0.4,up=0.6)):  #
                if (has_rev_compl(ss, 4) == False):
                    if (RECsite == True and ss.find(conser[0]) != -1):
                        self_dimer_s = s + right
                    elif (RECsite == True and ss.find(conser[1]) != -1):
                        self_dimer_s = left + s
                    else:
                        self_dimer_s = ss
                    if (has_rev_2chains_connect(self_dimer_s, 3) == False):
                        if (has_rev_2chains_gap(self_dimer_s, 8) == False):
                            if("AGTACT" not in ss):##Scal
                                    res=True
    return res


if __name__ == '__main__':
    conser = (
        "AAGCTT",
        "GAATTC"
    )
    re_conser = (
        "AAGCTT",
        "GAATTC"
    )
    sfp=sys.argv[1]
    fsfp=open(sfp,"r")
    seq_dic={ }
    start = time.perf_counter()

    for line in fsfp.readlines():
        line=line.strip('\n')
        seq_dic[line.split(' ')[0]]=line.split(' ')[1].split('-')[0]
    end = time.perf_counter()
    print("read seq file", end - start)

    for key in seq_dic.keys():
        res=True
        if(seq_dic[key]=="mw1"):
            res = check_seq(key, conser=conser, left='aagctt', right='ctcc', RECsite=True)
        if (seq_dic[key] == "mw2"):
            res = check_seq(key, conser=conser, left='ctcc', right='acta', RECsite=False)
        if (seq_dic[key] == "mw3"):
            res = check_seq(key, conser=conser, left='acta', right='ggta', RECsite=False)
        if (seq_dic[key] == "mw4"):
            res = check_seq(key, conser=conser, left="ggta", right="atag", RECsite=False)
        if (seq_dic[key] == "mw5"):
            res = check_seq(key, conser=conser, left="atag", right="aagt", RECsite=False)
        if (seq_dic[key] == "mw6"):
            res = check_seq(key, conser=conser, left="aagt", right="ctga", RECsite=False)
        if (seq_dic[key] == "mw7"):
            res = check_seq(key, conser=conser, left='ctga', right='caac', RECsite=False)
        if (seq_dic[key] == "mw8"):
            res = check_seq(key, conser=conser, left='caac', right='agga', RECsite=False)
        if (seq_dic[key] == "mw9"):
            res = check_seq(key, conser=conser, left='agga', right='tgtg', RECsite=False)
        if (seq_dic[key] == "mw10"):
            res = check_seq(key, conser=conser, left='tgtg', right='gaattc', RECsite=True)
        if(res==False):
            print("sequence not pass")
    print("single sequence done")
    for key1 in seq_dic.keys():
        for key2 in seq_dic.keys():
            if(key1!=key2):
                dis=distance.levenshtein(key1,key2)
                if(dis<6):
                    print(key1,key2,dis,seq_dic[key1],seq_dic[key2])
    print("double sequence done")
    fsfp.close()












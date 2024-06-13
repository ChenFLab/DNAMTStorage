#This script filter 200 4-nt linkers

def seq_quality(s,down=0.4,up=0.6):  ###GC content
    res=False
    ca = s.count('A')
    ct = s.count('T')
    cg = s.count('G')
    cc = s.count('C')
    score = round(float(float((cc + cg))/ float(ca + ct + cc + cg)),2)
    if(score-down>0.000001 and score-up<0.000001):
        res=True
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
                # print(s1,s2)
                res = True
            i = i + 1
    return res

def DNA_complement(sequence):  
    complement = {'C': 'G', 'G': 'C', 'T': 'A', 'A': 'T'}
    res = ''
    for chr in sequence:
        res = res + complement[chr]
    return res

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

def DNA_reverse(sequence):  
    return sequence[::-1]

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
                # print(s1,s2)
                res = True
            i = i + 1
    return res

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

def repeatable_permutation(n, s,down,up, cur=0, k=4,seeds=[]):  
    dna_map = ['A', 'T', 'G', 'C']
    if (cur == n):
        seq = ''.join(dna_map[elem] for elem in s)
        if (seq.find("AAA") == -1 and seq.find("TTT") == -1 and seq.find("CCC") == -1 and seq.find("GGG") == -1):
            if (seq_quality(seq,down=down,up=up)):
                if (has_rev_compl(seq, 4) == False):
                    if (has_rev_2chains_connect(seq, 3)== False):
                        if (has_rev_2chains_gap(seq, 3) == False):
                            bases = 0
                            if (seq.find('A') != -1):
                                bases = bases + 1
                            if (seq.find('T') != -1):
                                bases = bases + 1
                            if (seq.find('C') != -1):
                                bases = bases + 1
                            if (seq.find('G') != -1):
                                bases = bases + 1
                            if (bases >= 2):
                                seeds.append(seq)
    else:
        for i in range(0, k):
            s[cur] = i
            repeatable_permutation(n=n, s=s, down=down,up=up,cur=cur + 1, k=4)

    return seeds
def check(seq):
    if (seq.find("AAA") == -1 and seq.find("TTT") == -1 and seq.find("CCC") == -1 and seq.find("GGG") == -1):
        if (seq_quality(seq, down=0.2, up=0.8)):
            if (has_rev_compl(seq, 4) == False):
                if (has_rev_2chains_connect(seq, 3) == False):
                    if (has_rev_2chains_gap(seq, 3) == False):
                        bases = 0
                        if (seq.find('A') != -1):
                            bases = bases + 1
                        if (seq.find('T') != -1):
                            bases = bases + 1
                        if (seq.find('C') != -1):
                            bases = bases + 1
                        if (seq.find('G') != -1):
                            bases = bases + 1
                        if (bases >= 2):
                            print("ok")


if __name__=='__main__':
    conser = (
        "AAGCTT",
        "GAATTC"
    )
    re_conser = (
        "AAGCTT",
        "GAATTC"
    )
    f=open("linkers.txt","w")

    
    seeds=repeatable_permutation(n=4,s=[0,0,0,0],down=0.2,up=0.8)
    for se in seeds:
        print(se,file=f)











        
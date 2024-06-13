#This script scoring 5000 sets of linkers.
import random
def DNA_complement(sequence):  
    complement = {'C': 'G', 'G': 'C', 'T': 'A', 'A': 'T'}
    res = ''
    for chr in sequence:
        res = res + complement[chr]
    return res


def max_bp(s1, s2,ps1,ps2,n=10):
    s=""
    add_score=0
    m = [[0 for i in range(len(s2) + 1)] for j in range(len(s1) + 1)]
    #print(m)
    mmax = 0  # The length of the longest match
    p = 0  # The last position of the longest match in s1
    for i in range(len(s1)):
        for j in range(len(s2)):
            if s1[i] == s2[j]:
                m[i + 1][j + 1] = m[i][j] + 1
                if m[i + 1][j + 1] > mmax:
                    mmax = m[i + 1][j + 1]
                    p = i + 1
    s = s1[p - mmax:p]
    ls=len(s)
    l=len(s1)
    is1=s1.index(s)
    is2=s2.index(s)
    #print(s1,s2,mmax,n)

    if(n>mmax):
        add_score=add_score+2
        #print(s1, s2, is1, is2, s)
    elif(n==mmax and abs(ps2-ps1)>=2):
        add_score=add_score+1
        #print(s1,s2,is1,is2,s)
    else:
        return 0

    #print(s1,s2,mmax,n,cut_score)
    #print("maxscore",add_score)
    return add_score  

def max_bp_gap(s,smw,ps1,ps2,n=10):
    add_score=0
    s_rev = smw
    i = 0
    l = len(s)
    while (i < l ):
        count = 0
        sr1 = s_rev[0:l - i]
        k = 0
        tem1=""
        while (k < len(sr1)):
            if (sr1[k] == s[i + k]):
                tem1=tem1+sr1[k]
                count = count + 1
            else:
                tem1=tem1+" "
            k = k + 1
        tem1a = tem1.split(' ')
        str_len = list(map(len, tem1a))
        if (count != max(str_len)):  # 
            if (count < n):  # 
                add_score = add_score + 2
            elif (count == n and abs(ps2 - ps1) >= 2):  # 
                add_score = add_score + 1
            else:  # 
                return 0
        else:
            add_score = add_score + 2
        i = i + 1

    j = 0
    while (j < l ):
        tem2 = ""
        count = 0
        sr2 = s_rev[j:l]
        #print(sr2)
        k = 0
        while (k < len(sr2)):
            if (sr2[k] == s[k]):
                count = count + 1
                tem2 = tem2+sr2[k]
            else:
                tem2=tem2+" "
            k = k + 1
        #print(count)
        tem2a=tem2.split(' ')
        str_len = list(map(len, tem2a))
        #print("count2", count, abs(ps2 - ps1))
        j = j + 1
        if (count != max(str_len)):  
            if (count < n):  
                add_score = add_score + 2
                #print("count<n")
            elif (count == n and abs(ps2 - ps1) >= 2):  
                add_score = add_score + 1
                #print("count == n")
            else: 
                return 0
        else:
            add_score = add_score + 2
        #print("count != max(str_len)")
    #print("gap", add_score)


    return add_score

def check(s ,st,ps1,ps2):
    discard=False
    pair_score1=0
    pair_score2=0
    pair_score3=0
    pair_score4=0

    if (s != st and discard==False):
        bpsc1 = max_bp(s, st, ps1=ps1, ps2=ps2, n=3)
        if(bpsc1!=0):
            gbpsc1 = max_bp_gap(s, st, ps1=ps1, ps2=ps2, n=3)
            if(gbpsc1!=0):
                if(bpsc1==1 or gbpsc1==1):
                    pair_score1=0
                else:
                    pair_score1=1
            else:
                discard = True
        else:
            discard = True
    else:
        discard=True

    if (s != st[::-1]and discard==False):
        pair_score2 = 1
    else:
        discard = True


    if (s != DNA_complement(st)and discard==False):
        pair_score3 = 1
    else:
        discard = True


    if (s != DNA_complement(st[::-1]) and discard == False):
        bpsc4 = max_bp(s, DNA_complement(st[::-1]), ps1=ps1, ps2=ps2, n=3)
        if (bpsc4 != 0):
            gbpsc4 = max_bp_gap(s, DNA_complement(st[::-1]), ps1=ps1, ps2=ps2, n=3)
            if (gbpsc4 != 0):
                if (bpsc4 == 1 or gbpsc4 == 1):
                    pair_score4 = 0
                else:
                    pair_score4 = 1
            else:
                discard = True
        else:
            discard = True
    else:
        discard = True

    #print(pair_score1,pair_score2,pair_score3,pair_score4)
    return discard,pair_score1,pair_score2,pair_score3,pair_score4

def check_10(apl):
    apl.insert(0,"agct")
    apl.append("aatt")
    l=len(apl)
    #print("apl",l)
    score=0
    score1=0
    score2=0
    score3=0
    score4=0
    good_sticky=0
    for j in range(0, l):
        for k in range(j+1,l):
        
            disc,sdt1,sdt2,sdt3,sdt4=check(apl[j].upper(), apl[k].upper(),j,k)
            if(disc==True):
                score=0
                break
            else:
                score1 = score1 + sdt1
                score2 = score2 + sdt2
                score3 = score3 + sdt3
                score4 = score4 + sdt4
                score = score1+score2+score3+score4

        if(disc==True):
            score = 0
            break

    return disc,score,score1,score2,score3,score4




import os
if __name__=='__main__':
    print(os.getcwd())
    optim=["ctcc","atgc","ggta","atag","agac","ctga","caac","acga","tgtg"]
    fact=["CTCC","ACTA","GGTA","ATAG","AAGT","CTGA","CAAC","AGGA","TGTG"]

    print("optim",check_10(optim))
    print("fact", check_10(fact))

    coo=0
    for j in range(0, 11):
        for k in range(j + 1, 11):
            if((k-j)>=2):
                coo = coo + 1
            #print("saj",j,k)

    #print(coo)

    f = open("stick_list_4.txt", "r")
    f1 = open("sticky_composite4addnew2.txt", "w")
    sl = []
    for i in f.readlines():
        key = i.strip('\n')
        sl.append(key)

    for f in optim:
        #print(f)
        if(f.upper() not in sl):
            print("no",f)


    random_list = []
    count=0
    while (len(random_list) < 500):
        if ((count % 10000) == 0):
            print(count)
        L1 = tuple(random.sample(range(0, len(sl)), 9))
        # print(L1)
        if L1 not in random_list:
            count=count+1
            random_list.append(L1)
            apl = []
            for i in range(0, 9):
                apl.append(sl[L1[i]])
            # print(apl)

            disc,score,score1,score2,score3,score4 = check_10(apl)
            #print(score)
            #print(score,score1,score2,score3,score4,apl)
            if(disc==False):
                print(score,file=f1)
                #print(score)
    print(count)








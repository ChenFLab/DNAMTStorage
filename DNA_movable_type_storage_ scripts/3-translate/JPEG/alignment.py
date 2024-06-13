import sys
from Bio import pairwise2

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas import Series,DataFrame
import seaborn as sns
import palettable
from sklearn import datasets

seq1 = "ABCDEFGHIJ"

symbol = {"index1": "A", "index2": "B", "index3": "C", "index4": "D", "elem1": "E", "elem2": "F"
    , "elem3": "G", "elem4": "H", "elem5": "I", "check": "J"}

decode_dic_path = sys.argv[1]
nameID = decode_dic_path.split('/')[-1].split('.')[0]
fdecode_dic = open(decode_dic_path, "r")

pos_dup_path = nameID + "-pos-dup.res"
npos_ndup_path = nameID + "-npos-ndup.res"
pos_ndup_path = nameID + "-pos-ndup.res"
npos_dup_path = nameID + "-npos-dup.res"
fpos_dup = open(pos_dup_path, "w")
fnpos_ndup = open(npos_ndup_path, "w")
fpos_ndup = open(pos_ndup_path, "w")
fnpos_dup = open(npos_dup_path, "w")

syboml_dic_dict = {}
total_seq = 0
for line in fdecode_dic.readlines():
    if ("()" in line):
        print(line.strip('\n'), file=fpos_ndup)
        continue
    if (",)" in line):
        line = line.replace(',)', ")")
    line = line.strip()
    line = line.replace("\'", "")
    line = line.replace(")", "")
    line = line.replace("(", "")
    line = line.replace(" ", "")
    seq_list = line.split(',')
    i = 0
    l = len(seq_list)
    seq = ""
    while (i < l - 1):
        lll = (seq_list[i].split('-')[0])
        if(lll!=''):
            seq = seq + symbol[lll]
        i = i + 1
    if (seq in syboml_dic_dict.keys()):
        syboml_dic_dict[seq] = int(seq_list[-1]) + syboml_dic_dict[seq]
    else:
        syboml_dic_dict[seq] = int(seq_list[-1])
    total_seq = total_seq + int(seq_list[-1])

count_pos_dup = 0
count_npos_ndup = 0
count_npos_dup = 0
count_pos_ndup = 0
nw_count_pos_ndup = 0
pos_ndup_dict = {}
splice_count_dict = {}
for seq in syboml_dic_dict.keys():
    index1 = seq.find('A')
    if (index1 == -1):
        index1 = 0
    index2 = seq.find('B')
    if (index2 == -1):
        index2 = index1
    index3 = seq.find('C')
    if (index3 == -1):
        index3 = index2
    index4 = seq.find('D')
    if (index4 == -1):
        index4 = index3
    index5 = seq.find('E')
    if (index5 == -1):
        index5 = index4
    index6 = seq.find('F')
    if (index6 == -1):
        index6 = index5
    index7 = seq.find('G')
    if (index7 == -1):
        index7 = index6
    index8 = seq.find('H')
    if (index8 == -1):
        index8 = index7
    index9 = seq.find('I')
    if (index9 == -1):
        index9 = index8
    index10 = seq.find("J")
    if(index10==-1):
        index10=index9

    pos = (index1 <= index2 <= index3 <= index4 <= index5 <= index6 <= index7<=index8<=index9<=index10)  # 位置正确
    dup1 = seq.count('A')
    dup2 = seq.count('B')
    dup3 = seq.count('C')
    dup4 = seq.count('D')
    dup5 = seq.count('E')
    dup6 = seq.count('F')
    dup7 = seq.count('G')
    dup8= seq.count('H')
    dup9= seq.count('I')
    dup10=seq.count('J')
    ndup = (dup1 <= 1 and dup2 <= 1 and dup3 <= 1 and dup4 <= 1 and dup5 <= 1 and dup6 <= 1 and dup7 <= 1 and dup8<=1 and dup9<=1 and dup10<=1)  # 无重复

    ##位置错乱，无重复
    if (not pos and ndup):
        count_npos_ndup = count_npos_ndup + syboml_dic_dict[seq]
        print(seq, syboml_dic_dict[seq], file=fnpos_ndup)

    ##位置正确时，无重复

    if (pos and ndup):
        alignments = pairwise2.align.globalxx(seq, seq1)
        for alignment in alignments:
            #splice_count = len(alignment.seqA) - alignment.seqA.count('-')  ##最后有count段
            lack_count = alignment.seqA.count('-')  ##最后有count段
            if (lack_count in splice_count_dict.keys()):
                splice_count_dict[lack_count] = splice_count_dict[lack_count] + syboml_dic_dict[seq]
            else:
                splice_count_dict[lack_count] = syboml_dic_dict[seq]
            for s in seq1:
                if (s not in alignment.seqA):
                    # print(alignment.seqA)
                    key = str(lack_count) + '-' + s
                    if (key in pos_ndup_dict.keys()):
                        pos_ndup_dict[key] = pos_ndup_dict[key] + syboml_dic_dict[seq]
                    else:
                        pos_ndup_dict[key] = syboml_dic_dict[seq]

            break

        count_pos_ndup = count_pos_ndup + syboml_dic_dict[seq]
        print(len(seq), seq, syboml_dic_dict[seq], file=fpos_ndup)

    ##位置错乱且具有重复
    if (not pos and not ndup):
        count_npos_dup = count_npos_dup + syboml_dic_dict[seq]
        print(seq, syboml_dic_dict[seq], file=fnpos_dup)

    ##位置正确，有重复
    if (pos and not ndup):
        count_pos_dup = count_pos_dup + syboml_dic_dict[seq]
        print(seq, syboml_dic_dict[seq], file=fpos_dup)

print("全部条数", total_seq, file=fpos_ndup)
print("位置正确，重复", count_pos_dup, file=fpos_dup)
print("位置错乱，无重复", count_npos_ndup, file=fnpos_ndup)
print("位置错乱，重复", count_npos_dup, file=fnpos_dup)

print("位置正确，重复", count_pos_dup, file=fpos_ndup)
print("位置错乱，无重复", count_npos_ndup, file=fpos_ndup)
print("位置错乱，重复", count_npos_dup, file=fpos_ndup)
print("位置正确，无重复", count_pos_ndup, file=fpos_ndup)
lsn=[]
lsni=0
total_count=0
sorted_splice_count_list = sorted(splice_count_dict.items(), key=lambda item: item[0])
for elem in sorted_splice_count_list:
    print(elem[0], elem[1], file=fpos_ndup)
    lsn.append(elem[1])
    total_count=total_count+int(elem[1])
pp=plt.figure(dpi=120)
index=[1,2,3,4,5,6,7,8,9]
print(total_count)
print(index)
print(lsn)
#str1 = ("lack 1 seg", "上海", "武汉", "深圳", "重庆")#tick_label=str1
plt.bar(index, lsn[1:10],alpha=0.9,width=0.4)
for a, b in zip(index, lsn[1:10]):
    #plt.text(a, b + 0.05, '.1%' %(b/(total_count-lsn[6])*100), ha='center', va='bottom', fontsize=10)
    plt.text(a, b, str(round(b/(total_count-lsn[0])*100, 2))+"%", ha='center', va='bottom', fontsize=10)
plt.xticks(index, ('1', '2', '3', '4', '5', '6', '7', '8', '9'))
plt.title(nameID)
plt.xlabel('Number of missing segments')
plt.ylabel('Reads')
pp.savefig(nameID + "-count.pdf")



lcn = np.zeros((10,10))

lcdf = pd.DataFrame(lcn,
                  index=[ "0","1","2","3","4","5","6","7","8","9"],
                  columns = ['A', 'B', 'C', 'D', "E", "F", "G", "H", "I", "J"] )#

print(lcdf)
for key in pos_ndup_dict.keys():
    p = splice_count_dict[int(key.split('-')[0])]
    pp = pos_ndup_dict[key]
    sname=key.split('-')[0]
    cname=key.split('-')[1]
    lcrate=round(pp / p, 2)
    lcdf.loc[sname,cname]=lcrate
    print(key, pp, p, p - pp,lcrate , file=fpos_ndup)
lcdf.columns= ['index1', 'index2', 'index3',"index4",'elem1', 'elem2', 'elem3', 'elem4',"elem5","check"]


pp=plt.figure(dpi=120,figsize=(10,8))
print(lcdf.iloc[1:10,])
sns.heatmap(data=lcdf.iloc[1:10,],cmap=plt.get_cmap('Greens'),annot=True)
plt.title(nameID)
plt.ylabel('Number of missing segments',fontsize=12)
plt.xlabel('Missing segment location',fontsize=12) #x轴label的文本和字体大小
plt.xticks(fontsize=10,rotation=45) #x轴刻度的字体大小（文本包含在pd_data中了）
plt.yticks(fontsize=10) #y轴刻度的字体大小（文本包含在pd_data中了）
pp.savefig(nameID + ".pdf")


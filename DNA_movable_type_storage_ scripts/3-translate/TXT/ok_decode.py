# coding=utf-8
import sys
import numpy as np
from matplotlib import pyplot as plt
dic_seq_path=sys.argv[1]
temp_seq_path=sys.argv[2]
movable_word_path=sys.argv[3]
temp_seq=open(temp_seq_path,"r")
dic_seq=open(dic_seq_path,"r")
movable_word=open(movable_word_path,"r")
nameID=dic_seq_path.split('/')[-1].split('.')[0]
ok_dic=open(nameID+".okdic","w")
err_dic=open(nameID+".errdic","w")
ok_dict={}
err_dict={}
movable_word_dict={}
for word in movable_word.readlines():
    seq=word.strip().split('\t')[0]
    dic=word.strip().split('\t')[1]
    movable_word_dict[seq]=dic

def DNA_complement(sequence):  ###求互补
    complement = {'C': 'G', 'G': 'C', 'T': 'A', 'A': 'T','a':'t','t':'a','g':'c','c':'g'}
    res = ''
    for chr in sequence:
        res = res + complement[chr]
    return res


reads_seq_dict={}
reads_seq=temp_seq.readlines()
len_reads_seq=len(reads_seq)
i=0
while(i<len_reads_seq):
    reads_dic_tup = ()
    seq_elem_dict = {}#记录各个片段位置
    for seq in movable_word_dict.keys():
        if (seq in reads_seq[i]):
            seq_elem_dict[seq]=reads_seq[i].find(seq)
        elif (DNA_complement(seq)[::-1] in reads_seq[i]):
            seq_elem_dict[seq] = (reads_seq[i].find(DNA_complement(seq)[::-1])) * (-1)
    seq_elem_list=sorted(seq_elem_dict.items(), key=lambda item: item[1], reverse=True)

    for item in seq_elem_list:
        reads_dic_tup=(movable_word_dict[item[0]],)+reads_dic_tup

    ok_dict[reads_dic_tup] = 0

    i=i+1

for line in dic_seq.readlines():
    if("()" in line):
        print(line.strip('\n'),file=err_dic)
    else:
        line=line.replace('\'','')
        line_array=line.split(')')
        line_count=int(line_array[1].replace(',',''))
        line_tup_array=line_array[0].replace('(','').split(',')
        line_tup=()
        for elem in line_tup_array:
            line_tup=line_tup+(elem.strip(' '),)
        try:
            ok_dict[line_tup]=line_count+ok_dict[line_tup]
        except:
            if(line_tup in err_dict.keys()):
                err_dict[line_tup]=line_count+err_dict[line_tup]
            else:
                err_dict[line_tup] = line_count

            #print(line.strip('\n'),file=err_dic)


ok_total=0
for key in ok_dict.keys():
    ok_total=ok_total+ok_dict[key]
    print(ok_dict[key])
    print("(",key,",",ok_dict[key],")",file=ok_dic)
print(ok_total,file=ok_dic)
okme=np.zeros(len(ok_dict.keys()),dtype=int)
i=0
for key in ok_dict.keys():
    okme[i]=ok_dict[key]
    print(i,ok_dict[key],ok_total,okme[i])
    i=i+1

# 创建分组柱状图，需要自己控制x轴坐标
xticks = np.arange(len(ok_dict.keys()))
fig, ax = plt.subplots(figsize=(6, 4))
# 所有门店第一种产品的销量，注意控制柱子的宽度，这里选择0.25
ax.bar(xticks, okme, width=0.5, color="red")
ax.set_title("Correct reads Bar plot", fontsize=15)
ax.set_xlabel("ID")
ax.set_ylabel("Reads")
# 最后调整x轴标签的位置
#ax.set_xticks(xticks + 0.25)
#ax.set_xticklabels(range(1,len(ok_dict.keys())+1))
fig.savefig(nameID + "-okbarplot.pdf")


err_total=0
for key in err_dict.keys():
    err_total=err_total+err_dict[key]
    print("(",key,",",err_dict[key],")",file=err_dic)






#coding=utf8
import sys
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
dic_seq_path=sys.argv[1]
temp_seq_path=sys.argv[2]
movable_word_path=sys.argv[3]
temp_seq=open(temp_seq_path,"r")
dic_seq=open(dic_seq_path,"r")
movable_word=open(movable_word_path,"r")
nameID=dic_seq_path.split('/')[-1].split('.')[0]
err_dic=open(nameID+".arrangedic","w")
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

def comp(tem_list,seq):
    min_count = 11
    for tem in tem_list:
        count = 0
        for j in range(0,10):
            if(tem[j]!=seq[j]):
                count=count+1
        if(min_count>count):
            min_count=count

    arrange_pos_dic = {}  # 重组的位置
    i = 0
    for tem in tem_list:
        count = 0
        ap = []
        for j in range(0,10):
            if(tem[j]!=seq[j]):
                count=count+1
                ap.append(j)
        if(count==min_count):
            arrange_pos_dic[i]=ap
            #print(ap)
        i = i + 1

    return min_count,arrange_pos_dic

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

####从1-n单种和多种可能
sm=np.zeros((2,10),dtype=int)
print(sm)

####单种可能热图
si=np.zeros((10,10),dtype=int)
sirate=np.zeros((10,10),dtype=float)
print(si)
####多种可能热图
mu=np.zeros((10,10),dtype=int)
print(mu)

for key in err_dict.keys():
    if(len(key)!=10):
        continue
    min_count,arr=comp(ok_dict.keys(),key)#arr键值为最相似的模板序列编号
    print("(",key,",",err_dict[key],")",file=err_dic)
    print(min_count,arr,file=err_dic)
    bse = len(arr.keys())
    if(bse>1):
        sm[1][min_count-1]=sm[1][min_count-1]+err_dict[key]
    else:
        sm[0][min_count - 1] = sm[0][min_count - 1] + err_dict[key]
    for akey in arr.keys():
        for val in arr[akey]:
            if(bse==1):
                si[min_count - 1][val] = si[min_count - 1][val] + err_dict[key] * (1 / bse)
            else:
                mu[min_count - 1][val] = mu[min_count - 1][val] + err_dict[key] * (1 / bse)
for i in range(0,10):
    for j in range(0,10):
        if(sum(si[i])==0):
            continue
        sirate[i][j]=round(si[i][j]/sum(si[i]),2)

print(sm)
print(si)
print(sirate)
seg = [1, 2, 3, 4, 5,6,7,8,9,10]

# 创建分组柱状图，需要自己控制x轴坐标
xticks = np.arange(len(seg))

fig, ax = plt.subplots(figsize=(8, 6))
# 所有门店第一种产品的销量，注意控制柱子的宽度，这里选择0.25
ax.bar(xticks, sm[0], width=0.6, label="One", color="red",alpha=0.8)
# 所有门店第二种产品的销量，通过微调x轴坐标来调整新增柱子的位置
ax.bar(xticks, sm[1], width=0.6,bottom=sm[0], label="More", color="blue",alpha=0.8)
print(sm[0])
for a, b ,c in zip(seg, sm[0],sm[1]):
    #plt.text(a, b + 0.05, '.1%' %(b/(total_count-lsn[6])*100), ha='center', va='bottom', fontsize=10)
    if(b+c==0):
        continue
    ax.text(a-1-0.3, b+c+0.5, str(round(b/(c+b)*100,2))+"%",  fontsize=8)
    #ax.text(a-1+0.5, c/2+b, str(round(c / (c + b) * 100, 2)) + "%", fontsize=9)
ax.set_title("Recombination Bar plot", fontsize=15)
ax.set_xlabel("Number of recombinant segments")
ax.set_ylabel("Reads")
plt.ylim(0,max(sm[0])+50)
ax.legend()
# 最后调整x轴标签的位置
#ax.set_xticks(xticks + 0.25)
ax.set_xticks(xticks)
ax.set_xticklabels(seg)
fig.savefig(nameID + "-smarrange.pdf")


###########################创建热图
#psi,siax=plt.subplots(1,2,figsize=(6,4))

psi=plt.figure(figsize=(10,8),dpi=120)
lsi = pd.DataFrame(sirate,
                  index=[ "1","2","3","4","5","6","7","8","9","10"],#
                  columns=['index1', 'index2', 'index3',"index4",'elem1', 'elem2', 'elem3', 'elem4',"elem5","check"])#

print(lsi)
sns.heatmap(data=lsi,cmap=plt.get_cmap('OrRd'),annot=True)
plt.title(nameID)
plt.xlabel('Recombinant segment location',fontsize=12) #x轴label的文本和字体大小
plt.ylabel('Number of recombinant segments',fontsize=12) #y轴label的文本和字体大小
plt.xticks(fontsize=10,rotation=45) #x轴刻度的字体大小（文本包含在pd_data中了）
plt.yticks(fontsize=10) #y轴刻度的字体大小（文本包含在pd_data中了）
psi.savefig(nameID + "-siarrange.pdf")
#plt.title('title',fontsize=20) #图片标题文本和字体大小
"""###########################
#pmu,muax=plt.subplots(figsize=(6,4))
pmu=plt.figure(dpi=240)
lmu = pd.DataFrame(mu,
                  index=[ "1","2","3","4","5","6","7","8","9","10"],#
                   columns=['index1', 'index2', 'index3',"index4",'elem1', 'elem2', 'elem3', 'elem4',"elem5","check"])#

print(lmu)
sns.heatmap(data=lmu,cmap=plt.get_cmap('OrRd'),annot=True)
plt.title(nameID)
pmu.savefig(nameID + "-muarrange.pdf")
"""


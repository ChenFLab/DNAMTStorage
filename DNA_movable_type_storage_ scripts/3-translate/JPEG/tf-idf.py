#coding=utf-8
import math
import sys
from matplotlib import pyplot as plt
import  numpy as np
fallp=sys.argv[1]
fzerop=sys.argv[2]
total=sys.argv[3]
fall=open(fallp,"r")
fzero=open(fzerop,"r")
nameID = fzerop.split('/')[-1].split('.')[0]
word_fre={}
zero_list=[]
zero_fre={}
err={}
tf_idf={}
for line in fzero.readlines():
    if(line[0].isdigit()):
        line=line.strip('\n')
        line_array = line.strip().split(' ')
        if(line_array[0] not in zero_list):
            zero_list.append(line_array[0])
    else:
        continue



for line in fall.readlines():
    line_array = line.strip().split('\t')
    id=line_array[0]
    for item in line_array[1:]:
        if(id in zero_list):
            if(item in zero_fre.keys()):
                zero_fre[item]=zero_fre[item]+1
            else:
                zero_fre[item]=1

        if(item in word_fre.keys()):
            word_fre[item]=word_fre[item]+1
        else:
            word_fre[item]=1



for key in word_fre.keys():
    if(key not in zero_fre.keys()):
        zero_fre[key]=0
    tf_idf[key]=round(zero_fre[key]*math.log((int(total)-len(zero_list)+1)/(word_fre[key]-zero_fre[key]+1),2),2)##总文档/包含词条的文档
    err[key] = round(zero_fre[key]/word_fre[key]*100,2)



ferr=open(nameID+"-err.txt","w")
ftfidf=open(nameID+"-tfidf.txt","w")
ftf=open(nameID+"-zerotf.txt","w")
fttf=open(nameID+"-totaltf.txt","w")
for item in err.keys():
    print(item,err[item],file=ferr)
    print(item,tf_idf[item],file=ftfidf)
    print(item, word_fre[item], file=fttf)
    print(item, zero_fre[item], file=ftf)


err_dict_list = sorted(err.items(), key=lambda item: int(item[0]), reverse=False)
err_x=[]
err_y=[]
for item in err_dict_list:
    err_x.append(int(item[0]))
    err_y.append(float(item[1]))

print(err_dict_list)
xticks = np.arange(len(err_x))
fig, ax = plt.subplots(figsize=(8, 6))
# 所有门店第一种产品的销量，注意控制柱子的宽度，这里选择0.25
ax.bar(xticks, err_y, width=0.5, color="red",alpha=0.8)
# 所有门店第二种产品的销量，通过微调x轴坐标来调整新增柱子的位置
for a, b in zip(err_x, err_y):
    #plt.text(a, b + 0.05, '.1%' %(b/(total_count-lsn[6])*100), ha='center', va='bottom', fontsize=10)
    if(b<70):
        continue
    else:
        ax.text(a-0.25, b, b,  fontsize=9)
        ax.text(a -0.25, b+2, a, fontsize=9)

ax.set_title("movable type bar plot", fontsize=15)
ax.set_xlabel("movable type")
ax.set_ylabel("error rate(%)")
ax.legend()
# 最后调整x轴标签的位置
#ax.set_xticks(xticks + 0.25)
#ax.set_xticks(xticks)
#ax.set_xticklabels(err.keys())
fig.savefig(nameID + "-mwerr.pdf")


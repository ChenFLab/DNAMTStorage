#coding=utf-8
import sys
import math
em_path=sys.argv[1]
fm_path=sys.argv[2]

fem=open(em_path,"r")
ffm=open(fm_path,"r")
##
ID=1
em_dict={}
for line in fem.readlines():
    line=line.strip()
    line=line.replace(')','')
    em_dict[ID]=int(line.split(',')[-1])
    ID=ID+1



ID=1
fm_dict={}
for line in ffm.readlines():
    line=line.strip()
    line=line.replace(')','')
    fm_dict[ID]=int(line.split(',')[-1])
    ID=ID+1


merge_dict={}
merge_dict_path=""
for i in em_dict.keys():
    merge_dict[i]=em_dict[i]+fm_dict[i]

merfile_path=sys.argv[3]
fmerfile=open(merfile_path,"w")
fzero=open("zero"+merfile_path,'w')

zero_dict={}
groupno_dict={}
for item in merge_dict.keys():
    innerno=item%30
    if(innerno==0):
        innerno=30
    groupno=((item-innerno)//30)+1
    print(item,groupno,innerno,merge_dict[item],file=fmerfile)
    if(merge_dict[item]==0):
        print(item, groupno, innerno, merge_dict[item], file=fzero)
        if(groupno not in zero_dict.keys()):
            zero_dict[groupno]=1
        else:
            zero_dict[groupno]=zero_dict[groupno]+1

total_uncover=0
for gno in zero_dict.keys():
    if(zero_dict[gno]>=11):
        print("zero gno:",gno,"zero count:",zero_dict[gno],file=fzero)
        ino=1
        cannot_rev=0
        while(ino<=20):
            if(merge_dict[(gno-1) * 30+ino]==0):
                cannot_rev=cannot_rev+1
                print()
            ino=ino+1
        total_uncover=total_uncover+cannot_rev
        print("zero gno:",gno,"cannot_recover:", cannot_rev, file=fzero)
print("total_unrecover:",total_uncover,file=fzero)

fzero_uncover=open("uncover"+merfile_path,"w")
flag_uncover={}
for item in merge_dict.keys():###未恢复的冗余还是在里面
    innerno=item%30
    if(innerno==0):
        innerno=30
    groupno=((item-innerno)//30)+1

    if(merge_dict[item]==0 and zero_dict[groupno]>10):
        print(item,groupno,innerno,merge_dict[item],int(0),file=fzero_uncover)
    elif(merge_dict[item]==0 and zero_dict[groupno]<=10):
        print(item, groupno, innerno, merge_dict[item], int(1), file=fzero_uncover)
    else:
        print(item, groupno, innerno, merge_dict[item], int(1), file=fzero_uncover)










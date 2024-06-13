# coding=utf-8
import sys
dic_path1=sys.argv[1]
dic_path2=sys.argv[2]
merge_dic_path=sys.argv[3]
with open (dic_path1,"r") as fdic1,open (dic_path2,"r") as fdic2:
    dic1=fdic1.readlines()
    dic2=fdic2.readlines()
dicm=dic1+dic2
merge_dic={}
for line in dicm:
    line = line.replace('\n','')
    line = line.replace(' ','')
    if("()" in line):
        line = line.replace('(', '')
        line = line.replace(')', '')
        line_array = line.split(',')
        line_count=int(line_array[-1])
        if("()" in merge_dic.keys()):
            merge_dic["()"]=merge_dic["()"]+line_count
        else:
            merge_dic["()"]=line_count
    else:
        line = line.replace('(', '')
        line = line.replace(')', '')
        line=line.replace('\'','')
        line_array=line.split(',')
        line_count=int(line_array[-1])
        line_tup_array=line_array[:-1]
        line_tup=()
        for elem in line_tup_array:
            line_tup=line_tup+(elem.strip(' '),)

        if(line_tup in merge_dic.keys()):
            merge_dic[line_tup]=line_count+merge_dic[line_tup]
        else:
            merge_dic[line_tup] = line_count

fmerge_dic=open(merge_dic_path,"w")
for key in merge_dic.keys():
    print("(",key,",",merge_dic[key],")",file=fmerge_dic)









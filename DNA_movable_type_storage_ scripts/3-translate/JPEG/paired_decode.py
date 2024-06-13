import sys
import operator
import time

reads1_seq_path = sys.argv[1]
reads2_seq_path = sys.argv[2]
mw_path =sys.argv[3]
max_splice=sys.argv[4]
reads1_seq_path.split('/')
nameID=reads1_seq_path.split('/')[-1].split('.')[0]
em_res_path =nameID+"-em.res"
snp_path=nameID+"-snp.fa"
def DNA_complement(sequence):  ###求互补
    complement = {'C': 'G', 'G': 'C', 'T': 'A', 'A': 'T','a':'t','t':'a','g':'c','c':'g'}
    res = ''
    for chr in sequence:
        res = res + complement[chr]
    return res

####创建活字单元字典wh_temp_dict

def gen_next(s2):
    k = -1
    n = len(s2)
    j = 0
    next_list = [0 for i in range(n)]
    next_list[0] = -1                           #next数组初始值为-1
    while j < n-1:
        if k == -1 or s2[k] == s2[j]:
            k += 1
            j += 1
            next_list[j] = k                    #如果相等 则next[j+1] = k
        else:
            k = next_list[k]                    #如果不等，则将next[k]的值给k
    return next_list


def match(s1, s2, next_list):
    ans = -1
    i = 0
    j = 0
    while i < len(s1):
        if s1[i] == s2[j] or j == -1:
            i += 1
            j += 1
        else:
            j = next_list[j]
        if j == len(s2):
            ans = i - len(s2)
            break
    return ans







f_mw=open(mw_path,'r')
f_snp=open(snp_path,"w")
mw_lines=f_mw.readlines()
f_mw.close()
mw_dict={}
for word_line in mw_lines:
    word_line=word_line.strip('\n')
    mw_dict[word_line.split('\t')[0]]=word_line.split('\t')[1]

####创建reads1字典
reads1_dict={}
f_reads1_seq=open(reads1_seq_path,"r")
reads1_seq_lines=f_reads1_seq.readlines()
f_reads1_seq.close()
len_reads1=len(reads1_seq_lines)
i=0
while(i<len_reads1):
    reads1_dict[reads1_seq_lines[i].split(' ')[0].replace('@','>')]=reads1_seq_lines[i+1].strip('\n')
    i=i+4

####创建reads2字典
reads2_dict={}
f_reads2_seq=open(reads2_seq_path,"r")
reads2_seq_lines=f_reads2_seq.readlines()
f_reads2_seq.close()
len_reads2=len(reads2_seq_lines)
j=0
while(j<len_reads2):
    reads2_dict[reads2_seq_lines[j].split(' ')[0].replace('@','>')]=reads2_seq_lines[j+1].strip('\n')
    j=j+4

####分别解码并比对,统计正确reads数目，丢弃匹配不一致的序列，过滤出含有缺失值的序列,成对输出
f_em_res=open(em_res_path,"w")
f_snp=open(snp_path,"w")
em_dict={}
err_count=0
total=0

for read1_key in reads1_dict.keys():
    total=total+1
    read1=reads1_dict[read1_key]
    read2=reads2_dict[read1_key]
    read1_pos_dict={}###用来存储reads中的片段位置
    read2_pos_dict={}
    for word in mw_dict.keys():

        next_list = gen_next(word)

        if(read1.count(word)==1):
            # startk=time.perf_counter()
            # print("kmp", match(read1, word, next_list))
            # endk=time.perf_counter()
            # print(endk-startk)

            # startf=time.perf_counter()
            # print("find", read1.find(word))
            # endf=time.perf_counter()
            # print(endf-startf)

            # read1_pos_dict[mw_dict[word]]=read1.find(word)
            read1_pos_dict[mw_dict[word]]=match(read1, word, next_list)
        if(read2.count(word)==1):
            # print("kmp", match(read1, word, next_list))
            # print("find", read1.find(word))
            # read2_pos_dict[mw_dict[word]]=read2.find(word)
            read2_pos_dict[mw_dict[word]]=match(read2, word, next_list)

    for word in mw_dict.keys():

        c_word=DNA_complement(word)[::-1]
        next_list_c = gen_next(c_word)

        if (read1.count(c_word) == 1):
            # read1_pos_dict[mw_dict[word]] = (-1) * read1.find(c_word)
            read1_pos_dict[mw_dict[word]] = (-1) * match(read1, c_word, next_list_c)
        if (read2.count(c_word) == 1):
            # read2_pos_dict[mw_dict[word]] = (-1) * read2.find(c_word)
            read2_pos_dict[mw_dict[word]] = (-1) * match(read2, c_word, next_list_c)



    read1_dict_list = sorted(read1_pos_dict.items(), key=lambda item: item[1], reverse=False)
    read2_dict_list = sorted(read2_pos_dict.items(), key=lambda item: item[1], reverse=False)
    read1_list=[]
    read2_list=[]
    for itm1 in read1_dict_list:
        read1_list.append(itm1[0])
    for itm2 in read2_dict_list:
        read2_list.append(itm2[0])
 
    if(len(read1_list)==int(max_splice) and len(read2_list)==int(max_splice)):##判断长度是否不缺失
        if(operator.eq(read1_list,read2_list)): ##判断不缺失情况下内容是否一致
            read1_list_tup=tuple(read1_list)
            if(read1_list_tup in em_dict.keys()):
                em_dict[read1_list_tup]=em_dict[read1_list_tup]+2
            else:
                em_dict[read1_list_tup]=2
        else:
            err_count=err_count+2
    else:##进行缺失值需要snp矫正，snp矫正之后再通过overlap矫正
        print(read1_key + "-1", file=f_snp)
        print(reads1_dict[read1_key], file=f_snp)
        print(read1_key + "-2", file=f_snp)
        print(reads2_dict[read1_key], file=f_snp)

for key in em_dict.keys():
    print(key,em_dict[key],file=f_em_res)
print(total,err_count)
f_em_res.close()
f_snp.close()
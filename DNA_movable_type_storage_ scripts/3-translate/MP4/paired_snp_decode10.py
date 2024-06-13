####序列比对
####模糊匹配中全对的
import sys
snp_sam_path=sys.argv[1]
mw_file_path=sys.argv[2]
max_splice=sys.argv[3]
nameID=snp_sam_path.split('/')[-1].split('.')[0]
snp_decode_path=nameID+".snp_decode"
fsnp_sam=open(snp_sam_path,"r")
fmw_file=open(mw_file_path,"r")##与sam比对的mw文件
mw_len=len(fmw_file.readlines())/2
print(mw_len)
fsnp_decode=open(snp_decode_path,"w")
fm_dict={}##用来存放模糊匹配成功的序列
err_count=0###统计两条不匹配条数
seq_count=0



i =0
snp_sam=fsnp_sam.readlines()
len_snp_sam=len(snp_sam)
while(i<len_snp_sam):

    j=0##两条reads区域内指针
    read1_tup_pos={}
    read2_tup_pos={}
    read1_list=[]
    read2_list=[]
    lost= False
    duplicate = False
    unmatch =False
    while(j<(mw_len+2) and int(i+j)<len_snp_sam):##read1

        line = snp_sam[int(i + j)].split('\t')
        if (line[1] == '0'):
            try:
                dup = line[11].split(':')[2]
                if (dup == "R"):
                    read1_tup_pos[line[0] + "R"] = int(line[3])
                    duplicate=True
                else:
                    read1_tup_pos[line[0]] = int(line[3])
            except:
                read1_tup_pos[line[0]] = int(line[3])

        if (line[1] == '16'):
            try:
                dup = line[11].split(':')[2]
                if (dup == "R"):
                    read1_tup_pos[line[0] + "R"] = -1 * int(line[3])
                    duplicate = True
                else:
                    read1_tup_pos[line[0]] = -1 * int(line[3])
            except:
                read1_tup_pos[line[0]] = -1 * int(line[3])
        read1_dict_list = tuple(sorted(read1_tup_pos.items(), key=lambda item: item[1], reverse=False))
        read1_list = []

        for itm1 in read1_dict_list:
            read1_list.append(itm1[0])
        

        #print(j, snp_sam[int(i + j)])
        j=j+1

    while(j<2*(mw_len+2)and int(i+j)<len_snp_sam):##read2
        line = snp_sam[int(i + j)].split('\t')
        if (line[1] == '0'):
            try:
                dup = line[11].split(':')[2]
                if (dup == "R"):
                    read2_tup_pos[line[0] + "R"] = int(line[3])
                    duplicate = True
                else:
                    read2_tup_pos[line[0]] = int(line[3])
            except:
                read2_tup_pos[line[0]] = int(line[3])

        if (line[1] == '16'):
            try:
                dup = line[11].split(':')[2]
                if (dup == "R"):
                    read2_tup_pos[line[0] + "R"] = -1 * int(line[3])
                    duplicate = True
                else:
                    read2_tup_pos[line[0]] = -1 * int(line[3])
            except:
                read2_tup_pos[line[0]] = -1 * int(line[3])
        read2_dict_list = tuple(sorted(read2_tup_pos.items(), key=lambda item: item[1], reverse=False))
        read2_list = []
        for itm2 in read2_dict_list:
            read2_list.append(itm2[0])
        
        j=j+1

    if (len(read1_list) == int(max_splice) and len(read2_list) == int(max_splice) and duplicate==False and unmatch==False):  ##判断长度是否一致,且不缺失
        if (read1_list == read2_list):  ##判断长度一致情况下内容是否一致
            read1_list_tup = tuple(read1_list)
            if (read1_list_tup in fm_dict.keys()):
                fm_dict[read1_list_tup] = fm_dict[read1_list_tup] + 2
            else:
                fm_dict[read1_list_tup] = 2
        else:
            unmatch=True
            err_count=err_count+2

    elif((len(read1_list)<= int(max_splice) or len(read2_list) <= int(max_splice)) and duplicate==False and unmatch==False):  ##进行缺失值的初步overlap矫正
        target_read_dict = {"index1": 0, "index2": 0, "index3": 0, "index4": 0, "elem1": 0, "elem2": 0, "elem3": 0, "elem4": 0, "elem5": 0, "check": 0}
        if(len(target_read_dict)!=int(max_splice)):
            print("target read_dict length error")
        read1_dict = {}
        read2_dict = {}
        for elem in read1_list:

            read1_dict[elem.split('-')[0]] = elem.split('-')[1]
        for elem in read2_list:
            read2_dict[elem.split('-')[0]] = elem.split('-')[1]

        for key in target_read_dict:
            if (key in read1_dict.keys() and key in read2_dict.keys()):
                if(read1_dict[key] == read2_dict[key]):
                    target_read_dict[key] = read1_dict[key]
                else:
                    unmatch=True
            elif (key not in read1_dict.keys() and key in read2_dict.keys() ):
                target_read_dict[key] = read2_dict[key]
            elif (key in read1_dict.keys() and key not in read2_dict.keys() ):
                target_read_dict[key] = read1_dict[key]
            else:
                target_read_dict[key]=-1
                lost=True

        if (unmatch == True or duplicate==True):
            err_count = err_count + 2
        if(unmatch==False and duplicate==False):
            read_tup = ()
            for key in target_read_dict.keys():
                if(target_read_dict[key]!=-1):
                    read_tup = read_tup + (key + "-" + str(target_read_dict[key]),)
            if (read_tup in fm_dict.keys()):
                fm_dict[read_tup] = fm_dict[read_tup] + 2
            else:
                fm_dict[read_tup] = 2
    else:
        err_count = err_count + 2

    i=int(i+2*(mw_len+2))


fm_dict_list=sorted(fm_dict.items(),key=lambda item:item[1],reverse=True)
total=0
for items in fm_dict_list:
    total=total+int(items[1])
    print(items,file=fsnp_decode)
print(seq_count,"err:",err_count,"no_err:",total)



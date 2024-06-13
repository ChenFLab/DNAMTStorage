import encode_decode.getmap as getmap
import numpy as np
import encode_decode.erasure_code as ecode

import os
import random
index1_size = 64#mix 4*16
index2_size = 32
index3_size = 32
index4_size = 32
elem1_size = 64
elem2_size = 64
elem3_size = 64
elem4_size = 64
elem5_size = 64
row_check  = 64

def complementary(seq):
    l = len(seq)
    i = 0
    s = ''
    while (i < l):
        if (seq[i] == 'A'):
            s = s + 'T'
        elif (seq[i] == 'T'):
            s = s + 'A'
        elif (seq[i] == 'C'):
            s = s + 'G'
        elif(seq[i] == 'G'):
            s = s + 'C'
        elif (seq[i] == 't'):
            s = s + 'a'
        elif (seq[i] == 'c'):
            s = s + 'g'
        elif(seq[i] == 'g'):
            s = s + 'c'
        elif (seq[i] == 'a'):
            s = s + 't'
        i = i + 1
    return s

def get_inputfile_to_user(in_path):
    file = open(in_path, 'rb')
    buffer = file.read()
    content_list = []
    map_b_to_c = []
    count = 0
    for i in buffer:
        # by = bytes(i, encoding='utf8')
        # for j in by:
        print(chr(i))
        content_list.append(i)
        map_b_to_c.append(count)
        count += 1
    file.close()
    return content_list, buffer, map_b_to_c

def en_code_2216_10(in_path, k=6, m=2, k_r=3, m_r=1, size=2, size_r=1):
    ftoexcel=in_path[:-4] + '2216_10__v3excel.txt'
    ftojet=in_path[:-4] + '2216_10__v3jet.txt'
    ftoranid=in_path[:-4] + '2216_10__v3map.txt'

    fmap = "//Users//freya//PycharmProjects//DNAStorage20211108//encode_decode//map//221610mw.txt"
    toexcel=open(ftoexcel,"w")
    tojet = open(ftojet, "w")
    toranid = open(ftoranid, "w")


    fdirector = in_path[:-4] + '_v3连接指导文件.txt'
    fodirector = in_path[:-4] + '_v3origin连接指导文件.txt'
    director = open(fdirector, "w")
    odirector=open(fodirector,"w")
    outpath = in_path[:-4] + '_seq_v3.txt'
    fdw = open(outpath, 'w')

    # try:
    content_list, buffer, map_b_to_c = get_inputfile_to_user(in_path)
    # except:
    #    print("input_error!")
    #    return "input_error!", 0
    map1, map2, map3, map4, map5, map6, map7,map8,map9,map10 = getmap.getallmap_2216_10(fmap)

    content_list_size = len(content_list)
    if (content_list_size > 20000000):
        print("warning! the file size is greater than 20MB")
    # 第一组记录文件的大小
    # size_in_fileseq(fdw)

    extended_matrix_r = ecode.gf_gen_rs_matrix(k_r, m_r)
    extended_matrix = ecode.gf_gen_rs_matrix(k, m)


    g_tbls = ecode.ec_init_tables(k, m, extended_matrix[k:k + m])
    g_tbls_r = ecode.ec_init_tables(k_r, m_r, extended_matrix_r[k_r:k_r + m_r])

    count=0
    while (count < content_list_size):  ####分组k+m
        en_buf = np.zeros((k + m, size), dtype=int)  ##对每一组进行编码
        for i in range(k):
            if (count < content_list_size - 3):
                en_buf[i][0] = content_list[count]
                en_buf[i][1] = content_list[count + 1]
                en_buf[i][2] = content_list[count + 2]
                en_buf[i][3] = content_list[count + 3]
                #print(en_buf[i][0], en_buf[i][1],en_buf[i][2],en_buf[i][3])
                #print(chr(en_buf[i][0]),chr(en_buf[i][1]),chr(en_buf[i][2]),chr(en_buf[i][3]))
                count += 4
            elif (count == content_list_size - 3):
                en_buf[i][0] = content_list[count]
                en_buf[i][1] = content_list[count + 1]
                en_buf[i][2] = content_list[count + 2]
                en_buf[i][3] = 0
                count += 4
            elif (count == content_list_size - 2):
                en_buf[i][0] = content_list[count]
                en_buf[i][1] = content_list[count + 1]
                en_buf[i][2] = 0
                en_buf[i][3] = 0
                count += 4
            elif (count == content_list_size - 1):
                en_buf[i][0] = content_list[count]
                en_buf[i][1] = 0
                en_buf[i][2] = 0
                en_buf[i][3] = 0
                count += 4
            else:
                count += 4

        col_check = ecode.ec_encode_data(size, k, m, g_tbls, en_buf)
        for i in range(m):
            en_buf[k + i][0] = col_check[i][0]
            en_buf[k + i][1] = col_check[i][1]
            en_buf[k + i][2] = col_check[i][2]
            en_buf[k + i][3] = col_check[i][3]


        for i in range(k + m):
            oelem1 = en_buf[i][0] % 64
            oelem2 = en_buf[i][0] // 64 + en_buf[i][1] % 16 * 4
            oelem3 = en_buf[i][1] // 16 + en_buf[i][2] % 4 * 16
            oelem4 = en_buf[i][2] // 4
            oelem5 = en_buf[i][3] % 64
            group_no = (count - 1) // (4 * k)

            index4 = i
            index3=group_no%index3_size
            group_no1=group_no//index3_size
            oindex2=group_no1%index2_size
            group_no2=group_no1//index2_size
            oindex1=group_no2%(index1_size//4)*4 +en_buf[i][3]//64
            #print(group_no,group_no % (index1_size//4),oindex2,index3)

            en_buf_r = np.zeros((k_r + m_r, size_r), dtype=int)

            elem1 = (oelem1 & 0x1f) + (oindex1 & 0x20)
            elem2 = (oelem2 & 0x1f) + (oindex1 & 0x10) * 2
            elem3 = (oelem3 & 0x1f) + (oindex1 & 0x08) * 4
            elem4 = (oelem4 & 0x1f) + (oindex2 & 0x10) * 2
            elem5 = (oelem5 & 0x1f) + (oindex2 & 0x08) * 4
            index1 = (oelem1 & 0x20) + (oelem2 & 0x20) // 2 + (oelem3 & 0x20) // 4 + (oindex1 & 0x07)
            index2 = (oelem4 & 0x20) // 2 + (oelem5 & 0x20) // 4 + (oindex2 & 0x07)
            en_buf_r[0] = index1
            en_buf_r[1] = index2
            en_buf_r[2] = index3
            en_buf_r[3] = index4
            en_buf_r[4] = elem1
            en_buf_r[5] = elem2
            en_buf_r[6] = elem3
            en_buf_r[7] = elem4
            en_buf_r[8] = elem5

            row_check = ecode.ec_encode_data(size_r, k_r, m_r, g_tbls_r, en_buf_r)[0,0]
            en_buf_r[9] = row_check
            #####exchange

            #print( index1, index2, index3, index4, elem1, elem2, elem3, elem4, elem5, row_check,file=director)

            oen_buf_r = np.zeros((k_r + m_r, size_r), dtype=int)
            oen_buf_r[0] = oindex1
            oen_buf_r[1] = oindex2
            oen_buf_r[2] = index3
            oen_buf_r[3] = index4
            oen_buf_r[4] = oelem1
            oen_buf_r[5] = oelem2
            oen_buf_r[6] = oelem3
            oen_buf_r[7] = oelem4
            oen_buf_r[8] = oelem5
            orow_check = ecode.ec_encode_data(size_r, k_r, m_r, g_tbls_r, en_buf_r)[0, 0]
            oen_buf_r[9] = orow_check

            cont1=chr(en_buf[i][0])
            cont2=chr(en_buf[i][1])
            cont3=chr(en_buf[i][2])
            cont4=chr(en_buf[i][3])
            if (cont1 == '\n' ):
                cont1 = "huanhang"
            if (cont1 == '\r' ):
                cont1 = "huiche"

            if (cont2 == '\n'):
                cont3 = "huanhang"
            if (cont2 == '\r' ):
                cont2 = "huiche"

            if (cont3 == '\n'):
                cont3 = "huanhang"
            if (cont3 == '\r' ):
                cont3 = "huiche"

            if (cont4 == "\n"):
                cont4 = "huanhang"
            if (cont4 == '\r' ):
                cont4 = "huiche"
            if(i<20):
                print(index1, index2, index3,index4, elem1, elem2, elem3, elem4, elem5,row_check,file=director)

                #print( oindex1-en_buf[i][3]//64, oindex2, index3, index4,en_buf[i][0],en_buf[i][1],en_buf[i][2],en_buf[i][3],file=odirector)
                print( "block-id:",bin(group_no),"inner-id:",bin(i),"binary:",bin(en_buf[i][0]) ,bin(en_buf[i][1]),bin(en_buf[i][2]),bin(en_buf[i][3]),"content:",cont1,cont2,cont3,cont4,file=odirector)
            else:
                #print(oindex1-en_buf[i][3]//64, oindex2, index3, index4, en_buf[i][0] , en_buf[i][1],en_buf[i][2],en_buf[i][3],file=odirector)
                print( "block-id:",bin(group_no),"inner-id:",bin(i),"binary:",bin(en_buf[i][0]) ,bin(en_buf[i][1]),bin(en_buf[i][2]),bin(en_buf[i][3]),file=odirector)
                print(index1,index2, index3, index4, elem1, elem2, elem3, elem4, elem5,row_check,file=director)

            #print("change: ", index1, index2, index3, index4, elem1, elem2, elem3, elem4, elem5, row_check,file=director)
            #print("ori: ",oindex1,oindex2,index3,index4,oelem1,oelem2,oelem3,oelem4,oelem5,row_check,file=director)
            #print(group_no, i, en_buf[i][0], en_buf[i][1], en_buf[i][2], en_buf[i][3], file=director)

            #1 elem1
            left = 'aagctt'
            right = 'ctcc'
            insert1_insert_seqF = map1[index1]
            insert1_insert_seqR = complementary(insert1_insert_seqF)[::-1]
            insert1_seqF = "agctt" + insert1_insert_seqF
            insert1_seqR = complementary(right)[::-1] + insert1_insert_seqR+ 'a'
            insert1_id = "{:0>3}".format(index1)
            insert1F_info = insert1_id + "-1F " + insert1_seqF + " " + str(len(insert1_seqF))
            insert1R_info = insert1_id + "-1R " + insert1_seqR + " " + str(len(insert1_seqR))

            #2 elem2
            left = 'ctcc'
            right = 'acta'
            insert2_insert_seqF = map2[index2]
            insert2_insert_seqR = complementary(insert2_insert_seqF)[::-1]
            insert2_seqF = left + insert2_insert_seqF
            insert2_seqR = complementary(right)[::-1] + insert2_insert_seqR
            insert2_id = "{:0>3}".format(index2)
            insert2F_info = insert2_id + "-2F " + insert2_seqF + " " + str(len(insert2_seqF))
            insert2R_info = insert2_id + "-2R " + insert2_seqR + " " + str(len(insert2_seqR))

            #3 elem3
            left = 'acta'
            right = 'ggta'
            insert3_insert_seqF = map3[index3]
            insert3_insert_seqR = complementary(insert3_insert_seqF)[::-1]
            insert3_seqF = left + insert3_insert_seqF
            insert3_seqR = complementary(right)[::-1] + insert3_insert_seqR
            insert3_id = "{:0>3}".format(index3)
            insert3F_info = insert3_id + "-3F " + insert3_seqF + " " + str(len(insert3_seqF))
            insert3R_info = insert3_id + "-3R " + insert3_seqR + " " + str(len(insert3_seqR))

            #4 elem4
            left = 'ggta'
            right = 'atag'
            insert4_insert_seqF = map4[index4]
            insert4_insert_seqR = complementary(insert4_insert_seqF)[::-1]
            insert4_seqF = left + insert4_insert_seqF
            insert4_seqR = complementary(right)[::-1] + insert4_insert_seqR
            insert4_id = "{:0>3}".format(index4)
            insert4F_info = insert4_id + "-4F " + insert4_seqF + " " + str(len(insert4_seqF))
            insert4R_info = insert4_id + "-4R " + insert4_seqR + " " + str(len(insert4_seqR))

            #5 index1
            left = 'atag'
            right = 'aagt'
            insert5_insert_seqF = map5[elem1]
            insert5_insert_seqR = complementary(insert5_insert_seqF)[::-1]
            insert5_seqF = left + insert5_insert_seqF
            insert5_seqR = complementary(right)[::-1] + insert5_insert_seqR
            insert5_id = "{:0>3}".format(elem1)
            insert5F_info = insert5_id + "-5F " + insert5_seqF + " " + str(len(insert5_seqF))
            insert5R_info = insert5_id + "-5R " + insert5_seqR + " " + str(len(insert5_seqR))

            #6 check1
            left = 'aagt'
            right = 'ctga'
            insert6_insert_seqF = map6[elem2]
            insert6_insert_seqR = complementary(insert6_insert_seqF)[::-1]
            insert6_seqF = left + insert6_insert_seqF
            insert6_seqR = complementary(right)[::-1] + insert6_insert_seqR
            insert6_id = "{:0>3}".format(elem2)
            insert6F_info = insert6_id + "-6F " + insert6_seqF + " " + str(len(insert6_seqF))
            insert6R_info = insert6_id + "-6R " + insert6_seqR + " " + str(len(insert6_seqR))

            #7 elem3
            left = 'ctga'
            right = 'caac'
            insert7_insert_seqF = map7[elem3]
            insert7_insert_seqR = complementary(insert7_insert_seqF)[::-1]
            insert7_seqF = left + insert7_insert_seqF
            insert7_seqR = complementary(right)[::-1] + insert7_insert_seqR
            insert7_id = "{:0>3}".format(elem3)
            insert7F_info = insert7_id + "-7F " + insert7_seqF + " " + str(len(insert7_seqF))
            insert7R_info = insert7_id + "-7R " + insert7_seqR + " " + str(len(insert7_seqR))

            #8 index4
            left = 'caac'
            right = 'agga'
            insert8_insert_seqF = map8[elem4]
            insert8_insert_seqR = complementary(insert8_insert_seqF)[::-1]
            insert8_seqF = left + insert8_insert_seqF
            insert8_seqR = complementary(right)[::-1] + insert8_insert_seqR
            insert8_id = "{:0>3}".format(elem4)
            insert8F_info = insert8_id + "-8F " + insert8_seqF + " " + str(len(insert8_seqF))
            insert8R_info = insert8_id + "-8R " + insert8_seqR + " " + str(len(insert8_seqR))

            #9 index2
            left = 'agga'
            right = 'tgtg'
            insert9_insert_seqF = map9[elem5]
            insert9_insert_seqR = complementary(insert9_insert_seqF)[::-1]
            insert9_seqF = left + insert9_insert_seqF
            insert9_seqR = complementary(right)[::-1] + insert9_insert_seqR
            insert9_id = "{:0>3}".format(elem5)
            insert9F_info = insert9_id + "-9F " + insert9_seqF + " " + str(len(insert9_seqF))
            insert9R_info = insert9_id + "-9R " + insert9_seqR + " " + str(len(insert9_seqR))

            #10 check2
            left = 'tgtg'
            right = 'gaattc'

            insert10_insert_seqF = map10[row_check]
            insert10_insert_seqR = complementary(insert10_insert_seqF)[::-1]
            insert10_seqF = left + insert10_insert_seqF+ "g"
            insert10_seqR = 'aattc' + insert10_insert_seqR
            insert10_id = "{:0>3}".format(row_check)
            insert10F_info = insert10_id + "-10F " + insert10_seqF + " " + str(len(insert10_seqF))
            insert10R_info = insert10_id + "-10R " + insert10_seqR + " " + str(len(insert10_seqR))

            #print("change ",index1, index2, index3, index4,elem1,elem2,elem3,elem4,elem5,row_check,file=director)

            seq_excel = insert1F_info + " " + insert2F_info + " " + insert3F_info + " " + insert4F_info +\
                        " " + insert5F_info + " " + insert6F_info + " " + insert7F_info + " " + insert8F_info +\
                        " " + insert9F_info + " " + insert10F_info + "\n" +\
                        insert1R_info + " " + insert2R_info + " " + insert3R_info + " " + insert4R_info +\
                        " " + insert5R_info + " " + insert6R_info + " " + insert7R_info + " " + insert8R_info+\
                        " " + insert9R_info + " " + insert10R_info+'\n'
            seq_excel_content=str(hex(index1))+ " " +str(hex(index2))+" " +str(hex(index3))+" " +str(hex(index4))+" " +str(hex(elem1))+" " \
                              +str(hex(elem2))+" " +str(hex(elem3))+" " + str(hex(elem4))+" " +str(hex(elem5))+" " +str(hex(row_check))+"\n"
            toexcel.write(seq_excel_content)
            toexcel.write(seq_excel)
            seq ="a"+insert1_seqF + insert2_seqF + insert3_seqF + insert4_seqF + insert5_seqF + \
                  insert6_seqF + insert7_seqF + insert8_seqF + insert9_seqF + insert10_seqF+"aattc"
            fdw.write(seq + "\n")

            jets=insert1_id + "-1"+" " + insert2_id + "-2"+" " + insert3_id + "-3"+" " + insert4_id + "-4"\
                 +" " + insert5_id + "-5"+" " + insert6_id + "-6"+" " + insert7_id + "-7"+" " + insert8_id + "-8"\
                 +" " + insert9_id + "-9"+" " + insert10_id + "-10"+'\n'
            tojet.write(jets)


    fdw.close()
    toexcel.close()
    tojet.close()
    id_list=[]
    tojet = open(ftojet, "r")
    for line in tojet.readlines():
        id_list.append(line.strip('\n'))

    index = [i for i in range(1, len(id_list)+1)]
    random.shuffle(index)#changed
    i = 0
    print("randomid","id", file=toranid)
    for elem in index:
        print(i+1,elem,id_list[elem-1],file=toranid)
        i = i + 1

    return "OK", content_list_size

if __name__ == '__main__':
    os.chdir("//Users//freya//PycharmProjects//DNAStorage20211108")
    flag,file_size = en_code_2216_10(
        "encode_decode//data/2216-10//panda.jpg", k=20,
        m=10, k_r=9, m_r=1, size=4, size_r=1)
    print(file_size)

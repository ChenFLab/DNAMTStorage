import encode_decode.erasure_code as ecode
import numpy as np
import math
import sys
import os

max_groupno = 19200
file_group_size = 7458


def getallmap(read_map):
    map1 = []
    map2 = []
    map3 = []
    map4 = []
    map5 = []
    map6 = []
    map7 = []
    map8 = []
    map9 = []
    map10 = []

    fd = open(read_map, "r")
    count = 0
    line = fd.readline()
    while line != '':
        mark1 = line.strip('\n').split(' ')[1].split('-')[0]
        content = line.strip('\n').split(' ')[0]
        if mark1 == 'mw1':
            count += 1
            map1.append(content)
        elif mark1 == 'mw2':
            count += 1
            map2.append(content)
        elif mark1 == 'mw3':
            count += 1
            map3.append(content)
        elif mark1 == 'mw4':
            count += 1
            map4.append(content)
        elif mark1 == 'mw5':
            count += 1
            map5.append(content)
        elif mark1 == 'mw6':
            count += 1
            map6.append(content)
        elif mark1 == 'mw7':
            count += 1
            map7.append(content)
        elif mark1 == 'mw8':
            count += 1
            map8.append(content)
        elif mark1 == 'mw9':
            count += 1
            map9.append(content)
        elif mark1 == 'mw10':
            count += 1
            map10.append(content)
        line = fd.readline()
    return map1, map2, map3, map4, map5, map6, map7, map8, map9, map10


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
        else:
            s = s + 'C'
        i = i + 1
    return s[::-1]


def unit_from_overlap(overl, overr, map, seq):
    seq = seq.upper()
    if (overl in seq and overr in seq):
        for ss in map:
            newseq = overl + ss + overr
            if (newseq in seq):
                return map.index(ss)
    seq = complementary(seq)
    if (overl in seq and overr in seq):
        for ss in map:
            newseq = overl + ss + overr
            if (newseq in seq):
                return map.index(ss)
    return -1


def get_value_from_unit(unit_list):
    group_no = unit_list[2] + unit_list[1] * 32 + ((unit_list[0] & 0x3c) // 4) * 32 * 32
    ingroup_no = unit_list[3]
    elem1 = unit_list[4] + (unit_list[5] & 0x03) * 64
    elem2 = (unit_list[5] & 0x3c) // 4 + (unit_list[6] & 0x0f) * 16
    elem3 = (unit_list[6] & 0x30) // 16 + unit_list[7] * 4
    elem4 = unit_list[8] + (unit_list[0] & 0x03) * 64
    print(group_no, ingroup_no, elem1, elem2, elem3, elem4)
    return group_no, ingroup_no, elem1, elem2, elem3, elem4


def row_translate_2216_10(fasq, outpath, log_fd, fmap, file_size, k_r=9, m_r=1, k=20, m=10, size=4, size_r=1):
    log_fd.write("---------------------||This---Is---Without---No---Row---Translation---Start||---------------------\n")
    list_t = []
    map1, map2, map3, map4, map5, map6, map7, map8, map9, map10 = getallmap(fmap)

    group_size = math.ceil(file_size / (k * size))
    index_size = group_size * (k + m)
    for i in range(index_size):
        list_t.append([])
        list_t[i].append([])
        list_t[i].append([])
        list_t[i].append([])
        list_t[i].append([])
    for val in fasq:
        unit_list = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
        unit_list[0] = unit_from_overlap('AAGCTT', 'CTCC', map1, val)
        unit_list[1] = unit_from_overlap('CTCC', 'ACTA', map2, val)
        unit_list[2] = unit_from_overlap('ACTA', 'GGTA', map3, val)
        unit_list[3] = unit_from_overlap('GGTA', 'ATAG', map4, val)
        unit_list[4] = unit_from_overlap('ATAG', 'AAGT', map5, val)
        unit_list[5] = unit_from_overlap('AAGT', 'CTGA', map6, val)
        unit_list[6] = unit_from_overlap('CTGA', 'CAAC', map7, val)
        unit_list[7] = unit_from_overlap('CAAC', 'AGGA', map8, val)
        unit_list[8] = unit_from_overlap('AGGA', 'TGTG', map9, val)
        unit_list[9] = unit_from_overlap('TGTG', 'GAATTC', map10, val)

        if (unit_list.count(-1) == 0):
            ##exchange
            check = ecode.dna_check_2216_10(en_buf=unit_list, k=9, m=1, len=1)  ###重组问题
            if (check == True):

                print("change: ", unit_list[0], unit_list[1], unit_list[2], unit_list[3], unit_list[4], unit_list[5],
                      unit_list[6], unit_list[7], unit_list[8], unit_list[9])
                elem1 = (unit_list[4] & 0x1f) + (unit_list[0] & 0x20)
                elem2 = (unit_list[5] & 0x1f) + (unit_list[0] & 0x10) * 2
                elem3 = (unit_list[6] & 0x1f) + (unit_list[0] & 0x08) * 4
                elem4 = (unit_list[7] & 0x1f) + (unit_list[1] & 0x10) * 2
                elem5 = (unit_list[8] & 0x1f) + (unit_list[1] & 0x08) * 4

                index1 = (unit_list[4] & 0x20) + (unit_list[5] & 0x20) // 2 + (unit_list[6] & 0x20) // 4 + (
                        unit_list[0] & 0x07)
                index2 = (unit_list[7] & 0x20) // 2 + (unit_list[8] & 0x20) // 4 + (unit_list[1] & 0x07)
                unit_list[4] = elem1
                unit_list[5] = elem2
                unit_list[6] = elem3
                unit_list[7] = elem4
                unit_list[8] = elem5
                unit_list[0] = index1
                unit_list[1] = index2
                print("ori: ", unit_list[0], unit_list[1], unit_list[2], unit_list[3], unit_list[4], unit_list[5],
                      unit_list[6], unit_list[7], unit_list[8], unit_list[9])
                group_no, ingroup_no, elem1, elem2, elem3, elem4 = get_value_from_unit(unit_list)

                if (group_size > group_no and ingroup_no < k + m and unit_list[0] // 4 <= 16 and unit_list[1] <= 32 and
                        unit_list[2] <= 32 and unit_list[3] <= 32):
                    list_t[group_no * (k + m) + ingroup_no][0].append(elem1)
                    list_t[group_no * (k + m) + ingroup_no][1].append(elem2)
                    list_t[group_no * (k + m) + ingroup_no][2].append(elem3)
                    list_t[group_no * (k + m) + ingroup_no][3].append(elem4)
                else:
                    log_fd.write("error(pass check) :outof index range\n")
                    log_fd.write(
                        "group_no:" + str(group_no) + '\t' + 'index1:' + str(unit_list[1]) + '\t' + 'index2:' + str(
                            unit_list[5]) + '\t' + 'index3:' + str(unit_list[3]) + '\n')

        elif (unit_list.count(-1) == 1):
            i = unit_list.index(-1)
            if (i == k_r):
                print("change: ", unit_list[0], unit_list[1], unit_list[2], unit_list[3], unit_list[4], unit_list[5],
                      unit_list[6], unit_list[7], unit_list[8], unit_list[9])
                elem1 = (unit_list[4] & 0x1f) + (unit_list[0] & 0x20)
                elem2 = (unit_list[5] & 0x1f) + (unit_list[0] & 0x10) * 2
                elem3 = (unit_list[6] & 0x1f) + (unit_list[0] & 0x08) * 4
                elem4 = (unit_list[7] & 0x1f) + (unit_list[1] & 0x10) * 2
                elem5 = (unit_list[8] & 0x1f) + (unit_list[1] & 0x08) * 4

                index1 = (unit_list[4] & 0x20) + (unit_list[5] & 0x20) // 2 + (unit_list[6] & 0x20) // 4 + (
                        unit_list[0] & 0x07)
                index2 = (unit_list[7] & 0x20) // 2 + (unit_list[8] & 0x20) // 4 + (unit_list[1] & 0x07)
                unit_list[4] = elem1
                unit_list[5] = elem2
                unit_list[6] = elem3
                unit_list[7] = elem4
                unit_list[8] = elem5
                unit_list[0] = index1
                unit_list[1] = index2
                print("ori: ", unit_list[0], unit_list[1], unit_list[2], unit_list[3], unit_list[4], unit_list[5],
                      unit_list[6], unit_list[7], unit_list[8], unit_list[9])
                group_no, ingroup_no, elem1, elem2, elem3, elem4 = get_value_from_unit(unit_list)
                if (group_size > group_no and ingroup_no < k + m and unit_list[0] // 4 < 16 and unit_list[1] < 32 and
                        unit_list[2] < 32 and unit_list[3] < 32):
                    list_t[group_no * (k + m) + ingroup_no][0].append(elem1)
                    list_t[group_no * (k + m) + ingroup_no][1].append(elem2)
                    list_t[group_no * (k + m) + ingroup_no][2].append(elem3)
                    list_t[group_no * (k + m) + ingroup_no][3].append(elem4)

            else:
                extended_matrix = ecode.gf_gen_rs_matrix(9, 1)
                pnerrs = 1
                pnsrcerrs = 1
                err_list = [i]
                in_err = np.zeros(10, dtype=int)
                in_err[i] = 1

                dest = ecode.dna_decode(
                    extended_matrix,
                    en_buf=np.array(unit_list).reshape((10, 1)), src_in_err=in_err, src_err_list=err_list,
                    pnerrs=pnerrs, pnsrcerrs=pnsrcerrs, k=9, m=1, len=1)

                unit_list[i] = dest[0].tolist()[0]
                print("change: ", unit_list[0], unit_list[1], unit_list[2], unit_list[3], unit_list[4], unit_list[5],
                      unit_list[6], unit_list[7], unit_list[8], unit_list[9])
                elem1 = (unit_list[4] & 0x1f) + (unit_list[0] & 0x20)
                elem2 = (unit_list[5] & 0x1f) + (unit_list[0] & 0x10) * 2
                elem3 = (unit_list[6] & 0x1f) + (unit_list[0] & 0x08) * 4
                elem4 = (unit_list[7] & 0x1f) + (unit_list[1] & 0x10) * 2
                elem5 = (unit_list[8] & 0x1f) + (unit_list[1] & 0x08) * 4

                index1 = (unit_list[4] & 0x20) + (unit_list[5] & 0x20) // 2 + (unit_list[6] & 0x20) // 4 + (
                        unit_list[0] & 0x07)
                index2 = (unit_list[7] & 0x20) // 2 + (unit_list[8] & 0x20) // 4 + (unit_list[1] & 0x07)
                unit_list[4] = elem1
                unit_list[5] = elem2
                unit_list[6] = elem3
                unit_list[7] = elem4
                unit_list[8] = elem5
                unit_list[0] = index1
                unit_list[1] = index2
                print("ori: ", unit_list[0], unit_list[1], unit_list[2], unit_list[3], unit_list[4], unit_list[5],
                      unit_list[6], unit_list[7], unit_list[8], unit_list[9])

                group_no, ingroup_no, elem1, elem2, elem3, elem4 = get_value_from_unit(unit_list)
                if (group_size > group_no and ingroup_no < k + m and unit_list[0] // 4 < 16 and unit_list[1] < 32 and
                        unit_list[2] < 32 and unit_list[3] < 32):
                    list_t[group_no * (k + m) + ingroup_no][0].append(elem1)
                    list_t[group_no * (k + m) + ingroup_no][1].append(elem2)
                    list_t[group_no * (k + m) + ingroup_no][2].append(elem3)
                    list_t[group_no * (k + m) + ingroup_no][3].append(elem4)

    list_row = []
    fd = open(outpath, "w+")
    for i in range(index_size):
        list_row.append([])
        if (len(list_t[i][0]) != 0):
            max_count = max((list_t[i][0]), key=list_t[i][0].count)
            list_row[i].append(max_count)
        else:
            list_row[i].append(-1)
        if (len(list_t[i][1]) != 0):
            max_count = max((list_t[i][1]), key=list_t[i][1].count)
            list_row[i].append(max_count)
        else:
            list_row[i].append(-1)
        if (len(list_t[i][2]) != 0):
            max_count = max((list_t[i][2]), key=list_t[i][2].count)
            list_row[i].append(max_count)
        else:
            list_row[i].append(-1)
        if (len(list_t[i][3]) != 0):
            max_count = max((list_t[i][3]), key=list_t[i][3].count)
            list_row[i].append(max_count)
        else:
            list_row[i].append(-1)
    fd.close()
    return list_t, list_row


def column_translate_2216_10(list_row, outpath, file_size, log_fd, k=20, m=10, size=4):
    fd = open(outpath, "wb+")
    log_fd.write("---------------------||This---Is---Column---Translation---Start||---------------------\n")
    list_col = []
    count = 0
    err_num = 0
    pass_check = 0
    # file_size = 0
    group_no = math.ceil(file_size / (k * size))
    extended_matrix = ecode.gf_gen_rs_matrix(20, 10)
    for i in range(group_no):
        list_buf = list_row[i].tolist()
        count_1 = np.zeros(30, dtype=int).tolist()
        for j in range(k + m):  
            if (list_buf[j].count(-1) > 0):
                count_1[j] = -1

        ocount = count_1.count(-1)
        if (ocount == 0):
            check = ecode.dna_check_2216_10(en_buf=list_buf, k=20, m=10, len=4)
            if (check):
                pass_check += 1
        elif (ocount <= 10):
            pnerrs = 0
            pnsrcerrs = 0
            err_list = []
            in_err = np.zeros(30, dtype=int)
            ind = 0
            while (ind < len(count_1)):
                if (count_1[ind] == -1):
                    pnerrs = pnerrs + 1
                    if (ind < k):
                        pnsrcerrs = pnsrcerrs + 1
                    in_err[ind] = 1
                    err_list.append(ind)
                ind = ind + 1

            dest = ecode.dna_decode(
                extended_matrix,
                en_buf=np.array(list_buf).reshape((30, 4)), src_in_err=in_err, src_err_list=err_list,
                pnerrs=pnerrs, pnsrcerrs=pnsrcerrs, k=20, m=10, len=4)
            di = 0
            while (di < len(dest)):
                list_buf[err_list[di]] = dest[di].tolist()
                di = di + 1

        list_col.append(list_buf[0:k])
        for j in range(k):
            for i in range(size):
                if (count < file_size):
                    if (list_buf[j][i] < 0):
                        fd.write(int(0).to_bytes(length=1, byteorder='big'))
                        count += 1
                        err_num += 1
                    else:
                        fd.write(int(list_buf[j][i]).to_bytes(length=1, byteorder='big'))
                        count += 1

    log_fd.write("---------------------||This---Is---Column---Translation---End||---------------------\n")
    fd.close()
    return list_col, err_num


def get_input_dic_fromfq(in_path):
    fasta = []
    i = 0
    with open(in_path) as file:
        sequence = ""
        for line in file:
            if i % 1 == 0:
                fasta.append(line.rstrip('\n'))
            i += 1
            
    return fasta


def de_code_2216_10(in_path, output, log_path, fmap, file_size, k=6, m=2, k_r=9, m_r=1, len_c=4, len_r=1):
    try:
        fd = open(in_path, 'rb+', )
    except:
        print("input_error")
        return 1, "input_error"

    out_path = in_path[:-4] + '_row_translate.txt'

    fdw = open(output, 'w')
    fd_log = open(log_path, 'w')

    fastq = get_input_dic_fromfq(in_path)

    list_row_t, list_row = row_translate_2216_10(fastq, out_path, fmap=fmap, log_fd=fd_log, file_size=file_size)
    group_no = math.ceil(file_size / (k * len_c))
    index_size = group_no * (k + m)
    list_row = np.array(list_row)
    list_row = list_row.reshape((group_no, k + m, len_c))

    list_col, err_num = column_translate_2216_10(list_row, output, log_fd=fd_log, file_size=file_size)

    count = 0
    fd_log.write("Tranlate result:\n")
    for i in range(group_no):
        if (count < file_size):
            fd_log.write(str(i) + ':\t' + str(list_col[i]) + '\n')
            count += k
        else:
            count += k

    fd.close()
    fdw.close()
    fd_log.close()
    err_info = []

    err_info.append("translated stripes：" + str(index_size) + "\t have payloads：" + str(group_no * 4) + "\n")
    err_info.append("dna reads：" + str(len(fastq)) + "\n")

    return 0, err_info


def file_compare(origin_path, compare_path):
    o_fd = open(origin_path, 'rb')
    c_fd = open(compare_path, 'rb')

    buffer1 = o_fd.read()
    buffer2 = c_fd.read()

    mismatch = 0

    lenth = len(buffer1)
    for i in range(lenth):
        if (buffer1[i] != buffer2[i]):
            mismatch += 1
            print("mismatch", i, buffer1[i], buffer2[i])
    o_fd.close()
    c_fd.close()
    return mismatch


if __name__ == '__main__':
    # input = sys.argv[1] 11128,434,14490,18605
    os.chdir("//Users//freya//PycharmProjects//DNAStorage20211108//encode_decode//data/2216-10")
    input = "horse_seq_v3.txt"#sequencing data
    output = input[:-3] + "_translate.mp4"
    log_path = output[:-4] + "_log.txt"
    mark, err_info = de_code_2216_10(input,
                                     output,
                                     log_path,
                                     fmap="//Users//freya//PycharmProjects//DNAStorage20211108//encode_decode//map"
                                          "//2216-10-mw.txt",
                                     file_size=18605, k=20,
                                     m=10, k_r=9, m_r=1, len_c=4, len_r=1)
    print(mark, err_info)
    mismatch = file_compare("horse.mp4", "horse_seq_v3._translate.mp4")
    print(mismatch)

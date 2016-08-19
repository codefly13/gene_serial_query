#!/usr/bin/python
# -*- coding:utf-8 -*- 

import pandas as pd

def get_chromsome_fragment(data, input_filename, validate_filename):
    reverse_it = {'A':'T','C':'G', "T":"A", "G":"C", 'a':'t','c':'g', "t":"a", "g":"c"}

    for each_row in data.get_values():
        each_row = list(each_row)
        chromosome_id = each_row[0].strip()
        geneID = each_row[9].split("=")[-1]
        start = int( each_row[2] )
        end = int( each_row[3] )
        direction = each_row[4]

#         print(each_row)
        
        try :
#             with open("ysxhg19.fa") as f:
            with open(input_filename) as f:
                found = False
#                for each_line in f:
                while (True):
                    each_line = f.readline()
                    if each_line == "": break
                    each_line = each_line.strip()
#                     print('start',len(each_line),">=5",each_line[:4],"=" , ">chr",each_line[4:],"=" , chromosome_id)
                    if 'chr' in each_line and each_line.split("chr")[-1].strip() == chromosome_id:
                        found = True
                        first_line = f.readline()
                        first_line = first_line.strip()
                        line_len = len(first_line)

                        start_line = start // line_len + ((start%line_len) > 0) # 假设起始行是1行
                        end_line = end // line_len + ((end%line_len) > 0)

                        start_index = start%line_len
                        if start_index == 0:
                            start_index = line_len
                        start_index -= 1 # 转换为以0开始的索引

                        number = end - start + 1

                        fragment = ""
                        line_index = 1
                        if line_index >= start_line:
#                             print("index",line_index)
                            fragment += first_line
                        
                        for each_new_line in f:
                            if 'chr' in each_line and each_line.split("chr")[-1].strip() != chromosome_id:
#                                print("the range fragment is out of range")
                                print("起始位或终止位超出范围了".decode('utf-8').encode('gbk'))
                    
                            line_index += 1
                            if line_index < start_line:
                                continue
                            if line_index > end_line:
                                break
                            
#                             print("index",line_index)
                            fragment += each_new_line.strip() 
                            
#                         print('origin')
#                         print(fragment)
#                         print('uuuuuuuuuuuuuuuuu')
                        fragment = fragment[start_index:start_index+number]
                        fragment2 = fragment

                        """反向互补配对"""
                        if direction == '-':
                            new_fragment = ""
                            for b in fragment:
                                b = reverse_it.get(b,b)
                                new_fragment = b + new_fragment
                            fragment = new_fragment

                        print(fragment)

                        if fragment != "":
                            with open("result.tmp", "a") as result_file:
#                                 print(gene_id, fragment, file=result_file)
                                result_file.write(" ".join([str(each) for each in [gene_id, fragment]]) + "\n")
                            
                            with open("result2.tmp", "a") as result_file2:
#                                 print(start, end, chromosome_id, direction, fragment, file=result_file2)
                                result_file2.write(" ".join([str(each) for each in [start, end, chromosome_id, direction, fragment]]) + "\n")
                        
#                         with open("validate.tmp", "a") as validate_file:
                        with open(validate_filename, "a") as validate_file:
                            if fragment == "":
#                                 print("may occur error", file=validate_file)
                                validate_file.write("may occur error\n")
#                                 print(start,end,start_line,end_line,start_index, number, fragment2, file=validate_file)
                                validate_file.write(" ".join([str(each)for each in [start,end,start_line,end_line,start_index, number, fragment2]]) + "\n")
#                                 print("---", file=validate_file)
                                validate_file.write("---\n")
                            else:
#                                 print(*each_row, end=" ", file = validate_file)
                                validate_file.write(" ".join([str(each) for each in each_row]) + " ")
#                                 print(fragment, file = validate_file)
                                validate_file.write(str(fragment) + "\n")
#                                 print(start,end,start_line,end_line,start_index, number,fragment2, fragment, file=validate_file)
                                validate_file.write(" ".join([str(each) for each in [start,end,start_line,end_line,start_index, number,fragment2, fragment]]) + "\n")
                        break
                if found == False:
                    print("chromosome" , chromosome_id, "not found")
        except EOFError:
            print("eroor")
            pass
        

print("查询染色体片段".decode('utf-8').encode('gbk'))
gene_input_filename = raw_input("输入基因文件路径，如ysxgff.gff\n".decode('utf-8').encode('gbk'))
chromsome_input_filename = raw_input("输入染色体文件路径，如ysxhg19.fa\n".decode('utf-8').encode('gbk'))
validate_filename = raw_input("输入临时文件路径\n".decode('utf-8').encode('gbk'))
print("输入gene_id（q退出程序）".decode('utf-8').encode('gbk'))

while(True):
# \s*正则识别多个空格
#     gd = pd.read_csv('ysxgff.gff',sep="\s*", header = None, chunksize=10000, encoding = 'gbk', engine='python')
    gd = pd.read_csv(gene_input_filename,sep="\s*", header = None, chunksize=10000, encoding = 'gbk', engine='python')
    
    gene_id = str(raw_input())
    
    if gene_id == str("q") or gene_id == str("Q"):
        print("quit")
        break
    
    cnt = 0
    while(True):
        try:
            
            gdc = gd.get_chunk()
            gdc = gdc.astype(str)
            gdc2 = gdc.copy()
            gdc2.ix[:, 6] = gdc2.ix[:, 6].apply(lambda x: x.split("=")[-1])
            gdc2.ix[:, 7] = gdc2.ix[:, 7].apply(lambda x: x.split("=")[-1])
            gdc2.ix[:, 8] = gdc2.ix[:, 8].apply(lambda x: x.split("=")[-1])
            gdc2.ix[:, 9] = gdc2.ix[:, 9].apply(lambda x: x.split("=")[-1])
            
            
#             print( gdc2.apply(lambda x: x[9] == gene_id and x[1] == "CDS" , axis=1) )
            
            gdc = gdc[ gdc2.apply(lambda x: x[9] == gene_id and x[1] == "CDS" , axis=1) ]
#            if len(gdc) > 0:
#                print("len",len(gdc))
            
#             print("bool=>\n", gdc2[9]==gene_id ,'\n---\n', gdc2[1]=="CDS", "end")
            
#             print("select=>\n", gdc[gdc2[9]==gene_id])
            
#             print(cnt,"=>\n",gdc2)
            
#             gdc = gdc[gdc2[9]==gene_id]
            
            cnt += len(gdc)
                
#             print(cnt,"=>\n",gdc)
            
            if len(gdc) > 0:
                get_chromsome_fragment(gdc, chromsome_input_filename, validate_filename)
            
        except StopIteration:
            if cnt == 0:
                print("gene_id notFound")
                with open("result.tmp", "a") as result_file:
#                     print(gene_id, "notFound", file=result_file)
                    result_file.write(str(gene_id) + " " + "notFound\n")
            break
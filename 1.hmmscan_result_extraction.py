
import re
import os
import linecache
input_dir = f'C:\\Users\\P\\Desktop\\hmmscan匹配问题解决\\bacteria'
out_dir=f'C:\\Users\\P\\Desktop\\hmmscan匹配问题解决\\out\\'
filenames = os.listdir(input_dir)


for x in filenames:#############一个out文件
	file = open(f"{input_dir}\\{x}","r")###输入文件
	result_list=open(f'{out_dir}' + f'{x[:-4]}' + '.txt','w')##输出文件
	result_list.write('Protein_ID\tDomain_ID\tDomain_Name\tDomain_Description\tProtein_Evalue\tDomain_Evalue\tQuery_First\tQuery_End\tDomain_Sequence\tNormalized_Domain_Sequence\tHMM_From\tHMM_To\tHMM_File\n')
		

	list_num=[]
	list_num2=[]
	for num,value in enumerate(file):
	    if 'Query:' == value[:6]:
		        list_num.append(num)
	    elif '//' == value[:2]:
		        list_num2.append(num)
	file.close()
	for t in zip(list_num,list_num2):##########一种蛋白
	    
	    file_del=linecache.getlines(f'{input_dir}\\{x}')[t[0]:t[1]]###输入文件
	    name=file_del[0]
	    pro_id=name.split()[1]##########蛋白id
	    domain_d_nam={}
	    domain_d_dis={}
	    protein_edir={}
	    domain_sidr={}
	    domain_ddir={}
	    for line in file_del[5:]:
	        line_list=line.split()
	        if ('[No hits detected that satisfy reporting thresholds]' not in line ) and (len(line_list) !=0):
	            
	            
	            try :
	                e_value=line_list[0]
	                test_g=float(e_value)
	                domain_id=line_list[8]
	                protein_edir[domain_id]=e_value#################蛋白内对应domain对应evluea值
	            except ValueError:
	                pass

	        elif '[No hits detected that satisfy reporting thresholds]' in line:
	            break
	        elif len(line_list)==0:
	            
	            break                     ###########获取所有蛋白
	    pro_del=[]
	    for i in range(len(file_del)):
	        if '>>' ==file_del[i][:2]:
	            pro_del.append(i)
	        elif 'Internal pipeline statistics summary:' in file_del[i]:
	            pro_del.append(i)
	    
	    for i in range(len(pro_del)-1):#############一种蛋白一个结构域内循环
	        f_i=pro_del[i]
	        a_i=pro_del[i+1]
	        ltaa=[]
	        dorange=[]
	        for s_i in range(len(file_del[f_i:a_i])):
	            line=file_del[f_i:a_i][s_i]
	            
	            if '>>'==line[:2]:
	                line_list=line.split(' ',3)
	                domain_id=line_list[1]
	                domain_na_di=line_list[3]
	                domain_na_di_split=domain_na_di.split(', ')
	                domain_name=domain_na_di_split[0]
	                if len(domain_na_di_split)==2:
	                	            domain_dis=domain_na_di_split[1][:-2]
	                else:
	                	domain_dis=''
	                	domain_name=domain_na_di_split[0][:-2]
	                ltaa.append(domain_id)###########domain信息
	                
	                domain_d_nam[domain_id]=domain_name
	                domain_d_dis[domain_id]=domain_dis
	            else:
	                line_list=line.split()
	                try:
                        
	                   float(line_list[3]) 
	                   
	                       
	                   if (len(line_list) !=4 )and ( 'i' not in line_list[3] ) and ( 'I' not in line_list[3] ) and('n' not in line_list[3] ):
	                       domaid=ltaa[0]
	                       domain_ddir.setdefault(domaid,[]).append(line_list)
	                       
	                except IndexError:
	                   pass
	                    
	                except ValueError:
	                    
	                    if line_list[0]=='==':
	                        dorange.append(s_i)
	                        
	                    else:
	                        pass
	        dorange.append(len(file_del[f_i:a_i])-2)          
	                       
	        
	        file_del_del=file_del[f_i:a_i]    
	        for d_i in range(len(dorange)-1):#####一种蛋白一个结构域一个结构域
	            fd_i=dorange[d_i]
	            ad_i=dorange[d_i+1]
	            
	            seq=''
	            for dd_i in range(len(file_del_del[fd_i:ad_i])):
	                line=file_del_del[fd_i:ad_i][dd_i]
	                line_list=line.split()
	                #proid=ltaa[0]
	                if len(line_list)!=0:
	                    if line_list[0]==pro_id:
	                        seq=f'{seq}{line_list[2]}'
	            seq=re.sub(r'-','',seq)
	            proid=ltaa[0]
	            domain_sidr.setdefault(proid,[]).append(seq)
	                
	            
	    for key in protein_edir:
	        try:
	            result_edir=protein_edir[key]
	            result_ddir=domain_ddir[key]
	            result_sdir=domain_sidr[key]
	            
        	    result_name=domain_d_nam[key]
        	    result_dis=domain_d_dis[key]
	            for i in range(len(result_ddir)):
	                    rr_ddir=result_ddir[i]
	                    
	                    ss_sdir=result_sdir[i]
	                            
	                                
	                    result_list.write(f'{pro_id}\t{key}\t{result_name}\t{result_dis}\t{result_edir}\t{rr_ddir[4]}\t{rr_ddir[9]}\t{rr_ddir[10]}\t{ss_sdir}\t{ss_sdir}\t{rr_ddir[6]}\t{rr_ddir[7]}\t{ss_sdir}\tpfam-A\n')            
	                    #result_list.write(f'{key}\t{domain_list[1]}\t{domain_list[0]}\t{domain_list[2]}\t{result_edir}\t{rr_ddir[4]}\t{rr_ddir[9]}\t{rr_ddir[10]}\t{ss_sdir}\t{ss_sdir}\t{rr_ddir[6]}\t{rr_ddir[7]}\tpfam-A\n')
	        except KeyError:
	            pass                             
	result_list.close()
            

        
        
	
        
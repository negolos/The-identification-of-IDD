import pandas as pd
import re
import os


spotd_dir = f'D:\\Plant proteostasis 2nd edition\\program\\test_for_two'
doamin_info=pd.read_table('D:\\Plant proteostasis 2nd edition\\program\\domain_info.txt')
domain_result=open('D:\\Plant proteostasis 2nd edition\\program\\domain_result.txt','w')


filenames = os.listdir(spotd_dir)
def get_domain_site(g):
    star=g['Query_First']
    end=g['Query_end']
    psd=g['Domain_ID']
    dom_si=(star,end,psd)
    return dom_si
doamin_info['DOST']=doamin_info.apply(get_domain_site,axis=1)
domain_info_GR=doamin_info.groupby(['Protein_ID'])
def get_whole(y):
    q=[]
    for x in y:
        if x not in q:
        
            q.append(x)
    return q
dost_G=domain_info_GR['DOST'].agg(get_whole)
dost_G=dost_G.to_frame()
dost_G=dost_G.reset_index()

domain_result.write('Protein_ID\tDomain_ID\tStart\tEnd\tDSDR\tCDRN\n')
for x in filenames:
    file=pd.read_table(f'{spotd_dir}\\{x}',skiprows=2, header=None)
    file.columns =['Pos','AA','P(D)','Lab']
    proid=x[:-7]

    domain_selected=dost_G[dost_G['Protein_ID']==proid]

    for index,row in domain_selected.iterrows():
        site=row['DOST']
        for sit in site:
            Start_site=sit[0]
            End_site=sit[1]
            Domain_ID=sit[2]
            spotd_cut=file[Start_site-1:End_site]
            spotd_cut_num=spotd_cut.shape[0]
            spotd_list=spotd_cut['Lab'].tolist()
            spotd_str="".join(spotd_list)
            spotd_cut_D_num =spotd_cut[spotd_cut['Lab']=='D'].shape[0]
            DSDR=spotd_cut_D_num/spotd_cut_num*100
                

            
            g=re.findall(r'D+',spotd_str)
            CDRN=0
            if len(spotd_str) >=50:
                for x in g:
                    if len(x) >=20:
                        CDRN +=1
            else:
                for x in g:
                    if (len(x)/len(spotd_str))>=0.4:
                        CDRN +=1
            if End_site - Start_site >=20:
                domain_result.write(f'{proid}\t{Domain_ID}\t{Start_site}\t{End_site}\t{DSDR}\t{CDRN}\n')
domain_result.close()

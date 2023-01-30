import pandas as pd
import networkx as nx
import numpy as np
import itertools
import openpyxl as xl

try:
    #----------Read Input Data----------------
    feeder=pd.read_excel('GridData.xlsx', sheet_name='Lines')
    Generation=pd.read_excel('GridData.xlsx', sheet_name='Generation')
    Slack=pd.read_excel('GridData.xlsx', sheet_name='SlackBus')

    #--------------------------------- Caso con distintas contingencias -------------------------

    HM=feeder[feeder['HM']!=0].index.tolist() #Lines for maintenance 
    lines_under_maintenance=sum(feeder['HM']!=0) #Number of lines under maintenance
    bina=list(itertools.product(*[[1,0]]*lines_under_maintenance)) 
    tab_bin = pd.DataFrame(bina,columns=list(np.where(feeder['HM']!=0)[0])) #Posibles combinaciones
    pt_dic={}
    ld_dic={}

    data_pt={'Topology':[],'Nodes':[],'Lines':[],'PTDF':[]}
    table_pt=pd.DataFrame(data_pt)
    data_df={'Topology':[],'Line_l':[],'Line_m':[],'LODF':[]}
    table_df=pd.DataFrame(data_df)
    len_datapt=[0]
    len_datadf=[0]

    nodes=len(set.union(set(feeder['Node_i']),set(feeder['Node_j'])))
    lines_total=len(feeder)
    
    availability=[1]*len(tab_bin)
    for k in range(len(tab_bin)):
        li_ava=[]
        combi=list(tab_bin.iloc[k])
        for m in range(len(combi)):
            if combi[m]==0:
                li_ava+=[HM[m]]
            else:
                continue
        new_feeder=feeder.drop(li_ava)  
        l=len(new_feeder)
        n=len(set.union(set(new_feeder['Node_i']),set(new_feeder['Node_j'])))

        Bp=np.diag([1/(1j*feeder['X'][mq]) for mq in range(l)])#Primitive Susceptances Matrix
        A=np.zeros((n,l))
        cn=0
        cl=0
        for ldw in range(l):
            From=feeder['Node_i'][ldw]
            To=feeder['Node_j'][ldw]
            for m in range(n):
                if From==m+1:
                    A[m,ldw]=1
                elif To==m+1:
                    A[m,ldw]=-1
                    
        if 0 in np.sum(abs(A),axis=1) or len(A)!=nodes:
            print('the combination ', combi, 'leaves one node isolated from the system, therefore it will not be taken into account')
            availability[k]=0
            pt_dic[k]=np.zeros((len(A[0,:]),nodes))
            ld_dic[k]=np.zeros((len(A[0,:]),len(A[0,:])))
            len_datapt+=[len_datapt[-1]+len(new_feeder['Line'])+2]
            len_datadf+=[len_datadf[-1]+len(new_feeder['Line'])+2]
            continue
            
        #--Matrix Bbus--
        Bbus=A@Bp@A.T
        
        #--Building the X matrix selecting the slack node from the variable "Slack"--
        # Sets up to "0" the columns of the slack node
        slack_bus = Slack.iloc[0,0]
        Bbus[:,slack_bus]=0
        # Sets up to "0" the rows of the slack node
        Bbus[slack_bus,:]=0
        # Sets up to "1" the resistive value of the row and column of the slack node
        Bbus[slack_bus,slack_bus]=1
        
        X=np.linalg.inv(Bbus)
        X=abs(X)
        X[slack_bus,slack_bus]=0
        
        Al=A.T #Incidence Matrix Line-Node
        Bd=abs(Bp)

        #-------------Line-Node Sensibility Factors Calculation ---------------
        lodf=np.zeros((l,l))
        ptdf_rr=np.zeros((l,l))
        Inci=Al@np.identity(n)
        
        ptdf_rn=Bd@Al@X#PTDF Line-Node
        ptdf_rr=Bd@Al@X@Al.T#PTDF Line-Line
 
        #--LODFS CALCULATION--
        for s in range(l):
            for kl in range(l):
                if ptdf_rr[kl,kl]==1:
                    ptdf_rr[kl,kl]=0
                lodf[s,kl]=ptdf_rr[s,kl]*(1/(1-ptdf_rr[kl,kl]))
        
        np.fill_diagonal(lodf,0)    
        
        ptdf_new=np.zeros((lines_total,nodes))

        # for k5 in li_ava:
        #     ptdf_rn=np.insert(ptdf_rn, k5, np.zeros((1, nodes)), 0)
        #     lodf=np.insert(lodf, k5, np.zeros((1, len(new_feeder))), 0)
    
        # for k6 in li_ava:
        #     lodf=np.insert(lodf, k6, np.zeros((1, lines_total)), 1)
            
            
        pt_dic[k]=ptdf_rn
        ld_dic[k]=lodf
        pt_nf,pt_nc=ptdf_rn.shape
        df_nf,df_nc=lodf.shape
        
        fil_pf=[]
        row_pf=list(feeder['Line'])
        col_pf=list(range(1,nodes+1))*pt_nf
        for k1 in row_pf:
            fil_pf+=list(np.repeat(k1,pt_nc))
   
        fil_df=[]
        row_df=list(feeder['Line'])
        col_df=list(feeder['Line'])*df_nf
        for k2 in row_df:
            fil_df+=list(np.repeat(k2,df_nc))
          
    #     data_pt2={'Topology': [k]*len(fil_pf), feeder['Node_i']:col_pf,'Lines':fil_pf,'PTDF':list(ptdf_rn.flatten())}   
    #     data_df2={'Topology': [k]*len(fil_df),'Line_l':col_df,'Line_m':fil_df,'LODF':list(lodf.flatten())}
    #     table_pt2=pd.DataFrame(data_pt2)
    #     table_df2=pd.DataFrame(data_df2)
    #     table_pt=pd.concat([table_pt,table_pt2])
    #     table_df=pd.concat([table_df,table_df2])
    #     len_datapt+=[len_datapt[-1]+len(feeder['Line'])+2]
    #     len_datadf+=[len_datadf[-1]+len(feeder['Line'])+2]
        
     
    # tab_bin.insert(0, "Topology", list(range(len(tab_bin))), True)
    # tab_bin.insert(1,"Availability", availability, True)

    for k in pt_dic:
        df_ptdic=pd.DataFrame(pt_dic[k])
        df_lddic=pd.DataFrame(ld_dic[k]) 
    
    df_ptdic.set_index(feeder['Code_Line'],inplace=True, drop=True)
    gen_names = [Generation['Code_Gen']]
    df_ptdic.columns = gen_names
    df_ptdic.to_csv('PTDF.csv')



    #Calculation complete 
    print('The PTDF Matrix was calculated without errors')             

except:
    #Calculation with errors
    print('There are errors in the calculation')
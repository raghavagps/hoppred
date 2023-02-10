##############################################################################
#hoppred2.0 is developed for predicting allergenic and non allergenic        #
#protein from their primary sequence. It is developed by Prof G. P. S.       #
#Raghava's group. Please cite : hoppred 2.0                                  #
# ############################################################################
import argparse  
import warnings
import pickle
import os
import re
import sys
import numpy as np
import pandas as pd
import joblib

warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser(description='Please provide following arguments. Please make the suitable changes in the envfile provided in the folder.') 

## Read Arguments from command
parser.add_argument("-i", "--input", type=str, required=True, help="Input: protein or peptide sequence in FASTA format or single sequence per line in single letter code")
parser.add_argument("-o", "--output",type=str, help="Output: File for saving results by default outfile.csv")
parser.add_argument("-t","--threshold", type=float, help="Threshold: Value between 0 to 1 by default 0.3")
parser.add_argument("-m","--model",type=int, choices = [1, 2], help="Model: 1: Allergen, 2: Non-Allergen, by default 1")
parser.add_argument("-d","--display", type=int, choices = [1,2], help="Display: 1:Allergen, 2: All peptides, by default 1")
args = parser.parse_args()

std = list("ACDEFGHIKLMNPQRSTVWY")
PCP= pd.read_csv('Data/PhysicoChemical.csv', header=None)
AAindices = 'Data/aaind.txt'
AAIndex = pd.read_csv('Data/aaindex.csv',index_col='INDEX');
AAIndexNames = pd.read_csv('Data/AAIndexNames.csv',header=None);
dir_1 = Features


def aac_comp(file):
    filename, file_extension = os.path.splitext(file)
    f = open(dir_1+'/sam_allcomp.aac', 'w')
    sys.stdout = f
    df = pd.read_csv(file, header = None)
    zz = df.iloc[:,0]
    print("AAC_A,AAC_C,AAC_D,AAC_E,AAC_F,AAC_G,AAC_H,AAC_I,AAC_K,AAC_L,AAC_M,AAC_N,AAC_P,AAC_Q,AAC_R,AAC_S,AAC_T,AAC_V,AAC_W,AAC_Y,")
    for j in zz:
        for i in std:
            count = 0
            for k in j:
                temp1 = k
                if temp1 == i:
                    count += 1
                composition = (count/len(j))*100
            print("%.2f"%composition, end = ",")
        print("")
    f.truncate()

def dpc_comp(file):
    q = 1
    f = open(dir_1+'/sam_allcomp.dpc', 'w')
    sys.stdout = f
    df = pd.read_csv(file, header = None)
    df1 = pd.DataFrame(df[0].str.upper())
    zz = df1.iloc[:,0]
    for s in std:
        for u in std:
            print("DPC"+str(q)+"_"+s+u, end=',')
    print("")
    for i in range(0,len(zz)):
        for j in std:
            for k in std:
                count = 0
                temp = j+k
                for m3 in range(0,len(zz[i])-q):
                    b = zz[i][m3:m3+q+1:q]
                   # b.upper()

                    if b == temp:
                        count += 1
                    composition = (count/(len(zz[i])-(q)))*100
                print("%.2f" %composition, end = ',')
        print("")
    f.truncate()

#############################TPC_COMP##############################
def tpc_comp(file):
    filename, file_extension = os.path.splitext(file)
    f = open(dir_1+'/sam_allcomp.tpc', 'w')
    sys.stdout = f
    df = pd.read_csv(file, header = None)
    zz = df.iloc[:,0]
    for s in std:
       for u in std:
           for m in std:
               print('TPC_'+s+u+m, end=",")
    print("")
    for i in tqdm(range(0,len(zz))):
        for j in std:
            for k in std:
                for m1 in std:
                    count = 0
                    temp = j+k+m1
                    for m3 in range(0,len(zz[i])):
                        b = zz[i][m3:m3+3]
                        if b == temp:
                            count += 1
                        composition = (count/(len(zz[i])-2))*100
                    print("%.2f" %composition, end = ',')
        print("")
    f.truncate()

############################PhysicoChemical Properties###################################
PCP= pd.read_csv('Data/PhysicoChemical.csv', header=None)

headers = ['PCP_PC','PCP_NC','PCP_NE','PCP_PO','PCP_NP','PCP_AL','PCP_CY','PCP_AR','PCP_AC','PCP_BS','PCP_NE_pH','PCP_HB','PCP_HL','PCP_NT','PCP_HX','PCP_SC','PCP_SS_HE','PCP_SS_ST','PCP_SS_CO','PCP_SA_BU','PCP_SA_EX','PCP_SA_IN','PCP_TN','PCP_SM','PCP_LR','PCP_Z1','PCP_Z2','PCP_Z3','PCP_Z4','PCP_Z5'];

def encode(peptide):
    l=len(peptide);
    encoded=np.zeros(l);
    for i in range(l):
        if(peptide[i]=='A'):
            encoded[i] = 0;
        elif(peptide[i]=='C'):
            encoded[i] = 1;
        elif(peptide[i]=='D'):
            encoded[i] = 2;
        elif(peptide[i]=='E'):
            encoded[i] = 3;
        elif(peptide[i]=='F'):
            encoded[i] = 4;
        elif(peptide[i]=='G'):
            encoded[i] = 5;
        elif(peptide[i]=='H'):
            encoded[i] = 6;
        elif(peptide[i]=='I'):
            encoded[i] = 7;
        elif(peptide[i]=='K'):
            encoded[i] = 8;
        elif(peptide[i]=='L'):
            encoded[i] = 9;
        elif(peptide[i]=='M'):
            encoded[i] = 10;
        elif(peptide[i]=='N'):
            encoded[i] = 11;
        elif(peptide[i]=='P'):
            encoded[i] = 12;
        elif(peptide[i]=='Q'):
            encoded[i] = 13;
        elif(peptide[i]=='R'):
            encoded[i] = 14;
        elif(peptide[i]=='S'):
            encoded[i] = 15;
        elif(peptide[i]=='T'):
            encoded[i] = 16;
        elif(peptide[i]=='V'):
            encoded[i] = 17;
        elif(peptide[i]=='W'):
            encoded[i] = 18;
        elif(peptide[i]=='Y'):
            encoded[i] = 19;
        else:
            print('Wrong residue!');
    return encoded;
def lookup(peptide,featureNum):
    l=len(peptide);
    peptide = list(peptide);
    out=np.zeros(l);
    peptide_num = encode(peptide);

    for i in range(l):
        out[i] = PCP[peptide_num[i]][featureNum];
    return sum(out);
def pcp_1(file):

    if(type(file) == str):
        seq = pd.read_csv(file,header=None);
        #seq=seq.T
        seq[0].values.tolist()
        seq=seq[0];
    else:
        seq  = file;

    l = len(seq);

    rows = PCP.shape[0]; # Number of features in our reference table
    col = 20 ; # Denotes the 20 amino acids

    seq=[seq[i].upper() for i in range(l)]
    sequenceFeature = [];
    sequenceFeature.append(headers); #To put property name in output csv

    for i in range(l): # Loop to iterate over each sequence
        nfeatures = rows;
        sequenceFeatureTemp = [];
        for j in range(nfeatures): #Loop to iterate over each feature
            featureVal = lookup(seq[i],j)
            if(len(seq[i])!=0):
                sequenceFeatureTemp.append(round(featureVal/len(seq[i]),3));
            else:
                sequenceFeatureTemp.append('NaN')

        sequenceFeature.append(sequenceFeatureTemp);

    out = pd.DataFrame(sequenceFeature);
    file = open(dir_1+'/sam_allcomp.pcp','w')
    with file:
        writer = csv.writer(file);
        writer.writerows(sequenceFeature);
    return sequenceFeature;


#################################DDOR####################################
def DDOR(file) :
    df = pd.read_csv(file, header = None)
    df1 = pd.DataFrame(df[0].str.upper())
    f = open(dir_1+'/sam_allcomp.ddr','w')
    sys.stdout = f
    for i in std:
        print('DDR_'+i, end=",")
    print("")
    for i in range(0,len(df1)):
        s = df1[0][i]
        p = s[::-1]
        for j in std:
            zz = ([pos for pos, char in enumerate(s) if char == j])
            pp = ([pos for pos, char in enumerate(p) if char == j])
            ss = []
            for i in range(0,(len(zz)-1)):
                ss.append(zz[i+1] - zz[i]-1)
            if zz == []:
                ss = []
            else:
                ss.insert(0,zz[0])
                ss.insert(len(ss),pp[0])
            cc1=  (sum([e for e in ss])+1)
            cc = sum([e*e for e in ss])
            zz2 = cc/cc1
            print("%.2f"%zz2,end=",")
        print("")
    f.truncate()

################################Shannon_Entropy residue#############################################
def SE_residue_level(filename):
    data=list((pd.read_csv(filename,sep=',',header=None)).iloc[:,0])
    data2=list((pd.read_csv(filename,sep=',',header=None)).iloc[:,0])
    Val=np.zeros(len(data))
    GH=[]
    for i in range(len(data)):
        my_list={'A':0,'C':0,'D':0,'E':0,'F':0,'G':0,'H':0,'I':0,'K':0,'L':0,'M':0,'N':0,'P':0,'Q':0,'R':0,'S':0,'T':0,'V':0,'W':0,'Y':0}
        data1=''
        data1=str(data[i])
        data1=data1.upper()
        allowed = set(('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'))
        is_data_invalid = set(data1).issubset(allowed)
        if is_data_invalid==False:
            print("Error: Please check for invalid inputs in the sequence.","\nError in: ","Sequence number=",i+1,",","Sequence = ",data[i],",","\nNOTE: Spaces, Special characters('[@_!#$%^&*()<>?/\|}{~:]') and Extra characters(BJOUXZ) should not be there.")
            return
        seq=data[i]
        seq=seq.upper()
        num, length = Counter(seq), len(seq)
        num=dict(sorted(num.items()))
        C=list(num.keys())
        F=list(num.values())
        for key, value in my_list.items():
             for j in range(len(C)):
                if key == C[j]:
                    my_list[key] = round(((F[j]/length)* math.log(F[j]/length, 2)),3)
        GH.append(list(my_list.values()))
    file= open(dir_1+'/sam_allcomp.ser','w', newline='')#output file
    with file:
        writer=csv.writer(file);
        writer.writerow(('SER_A','SER_C','SER_D','SER_E','SER_F','SER_G','SER_H','SER_I','SER_K','SER_L','SER_M','SER_N','SER_P','SER_Q','SER_R','SER_S','SER_T','SER_V','SER_W','SER_Y'));
        writer.writerows(GH);
    return(GH)

########################Shannon_Entropy Whole protein######################################
def entropy_single(seq):
    seq=seq.upper()
    num, length = Counter(seq), len(seq)
    return -sum( freq/length * math.log(freq/length, 2) for freq in num.values())

def SE(filename):
    data=list((pd.read_csv(filename,sep=',',header=None)).iloc[:,0])
#     print(data)
    Val=[]
    header=["SEP"]
    for i in range(len(data)):
        data1=''
        data1=str(data[i])
        data1=data1.upper()
        allowed = set(('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'))
        is_data_invalid = set(data1).issubset(allowed)
        if is_data_invalid==False:
            print("Error: Please check for invalid inputs in the sequence.","\nError in: ","Sequence number=",i+1,",","Sequence = ",data[i],",","\nNOTE: Spaces, Special characters('[@_!#$%^&*()<>?/\|}{~:]') and Extra characters(BJOUXZ) should not be there.")
            return
        Val.append(round((entropy_single(str(data[i]))),3))
        #print(Val[i])
        file= open(dir_1+'/sam_allcomp.sep1','w', newline='\n')#output file
        with file:
            writer=csv.writer(file,delimiter='\n');
            writer.writerow(header)
            writer.writerow(Val);
    return Val

##########################################CTC###################################
x = [1, 2, 3, 4, 5, 6,7]
p=[]
Y=[]
LS=[]


for i in range(len(x)):
    p=itertools.product(x,repeat=3)
    p=list(p)

def concatenate_list_data(list):
    result= ''
    for element in list:
        result += str(element)
    return result

for i in range(len(p)):
    LS.append(concatenate_list_data(p[i]))

def repstring(string):
    string=string.upper()
    char={"A":"1","G":"1","V":"1","I":"2","L":"2","F":"2","P":"2","Y":"3","M":"3","T":"3","S":"3","H":"4","N":"4","Q":"4","W":"4","R":"5","K":"5","D":"6","E":"6","C":"7"}
    string=list(string)
    for index,item in enumerate(string):
        for key,value in char.items():
            if item==key:
                string[index]=value
    return("".join(string))

def occurrences(string, sub_string):
    count=0
    beg=0
    while(string.find(sub_string,beg)!=-1) :
        count=count+1
        beg=string.find(sub_string,beg)
        beg=beg+1
    return count


def CTC(filename):
    df = pd.DataFrame(columns=['Sequence','Triad:Frequency'])
    data=list((pd.read_csv(filename,sep=',',header=None)).iloc[:,0])
    for i in range(len(data)):
        data1=''
        data1=str(data[i])
        data1=data1.upper()
        allowed = set(('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'))
        is_data_invalid = set(data1).issubset(allowed)
        if is_data_invalid==False:
            print("Errror: Please check for invalid inputs in the sequence.","\nError in: ","Sequence number=",i+1,",","Sequence = ",data[i],",","\nNOTE: Spaces, Special characters('[@_!#$%^&*()<>?/\|}{~:]') and Extra characters(BJOUXZ) should not be there.")
            return
        df.at[i,'Sequence'] = data[i]
        Y.append("".join(repstring(str(data[i]))))
    val2=[[]]
    for f in range(len(LS)):
        val2[0]=val2[0]+["CTC_"+str(LS[f])]
    for j in range(len(data)):
        MM=[]
        for m in range(len(LS)):
            MM=MM+[occurrences(Y[j],LS[m])]
        Min_MM=min(MM)
        Max_MM=max(MM)
        if (Max_MM==0):
            print("Errror: Splits/ Sequence length should be greater than equal to 3")
            return
        val=[]
#         val.append(data[j])
        for k in range(len(LS)):
            val=val+[round(((occurrences(Y[j],LS[k])-Min_MM)/Max_MM),3)]
        val2.append(val)
#     print(val2)
    #file= open(sys.argv[2],'w', newline='')#output file
    file= open(dir_1+'/sam_allcomp.ctc','w', newline='')
    with file:
        writer=csv.writer(file);
        writer.writerows(val2);
    return val2

######################################CTD###################################
def ctd(file):
    attr=pd.read_csv("Data/aa_attr_group.csv", sep="\t")
    filename, file_extension = os.path.splitext(file)
    df1 = pd.read_csv(file, header = None)
    df = pd.DataFrame(df1[0].str.upper())
    n = 0
    stt1 = []
    m = 1
    for i in range(0,len(attr)) :
        st =[]
        stt1.append([])
        for j in range(0,len(df)) :
            stt1[i].append([])
            for k in range(0,len(df[0][j])) :
                while m < 4 :
                    while n < len(attr.iloc[i,m]) :
                        if df[0][j][k] == attr.iloc[i,m][n] :
                            st.append(m)
                            stt1[i][j].append(m)
                        n += 2
                    n = 0
                    m += 1
                m = 1
#####################Composition######################
    f = open(dir_1+"/compout_1", 'w')
    sys.stdout = f
    std = [1,2,3]
    print("1,2,3,")
    for p in range (0,len(df)) :
        for ii in range(0,len(stt1)) :
            #for jj in stt1[ii][p]:
            for pp in std :
                count = 0
                for kk in stt1[ii][p] :
                    temp1 = kk
                    if temp1 == pp :
                        count += 1
                    composition = (count/len(stt1[ii][p]))*100
                print("%.2f"%composition, end = ",")
            print("")
    f.truncate()

#################################Transition#############
    tt = []
    tr=[]
    kk =0
    for ii in range(0,len(stt1)) :
        tt = []
        tr.append([])
        for p in range (0,len(df)) :
            tr[ii].append([])
            while kk < len(stt1[ii][p]) :
                if kk+1 <len(stt1[ii][p]):
                #if  stt1[ii][p][kk] < stt1[ii][p][kk+1] or stt1[ii][p][kk] > stt1[ii][p][kk+1]: # condition for adjacent values
                    tt.append(stt1[ii][p][kk])
                    tt.append(stt1[ii][p][kk+1])
                    tr[ii][p].append(stt1[ii][p][kk])
                    tr[ii][p].append(stt1[ii][p][kk+1])

                kk += 1
            kk = 0

    pp = 0
    xx = []
    xxx = []
    for mm in range(0,len(tr)) :
        xx = []
        xxx.append([])
        for nn in range(0,len(tr[mm])):
            xxx[mm].append([])
            while pp < len(tr[mm][nn]) :
                xx .append(tr[mm][nn][pp:pp+2])
                xxx[mm][nn].append(tr[mm][nn][pp:pp+2])
                pp+=2
            pp = 0

    f1 = open(dir_1+"/compout_2", 'w')
    sys.stdout = f1
    std1 = [[1,1],[1,2],[1,3],[2,1],[2,2],[2,3],[3,1],[3,2],[3,3]]
    print("1->1,1->2,1->3,2->1,2->2,2->3,3->1,3->2,3->3,")
    for rr in range(0,len(df)) :
        for qq in range(0,len(xxx)):
            for tt in std1 :
                count = 0
                for ss in xxx[qq][rr] :
                    temp2 = ss
                    if temp2 == tt :
                        count += 1
                print(count, end = ",")
            print("")
    f1.truncate()

    #################################Distribution#############
    c_11 = []
    c_22 = []
    c_33 = []
    zz = []
    #print("0% 25% 50% 75% 100%")
    for x in range(0,len(stt1)) :
        #c_11.append([])
        c_22.append([])
        #c_33.append([])
        yy_c_1 = []
        yy_c_2 = []
        yy_c_3 = []
        ccc = []

        k = 0
        j = 0
        for y in range(0,len(stt1[x])):
            #c_11[x].append([])
            c_22[x].append([])
            for i in range(1,4) :
                cc = []
                c1 = [index for index,value in enumerate(stt1[x][y]) if value == i]
                c_22[x][y].append(c1)
    cc = []
    for ss in range(0,len(df)):
        for uu in range(0,len(c_22)):
            for mm in range(0,3):
                for ee in range(0,101,25):
                    k = (ee*(len(c_22[uu][ss][mm])))/100
                    cc.append(math.floor(k))
    f2 = open(dir_1+'/compout_3', 'w')
    sys.stdout = f2
    print("0% 25% 50% 75% 100%")
    for i in range (0,len(cc),5):
        print(*cc[i:i+5])
    f2.truncate()
    head = []
    header1 = ['CeTD_HB','CeTD_VW','CeTD_PO','CeTD_PZ','CeTD_CH','CeTD_SS','CeTD_SA']
    for i in header1:
        for j in range(1,4):
            head.append(i+str(j))
    df11 = pd.read_csv(dir_1+"/compout_1")
    df_1 = df11.iloc[:,:-1]
    zz = pd.DataFrame()
    for i in range(0,len(df_1),7):
        zz = zz.append(pd.DataFrame(pd.concat([df_1.loc[i],df_1.loc[i+1],df_1.loc[i+2],df_1.loc[i+3],df_1.loc[i+4],df_1.loc[i+5],df_1.loc[i+6]],axis=0)).transpose()).reset_index(drop=True)
    zz.columns = head
    #zz.to_csv(filename+".ctd_comp", index=None, encoding='utf-8')
    head2 = []
    header2 = ['CeTD_11','CeTD_12','CeTD_13','CeTD_21','CeTD_22','CeTD_23','CeTD_31','CeTD_32','CeTD_33']
    for i in header2:
        for j in ('HB','VW','PO','PZ','CH','SS','SA'):
            head2.append(i+'_'+str(j))
    df12 = pd.read_csv(dir_1+"/compout_2")
    df_2 = df12.iloc[:,:-1]
    ss = pd.DataFrame()
    for i in range(0,len(df_2),7):
        ss = ss.append(pd.DataFrame(pd.concat([df_2.loc[i],df_2.loc[i+1],df_2.loc[i+2],df_2.loc[i+3],df_2.loc[i+4],df_2.loc[i+5],df_2.loc[i+6]],axis=0)).transpose()).reset_index(drop=True)
    ss.columns = head2
    #ss.to_csv(filename+".ctd_trans", index=None, encoding='utf-8')
    head3 = []
    header3 = ['CeTD_0_p','CeTD_25_p','CeTD_50_p','CeTD_75_p','CeTD_100_p']
    header4 = ['HB','VW','PO','PZ','CH','SS','SA']
    for j in range(1,4):
        for k in header4:
            for i in header3:
                head3.append(i+'_'+k+str(j))
    df_3 = pd.read_csv(dir_1+"/compout_3", sep=" ")
    rr = pd.DataFrame()
    for i in range(0,len(df_3),21):
        rr = rr.append(pd.DataFrame(pd.concat([df_3.loc[i],df_3.loc[i+1],df_3.loc[i+2],df_3.loc[i+3],df_3.loc[i+4],df_3.loc[i+5],df_3.loc[i+6],df_3.loc[i+7],df_3.loc[i+8],df_3.loc[i+9],df_3.loc[i+10],df_3.loc[i+11],df_3.loc[i+12],df_3.loc[i+13],df_3.loc[i+14],df_3.loc[i+15],df_3.loc[i+16],df_3.loc[i+17],df_3.loc[i+18],df_3.loc[i+19],df_3.loc[i+20]],axis=0)).transpose()).reset_index(drop=True)
    rr.columns = head3
    cotrdi= pd.concat([zz,ss,rr],axis=1)
    cotrdi.to_csv(dir_1+'/sam_allcomp.ctd', index=None, encoding='utf-8')

###################################qos#######################################
def qos(file,gap,w=0.1):
    ff = []
    filename, file_extension = os.path.splitext(file)
    df = pd.read_csv(file, header = None)
    df2 = pd.DataFrame(df[0].str.upper())
    for i in range(0,len(df2)):
        ff.append(len(df2[0][i]))
    if min(ff) < gap:
        print("Error: All sequences' length should be higher than :", gap)
    else:
        mat1 = pd.read_csv("Data/Schneider-Wrede.csv", index_col = 'Name')
        mat2 = pd.read_csv("Data/Grantham.csv", index_col = 'Name')
        s1 = []
        s2 = []
        for i in range(0,len(df2)):
            for n in range(1, gap+1):
                sum1 = 0
                sum2 = 0
                for j in range(0,(len(df2[0][i])-n)):
                    sum1 = sum1 + (mat1[df2[0][i][j]][df2[0][i][j+n]])**2
                    sum2 = sum2 + (mat2[df2[0][i][j]][df2[0][i][j+n]])**2
                s1.append(sum1)
                s2.append(sum2)
        zz = pd.DataFrame(np.array(s1).reshape(len(df2),gap))
        zz["sum"] = zz.sum(axis=1)
        zz2 = pd.DataFrame(np.array(s2).reshape(len(df2),gap))
        zz2["sum"] = zz2.sum(axis=1)
        c1 = []
        c2 = []
        c3 = []
        c4 = []
        h1 = []
        h2 = []
        h3 = []
        h4 = []
        for aa in std:
            h1.append('QSO'+str(gap)+'_SC_' + aa)
        for aa in std:
            h2.append('QSO'+str(gap)+'_G_' + aa)
        for n in range(1, gap+1):
            h3.append('SC' + str(n))
        h3 = ['QSO'+str(gap)+'_'+sam for sam in h3]
        for n in range(1, gap+1):
            h4.append('G' + str(n))
        h4 = ['QSO'+str(gap)+'_'+sam for sam in h4]
        for i in range(0,len(df2)):
            AA = {}
            for j in std:
                AA[j] = df2[0][i].count(j)
                c1.append(AA[j] / (1 + w * zz['sum'][i]))
                c2.append(AA[j] / (1 + w * zz2['sum'][i]))
            for k in range(0,gap):
                c3.append((w * zz[k][i]) / (1 + w * zz['sum'][i]))
                c4.append((w * zz[k][i]) / (1 + w * zz['sum'][i]))
        pp1 = np.array(c1).reshape(len(df2),len(std))
        pp2 = np.array(c2).reshape(len(df2),len(std))
        pp3 = np.array(c3).reshape(len(df2),gap)
        pp4 = np.array(c4).reshape(len(df2),gap)
        zz5 = round(pd.concat([pd.DataFrame(pp1, columns = h1),pd.DataFrame(pp2,columns = h2),pd.DataFrame(pp3, columns = h3),pd.DataFrame(pp4, columns = h4)], axis = 1),4)
        zz5.to_csv(dir_1+'/sam_allcomp.qso', index = None, encoding = 'utf-8')

#########################Shanon entropy for PCP#################################
def lookup(peptide,featureNum):
    l=len(peptide);
    peptide = list(peptide);
    out=np.zeros(l);
    peptide_num = encode(peptide);
    for i in range(l):
        out[i] = PCP[peptide_num[i]][featureNum];
    return sum(out);
def pcp(file):
    SEP_headers = ['SEP_PC','SEP_NC','SEP_NE','SEP_PO','SEP_NP','SEP_AL','SEP_CY','SEP_AR','SEP_AC','SEP_BS','SEP_NE_pH','SEP_HB','SEP_HL','SEP_NT','SEP_HX','SEP_SC','SEP_SS_HE','SEP_SS_ST','SEP_SS_CO','SEP_SA_BU','SEP_SA_EX','SEP_SA_IN','SEP_TN','SEP_SM','SEP_LR']
    if(type(file) == str):
        seq = pd.read_csv(file,header=None, sep=',');
        seq=seq.T
        seq[0].values.tolist()
        seq=seq[0];
    else:
        seq  = file;
    l = len(seq);
    rows = PCP.shape[0]; # Number of features in our reference table
    col = 20 ; # Denotes the 20 amino acids
    seq=[seq[i].upper() for i in range(l)]
    sequenceFeature = [];
    sequenceFeature.append(SEP_headers); #To put property name in output csv
    
    for i in range(l): # Loop to iterate over each sequence
        nfeatures = rows;
        sequenceFeatureTemp = [];
        for j in range(nfeatures): #Loop to iterate over each feature
            featureVal = lookup(seq[i],j)   
            if(len(seq[i])!=0):
                sequenceFeatureTemp.append(featureVal/len(seq[i]))
            else:
                sequenceFeatureTemp.append('NaN')
        sequenceFeature.append(sequenceFeatureTemp);
    out = pd.DataFrame(sequenceFeature);
    return sequenceFeature;
def phyChem(file,mode='all',m=0,n=0):
    if(type(file) == str):
        seq1 = pd.read_csv(file,header=None, sep=',');
        seq1 = pd.DataFrame(seq1[0].str.upper())
        seq=[]
        [seq.append(seq1.iloc[i][0]) for i in range(len(seq1))]
    else:
        seq  = file;
    l = len(seq);
    newseq = [""]*l; # To store the n-terminal sequence
    for i in range(0,l):
        l = len(seq[i]);
        if(mode=='NT'):
            n=m;
            if(n!=0):
                newseq[i] = seq[i][0:n];
            elif(n>l):
                print('Warning! Sequence',i,"'s size is less than n. The output table would have NaN for this sequence");
            else:
                print('Value of n is mandatory, it cannot be 0')
                break;
        elif(mode=='CT'):
            n=m;
            if(n!=0):
                newseq[i] = seq[i][(len(seq[i])-n):]
            elif(n>l):
                print('WARNING: Sequence',i+1,"'s size is less than the value of n given. The output table would have NaN for this sequence");
            else:
                print('Value of n is mandatory, it cannot be 0')
                break;
        elif(mode=='all'):
            newseq = seq;
        elif(mode=='rest'):
            if(m==0):
                print('Kindly provide start index for rest, it cannot be 0');
                break;
            else:
                if(n<=len(seq[i])):
                    newseq[i] = seq[i][m-1:n+1]
                elif(n>len(seq[i])):
                    newseq[i] = seq[i][m-1:len(seq[i])]
                    print('WARNING: Since input value of n for sequence',i+1,'is greater than length of the protein, entire sequence starting from m has been considered')
        else:
            print("Wrong Mode. Enter 'NT', 'CT','all' or 'rest'");        
    output = pcp(newseq);
    return output


def shannons(filename):
    SEP_headers = ['SEP_PC','SEP_NC','SEP_NE','SEP_PO','SEP_NP','SEP_AL',
                   'SEP_CY','SEP_AR','SEP_AC','SEP_BS','SEP_NE_pH','SEP_HB',
                   'SEP_HL','SEP_NT','SEP_HX','SEP_SC','SEP_SS_HE','SEP_SS_ST',
                   'SEP_SS_CO','SEP_SA_BU','SEP_SA_EX','SEP_SA_IN','SEP_TN','SEP_SM','SEP_LR']
    if(type(filename) == str):
        seq1 = pd.read_csv(filename,header=None, sep=',');
        seq1 = pd.DataFrame(seq1[0].str.upper())
    else:
        seq1  = filename;
    seq=[]
    [seq.append(seq1.iloc[i][0]) for i in range(len(seq1))]
    comp = phyChem(seq);
    new = [comp[i][0:25] for i in range(len(comp))]
    entropy  = [];
    entropy.append(SEP_headers[0:25])
    for i in range(1,len(new)):
        seqEntropy = [];
        for j in range(len(new[i])):
            p = new[i][j]; 
            if((1-p) == 0. or p==0.):
                temp = 0;#to store entropy of each sequence
            else:
                temp = -(p*math.log2(p)+(1-p)*math.log2(1-p));
            seqEntropy.append(round(temp,3));
        entropy.append(seqEntropy);

    file = open(dir_1+'/sam_allcomp.sep','w')
    with file:
        writer = csv.writer(file);
        writer.writerows(entropy);
    return entropy;	


def feature_gen(file):
    aac_comp(file)
    dpc_comp(file)
    tpc_comp(file)
    pcp_1(file)
    DDOR(file)
    CTC(file)
    SE_residue_level(file)
    SE(file)
    ctd(file)
    qos(file, 1)
    shannons(file)

    df1 = pd.read_csv(os.path.join(dir_1, "sam_allcomp.aac"))
    df2 = pd.read_csv(os.path.join(dir_1, "sam_allcomp.dpc"))
    df3 = pd.read_csv(os.path.join(dir_1, "sam_allcomp.tpc"))
    df6 = pd.read_csv(os.path.join(dir_1, "sam_allcomp.pcp"))
    df9 = pd.read_csv(os.path.join(dir_1, "sam_allcomp.ddr"))
    df10 = pd.read_csv(os.path.join(dir_1, "sam_allcomp.sep"))
    df11 = pd.read_csv(os.path.join(dir_1, "sam_allcomp.ser"))
    df13 = pd.read_csv(os.path.join(dir_1, "sam_allcomp.ctc"))
    df14 = pd.read_csv(os.path.join(dir_1, "sam_allcomp.ctd"))
    df17 = pd.read_csv(os.path.join(dir_1, "sam_allcomp.qso"))

    df19 = pd.concat([df1.iloc[:,:-1],df2.iloc[:,:-1],df3.iloc[:,:-1],df6,df9.iloc[:,:-1],df10,df11,df13,df14,df17],axis=1)
    return df19
    # df19.to_csv(features.csv, index=None)
    # filelist=glob.glob(dir_1+"/sam_allcomp*")
    # for file_2 in filelist:
    #     os.remove(file_2)
    # For Top Features

def top_feat(file):
    df5 = pd.read_csv(file)
    feat = ['DPC1_CF',	'TPC_AQI',	'TPC_CFN',	'TPC_CVL',	'TPC_DLN',	'TPC_DML',	'TPC_EQS',	'TPC_FRP',	'TPC_GCK',	'TPC_GLM',	'TPC_GNF',	'TPC_HLC',	'TPC_KCC',	'TPC_KYS',	'TPC_LAN',	'TPC_LGM',	'TPC_LLF',	'TPC_LMG',	'TPC_LNS',	'TPC_MAY',	'TPC_NTP',	'TPC_RGL',	'TPC_RKY',	'TPC_RRP',	'TPC_SIR',	'TPC_THR',	'TPC_TTG',	'TPC_VCG',	'TPC_VSF',	'PCP_AR',	'DDR_C',	'DDR_K',	'SER_T',	'SEP_LR',	'CTC_174',	'CTC_374',	'CTC_477',	'CTC_677',	'CTC_713',	'CeTD_11_PO',	'CeTD_11_SS',	'CeTD_11_SA',	'CeTD_12_HB',	'CeTD_13_SA',	'CeTD_21_HB',	'CeTD_21_PZ',	'CeTD_23_PO',	'CeTD_100_p_PO1',	'CeTD_100_p_PZ3',	'QSO1_SC_V']
    df6 = df5[feat]
    return df6

def prediction(inputfile,model,out):
    df = pd.DataFrame()
    a=[]
    file_name = inputfile
    file_name1 = out
    file_name2 = model
    clf = joblib.load(file_name2)
    data_test = np.loadtxt(file_name, delimiter=',')
    X_test = data_test
    y_p_score1=clf.predict_proba(X_test)
    y_p_s1=y_p_score1.tolist()
    df = pd.DataFrame(y_p_s1)
    df_1 = df.iloc[:,-1]
    df_1.to_csv(file_name1, index=None, header=False)

def class_assignment(file1,thr,out):
    df1 = pd.read_csv(file1, header=None)
    df1.columns = ['ML Score']
    cc = []
    for i in range(0,len(df1)):
        if df1['ML Score'][i]>=float(thr):
            cc.append('Peptide Hormone')
        else:
            cc.append('Non-Hormonal Peptide')
    df1['Prediction'] = cc
    df1 =  df1.round(3)
    df1.to_csv(out, index=None)

def MERCI_Processor(merci_file,merci_processed,name):
    hh =[]
    jj = []
    kk = []
    qq = []
    filename = merci_file
    df = pd.DataFrame(name)
    zz = list(df[0])
    check = '>'
    with open(filename) as f:
        l = []
        for line in f:
            if not len(line.strip()) == 0 :
                l.append(line)
            if 'COVERAGE' in line:
                for item in l:
                    if item.lower().startswith(check.lower()):
                        hh.append(item)
                l = []
    if hh == []:
        ff = [w.replace('>', '') for w in zz]
        for a in ff:
            jj.append(a)
            qq.append(np.array(['0']))
            kk.append('Non-Allergen')
    else:
        ff = [w.replace('\n', '') for w in hh]
        ee = [w.replace('>', '') for w in ff]
        rr = [w.replace('>', '') for w in zz]
        ff = ee + rr
        oo = np.unique(ff)
        df1 = pd.DataFrame(list(map(lambda x:x.strip(),l))[1:])
        df1.columns = ['Name']
        df1['Name'] = df1['Name'].str.strip('(')
        df1[['Seq','Hits']] = df1.Name.str.split("(",expand=True)
        df2 = df1[['Seq','Hits']]
        df2.replace(to_replace=r"\)", value='', regex=True, inplace=True)
        df2.replace(to_replace=r'motifs match', value='', regex=True, inplace=True)
        df2.replace(to_replace=r' $', value='', regex=True,inplace=True)
        total_hit = int(df2.loc[len(df2)-1]['Seq'].split()[0])
        for j in oo:
            if j in df2.Seq.values:
                jj.append(j)
                qq.append(df2.loc[df2.Seq == j]['Hits'].values)
                kk.append('Allergen')
            else:
                jj.append(j)
                qq.append(np.array(['0']))
                kk.append('Non-Allergen')
    df3 = pd.concat([pd.DataFrame(jj),pd.DataFrame(qq),pd.DataFrame(kk)], axis=1)
    df3.columns = ['Name','Hits','Prediction']
    df3.to_csv(merci_processed,index=None)

def Merci_after_processing(merci_processed,final_merci):
    df5 = pd.read_csv(merci_processed)
    df5 = df5[['Name','Hits']]
    df5.columns = ['Subject','Hits']
    kk = []
    for i in range(0,len(df5)):
        if df5['Hits'][i] > 0:
            kk.append(0.5)
        else:
            kk.append(0)
    df5["MERCI Score"] = kk
    df5 = df5[['Subject','MERCI Score']]
    df5.to_csv(final_merci, index=None)

def BLAST_processor(blast_result,blast_processed,name1):
    if os.stat(blast_result).st_size != 0:
        df1 = pd.read_csv(blast_result, sep="\t",header=None)
        df2 = df1.iloc[:,:2]
        df2.columns = ['Subject','Query']
        df3 = pd.DataFrame()
        for i in df2.Subject.unique():
            df3 = df3.append(df2.loc[df2.Subject==i][0:5]).reset_index(drop=True)
        cc= []
        for i in range(0,len(df3)):
            cc.append(df3['Query'][i].split("_")[0])
        df3['label'] = cc
        dd = []
        for i in range(0,len(df3)):
            if df3['label'][i] == 'P':
                dd.append(1)
            else:
                dd.append(-1)
        df3["vote"] = dd
        ff = []
        gg = []
        for i in df3.Subject.unique():
            ff.append(i)
            gg.append(df3.loc[df3.Subject==i]["vote"].sum())
        df4 = pd.concat([pd.DataFrame(ff),pd.DataFrame(gg)],axis=1)
        df4.columns = ['Subject','Blast_value']
        hh = []
        for i in range(0,len(df4)):
            if df4['Blast_value'][i] >0:
                hh.append(0.5)
            elif df4['Blast_value'][i] == 0:
                hh.append(0)
            else:
                hh.append(-0.5)
        df4['BLAST Score'] = hh
        df4 = df4[['Subject','BLAST Score']]
    else:
        ss = []
        vv = []
        for j in seqid:
            ss.append(j)
            vv.append(0)
        df4 = pd.concat([pd.DataFrame(ss),pd.DataFrame(vv)],axis=1)
        df4.columns = ['Subject','BLAST Score']
    df4.to_csv(blast_processed, index=None)

def hybrid(ML_output,name1,merci_output,blast_output,threshold,final_output):
    df6_2 = pd.read_csv(ML_output,header=None)
    df6_1 = pd.DataFrame(name1)
    df5 = pd.read_csv(merci_output)
    df4 = pd.read_csv(blast_output)
    df6 = pd.concat([df6_1,df6_2],axis=1)
    df6.columns = ['Subject','ML Score']
    df6['Subject'] = df6['Subject'].str.replace('>','')
    df7 = pd.merge(df6,df5, how='outer',on='Subject')
    df8 = pd.merge(df7,df4, how='outer',on='Subject')
    df8.fillna(0, inplace=True)
    df8['Hybrid Score'] = df8.sum(axis=1)
    df8 = df8.round(3)
    ee = []
    for i in range(0,len(df8)):
        if df8['Hybrid Score'][i] > float(threshold):
            ee.append('Allergen')
        else:
            ee.append('Non-Allergen')
    df8['Prediction'] = ee
    df8.to_csv(final_output, index=None)

print('##############################################################################')
print('# The program Hoppred is developed for predicting peptide hormones #')
print("# protein from their primary sequence, developed by Prof G. P. S. Raghava's group. #")
print('# ############################################################################')

# Parameter initialization or assigning variable for command level arguments

Sequence= args.input        # Input variable 
 
# Output file 
 
if args.output == None:
    result_filename= "outfile.csv" 
else:
    result_filename = args.output
         
# Threshold 
if args.threshold == None:
        Threshold = 0.5
else:
        Threshold= float(args.threshold)
# Model
if args.model == None:
        Model = int(1)
else:
        Model = int(args.model)
# Display
if args.display == None:
        dplay = int(1)
else:
        dplay = int(args.display)

print('Summary of Parameters:')
print('Input File: ',Sequence,'; Model: ',Model,'; Threshold: ', Threshold)
print('Output File: ',result_filename,'; Display: ',dplay)

#------------------ Read input file ---------------------
f=open(Sequence,"r")
len1 = f.read().count('>')
f.close()

with open(Sequence) as f:
        records = f.read()
records = records.split('>')[1:]
seqid = []
seq = []
for fasta in records:
    array = fasta.split('\n')
    name, sequence = array[0].split()[0], re.sub('[^ARNDCQEGHILKMFPSTWYV-]', '', ''.join(array[1:]).upper())
    seqid.append(name)
    seq.append(sequence)
if len(seqid) == 0:
    f=open(Sequence,"r")
    data1 = f.readlines()
    for each in data1:
        seq.append(each.replace('\n',''))
    for i in range (1,len(seq)+1):
        seqid.append("Seq_"+str(i))

seqid_1 = list(map(">{}".format, seqid))
CM = pd.concat([pd.DataFrame(seqid_1),pd.DataFrame(seq)],axis=1)
CM.to_csv("Sequence_1",header=False,index=None,sep="\n")
f.close()

#======================= Prediction Module start from here =====================
if Model==1:
    aac_comp(seq,'seq.aac')
    os.system("perl -pi -e 's/,$//g' seq.aac")
    prediction('seq.aac','rf_model','seq.pred')
    class_assignment('seq.pred',Threshold,'seq.out')
    df1 = pd.DataFrame(seqid)
    df2 = pd.DataFrame(seq)
    df3 = pd.read_csv("seq.out")
    df3 = round(df3,3)
    df4 = pd.concat([df1,df2,df3],axis=1)
    df4.columns = ['ID','Sequence','ML_Score','Prediction']
    if dplay == 1:
        df4 = df4.loc[df4.Prediction=="Allergen"]
    else:
        df4 = df4
    df4.to_csv(result_filename, index=None)
    os.remove('seq.aac')
    os.remove('seq.pred')
    os.remove('seq.out')
else:
    if os.path.exists('envfile'):
        with open('envfile', 'r') as file:
            data = file.readlines()
        output = []
        for line in data:
            if not "#" in line:
                output.append(line)
        if len(output)==4: 
            paths = []
            for i in range (0,len(output)):
                paths.append(output[i].split(':')[1].replace('\n',''))
            blastp = paths[0]
            blastdb = paths[1]
            merci = paths[2]
            motifs = paths[3]
        else:
            print("####################################################################################")
            print("Error: Please provide paths for BLAST, MERCI and required files", file=sys.stderr)
            print("####################################################################################")
            sys.exit()
 
    else:
        print("####################################################################################")
        print("Error: Please provide the '{}', which comprises paths for BLAST and MERCI".format('envfile'), file=sys.stderr)
        print("####################################################################################")
        sys.exit()
    aac_comp(seq,'seq.aac')
    os.system("perl -pi -e 's/,$//g' seq.aac")
    prediction('seq.aac','rf_model','seq.pred')
    os.system(blastp + " -task blastp -db " + blastdb + " -query " + "Sequence_1"  + " -out RES_1_6_6.out -outfmt 6 -evalue 0.1")
    os.system(merci + " -p " + "Sequence_1" +  " -i " + motifs + " -o merci.txt")
    MERCI_Processor('merci.txt','merci_output.csv',seqid)
    Merci_after_processing('merci_output.csv','merci_hybrid.csv')
    BLAST_processor('RES_1_6_6.out','blast_hybrid.csv',seqid)
    hybrid('seq.pred',seqid,'merci_hybrid.csv','blast_hybrid.csv',Threshold,'final_output')
    df44 = pd.read_csv('final_output')
    if dplay == 1:
        df44 = df44.loc[df44.Prediction=="Allergen"]
    else:
        df44 = df44
    df44 = round(df44,3)
    df44.to_csv(result_filename, index=None)
    os.remove('seq.aac')
    os.remove('seq.pred')
    os.remove('final_output')
    os.remove('RES_1_6_6.out')
    os.remove('merci_output.csv')
    os.remove('merci_hybrid.csv')
    os.remove('blast_hybrid.csv')
    os.remove('merci.txt')
    os.remove('Sequence_1')

print('\n======= Thanks for using hoppred2.0. Your results are stored in file :',result_filename,' =====\n\n')
print('Please cite: hoppred2.0\n')

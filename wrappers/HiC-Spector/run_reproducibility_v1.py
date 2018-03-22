import sys
import scipy
import numpy
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import eigsh
import gzip

def Parse_matrix(file1,file2):
  max_index=0
  max_index_temp=0
  with open(file1) as input_file:
    for line in input_file:
        x,y,z=map(int, line.split())
        max_index_temp=max(x,y)
        if max_index_temp>max_index:
               max_index=max_index_temp
  
  with open(file2) as input_file:
    for line in input_file:
        x,y,z=map(int, line.split())
        max_index_temp=max(x,y)
        if max_index_temp>max_index:
               max_index=max_index_temp


  M1=lil_matrix((max_index,max_index))
  M2=lil_matrix((max_index,max_index))
  with open(file1) as input_file:
    for line in input_file:
        x,y,z=map(int, line.split())
        M1[x-1,y-1]=z
        M1[y-1,x-1]=z
  with open(file2) as input_file:
    for line in input_file:
        x,y,z=map(int, line.split())
        M2[x-1,y-1]=z
        M2[y-1,x-1]=z
  return M1, M2 


def Parse_matrix_lieberman(file1,file2, resolution):
  max_index=0
  max_index_temp=0
  with gzip.open(file1) as input_file:
    for line in input_file:
      if line[0]!="#":  
        x,y,z=map(float, line.split())
        x=int(x)
        y=int(y)
        z=int(z)
        max_index_temp=max(x,y)
        if max_index_temp>max_index:
               max_index=max_index_temp
  
  with gzip.open(file2) as input_file:
    for line in input_file:
      
      if line[0]!="#":  
        x,y,z=map(float, line.split())
        x=int(x)
        y=int(y)
        z=int(z)
        max_index_temp=max(x,y)
        if max_index_temp>max_index:
               max_index=max_index_temp

  max_index=max_index/resolution+1
  
  M1=lil_matrix((max_index,max_index))
  M2=lil_matrix((max_index,max_index))
  with gzip.open(file1) as input_file:
    for line in input_file:
      if line[0]!="#":  
        x,y,z=map(float, line.split())
        x=int(x)
        y=int(y)
        z=int(z)
        M1[x/resolution,y/resolution]=z
        M1[y/resolution,x/resolution]=z
  with gzip.open(file2) as input_file:
    for line in input_file:
      if line[0]!="#":  
        x,y,z=map(float, line.split())
        x=int(x)
        y=int(y)
        z=int(z)
        M2[x/resolution,y/resolution]=z
        M2[y/resolution,x/resolution]=z
  return M1, M2 


def Parse_matrix_hic(file1, file2, chrn, resolution):
    
    Table1=straw.straw("NONE",file1, chrn, chrn,"BP",resolution)
    Table2=straw.straw("NONE",file2, chrn, chrn,"BP",resolution)
    max_index=max(max(Table1[0]),max(Table1[1]),max(Table2[0]),max(Table2[1]))
    max_index=max_index/resolution
    M1=lil_matrix((max_index+1,max_index+1))
    M2=lil_matrix((max_index+1,max_index+1))
                     
    for i in range(len(Table1[0])):
        M1[Table1[0][i]/resolution,Table1[1][i]/resolution]=Table1[2][i]
        M1[Table1[1][i]/resolution,Table1[0][i]/resolution]=Table1[2][i]
                                                             
    for i in range(len(Table2[0])):
        M2[Table2[0][i]/resolution,Table2[1][i]/resolution]=Table2[2][i]
        M2[Table2[1][i]/resolution,Table2[0][i]/resolution]=Table2[2][i]
    return M1, M2


def  get_Laplacian(M):
     S=M.sum(1)
     i_nz=numpy.where(S>0)[0]
     S=S[i_nz]
     M=(M[i_nz].T)[i_nz].T
     S=1/numpy.sqrt(S)
     M=S*M
     M=(S*M.T).T
     n=numpy.size(S)
     M=numpy.identity(n)-M
     M=(M+M.T)/2
     return M

def evec_distance(v1,v2):
    d1=numpy.dot(v1-v2,v1-v2)
    d2=numpy.dot(v1+v2,v1+v2)
    if d1<d2:
         d=d1
    else:
        d=d2
    return numpy.sqrt(d)

def get_ipr(evec):
      ipr=1.0/(evec*evec*evec*evec).sum()
      return ipr


def get_reproducibility(M1,M2,num_evec):
   k1=numpy.sign(M1.A).sum(1)
   d1=numpy.diag(M1.A)
   kd1=~((k1==1)*(d1>0))
   k2=numpy.sign(M2.A).sum(1)
   d2=numpy.diag(M2.A)
   kd2=~((k2==1)*(d2>0))
   iz=numpy.nonzero((k1+k2>0)*(kd1>0)*(kd2>0))[0]
   M1b=(M1[iz].A.T)[iz].T
   M2b=(M2[iz].A.T)[iz].T

   i_nz1=numpy.where(M1b.sum(1)>0)[0]
   i_nz2=numpy.where(M2b.sum(1)>0)[0]
   i_z1=numpy.where(M1b.sum(1)==0)[0]
   i_z2=numpy.where(M2b.sum(1)==0)[0]
   
   M1b_L=get_Laplacian(M1b)
   M2b_L=get_Laplacian(M2b)
   
   a1, b1=eigsh(M1b_L,k=num_evec,which="SM")
   a2, b2=eigsh(M2b_L,k=num_evec,which="SM")
   
   b1_extend=numpy.zeros((numpy.size(M1b,0),num_evec))
   b2_extend=numpy.zeros((numpy.size(M2b,0),num_evec))
   for i in range(num_evec):
       b1_extend[i_nz1,i]=b1[:,i]
       b2_extend[i_nz2,i]=b2[:,i]
   
   #ipr_cut=5
   #ipr1=numpy.zeros(num_evec)
   #ipr2=numpy.zeros(num_evec)
   #for i in range(num_evec):
   #    ipr1[i]=get_ipr(b1_extend[:,i])
   #    ipr2[i]=get_ipr(b2_extend[:,i])
  
   #b1_extend_eff=b1_extend[:,ipr1>ipr_cut]
   #b2_extend_eff=b2_extend[:,ipr2>ipr_cut]
   #num_evec_eff=min(numpy.size(b1_extend_eff,1),numpy.size(b2_extend_eff,1))
   b1_extend_eff=b1_extend
   b2_extend_eff=b2_extend
   num_evec_eff=num_evec
  
   evd=numpy.zeros(num_evec_eff)
   for i in range(num_evec_eff):
       evd[i]=evec_distance(b1_extend_eff[:,i],b2_extend_eff[:,i])
   
   Sd=evd.sum()
   l=numpy.sqrt(2)
   evs=abs(l-Sd/num_evec_eff)/l
   #print("size of maps: %d" %(numpy.size(M1,0)))
   print("reproducibility score: %6.3f " %(evs))
   #print("num_evec_eff: %d" %(num_evec_eff))

def main():
    #num_evec=20;
    if len(sys.argv)==4 and sys.argv[1]=="-F":
        M1, M2=Parse_matrix(sys.argv[2],sys.argv[3])
        get_reproducibility(M1,M2,num_evec)
    elif len(sys.argv)==6 and sys.argv[1]=="-f":
        M1, M2=Parse_matrix_hic(sys.argv[2],sys.argv[3],sys.argv[4],int(sys.argv[5]))
        get_reproducibility(M1,M2,num_evec)
    elif len(sys.argv)==7 and sys.argv[1]=="t":
        M1, M2=Parse_matrix_lieberman(sys.argv[2],sys.argv[3],int(sys.argv[5]))
        sys.stdout = open(sys.argv[4], 'w')
        get_reproducibility(M1,M2,int(sys.argv[6]))
    else:
      print('To use lieberman matrix files as the input:')
      print('python run_reproducibility.py t matrix_file1 matrix_file2 out resolution[int] num_evec')


if __name__ == '__main__':
    main()

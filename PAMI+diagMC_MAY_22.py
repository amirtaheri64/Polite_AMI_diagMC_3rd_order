'''
Polite AMI+diagMC for the
self-energy up to the second order

Last updates: MAY-22-2019
'''

###################################
##########Import libraries#########
###################################
#Begin
from random import random
from random import randint
import pickle
from math import *
import numpy as np
from Symbolic_multi_AMI_new import *
#End
###################################

###################################
######Proposal probabilities#######
###################################
#Begin
def W(ORDER): 
  if ORDER==1:
    return 1.0
  else:
    val=1.0/(2*pi)**(2*ORDER)
  return val
#End
###################################

###################################
###########Update_order############
###################################
#Begin
def update_order(ORDER):  # Change order by 1
  r=random()
  #if flag_re==1:
  if r<1.0/2.0: 
    val=ORDER-1
  if r>=1.0/2.0:
    val=ORDER+1
  return val
  #if flag_re==0:
    #if r<1.0/2.0: 
      #val=ORDER-1
    #if r>=1.0/2.0:
      #val=ORDER+1
    #return val
#End
###################################

###################################
###########Update_NUMBER###########
###################################
#Begin
def update_NUMBER(ORDER):
  if ORDER==1 or ORDER==2:
    return 1
  if ORDER==3:
    return 1
  if ORDER==4:
    return randint(1,3)
#End
###################################  

###################################
#########Number of terms###########
###################################
#Begin
def num_terms(ORDER): 
  if ORDER==1 or ORDER==2:
    return 1.0
  if ORDER==3:
    return 1.0
  if ORDER==4:
    return 3.0
#End
###################################
  

###################################
############Update momenta#########
###################################
#Begin
def Integrand(M,NUMBER):   # For now number=1
  out=[]  # To store output
  
  ##################################
  ##########First order term########
  ##################################
  #Begin
  if M==1:   # First order term is our Normalization
    out.append(1.0)
    out.append(1.0)
    return 1.0
  #if M==0 or M==3:
    #out.append(0)
    #out.append(0)
    #return out
  #End
  ###################################

  ###################################
  ##########Read AMI arrays##########
  ###################################
  #Begin
  file_name_S = 'S_m_'+str(M)+'_num_'+str(NUMBER)
  with open(file_name_S,"rb") as fp:
    s_list=pickle.load(fp)
  #print 'S_list = ', s_list

  file_name_P_freq = 'P_freq_m_'+str(M)+'_num_'+str(NUMBER)
  with open(file_name_P_freq,"rb") as fp:
    p_list_freq=pickle.load(fp)
  #print 'P_list_freq = ', p_list_freq

  file_name_P_mnta = 'P_mnta_m_'+str(M)+'_num_'+str(NUMBER)
  with open(file_name_P_mnta,"rb") as fp:
    p_list_mnta=pickle.load(fp)
  #print 'P_list_mnta = ', p_list_mnta

  file_name_R_freq = 'R_freq_m_'+str(M)+'_num_'+str(NUMBER)
  with open(file_name_R_freq,"rb") as fp:
    r_freq=pickle.load(fp)
  #print 'R_freq = ', r_freq

  file_name_R_mnta = 'R_mnta_m_'+str(M)+'_num_'+str(NUMBER)
  with open(file_name_R_mnta,"rb") as fp:
    r_mnta=pickle.load(fp)
  #print 'R_mnta = ', r_mnta

  file_name_G_sym = 'G_sym_m_'+str(M)+'_num_'+str(NUMBER)
  with open(file_name_G_sym,"rb") as fp:
    g_sym=pickle.load(fp)
  #print 'G_sym = ', g_sym
  #End
  ###################################

  ###################################
  #######Read number of loops########
  ###################################
  #Begin
  file_name_F='f_m_'+str(M)+'_num_'+str(NUMBER)
  with open(file_name_F,"rb") as fp:
    F=pickle.load(fp)  # F: number of loops
  #End
  ###################################
  
  ###################################
  #####Initialize moemntum array#####
  ###################################
  #Begin
  k_init=[ [0]*2 for i in range (0,M+1) ]  # initialize momenta
  k_init[M][0] = px_ext
  k_init[M][1] = py_ext
  #print k_init
  #End
  ###################################
  
  ###################################
  ############AMI algebra############ 
  ###################################
  #Begin
  for i in range (0,M):  # Pick momenta randomly
    for j in range (0,2):
      k_init[i][j]= (2*pi)*random()-pi
  #k_init[0][0]=0
  #k_init[0][1]=0
  #k_init[1][0]=pi/4.0
  #k_init[1][1]=pi/3.0
  #k_init[2][0]=pi
  #k_init[2][1]=pi/3.0
  #print 'k_init = ', k_init
  poles = Poles(p_list_freq, p_list_mnta,g_sym[0],k_init,t,mu)
  for i in range (0,len(poles)):
    for j in range (0,len(poles[i])):
      if poles[i][j]==None:
        j = len(poles[i])+1
      else:
        for l in range (0,len(poles[i][j])):
          #print poles[i][j][k]
          poles[i][j][l] = f(poles[i][j][l],beta)  # Find the numerical values 
  #print 'poles = ', poles

  ########################################
  ################AMI Algebra#############
  ########################################

  res1 = dot_num(s_list[1][0], poles[1][0])
  temp=res1
  for it_m in range (2,M+1):
    o = 0
    for i in range(0, len(s_list[it_m])):
      if s_list[it_m][i] != None:
        o = o + 1
      else:
        i = len(poles[it_m]) + 1
    
    S_new = [None]*(o)
    P_new = [None]*(o)
    for i in range(0, len(s_list[it_m])):
      if s_list[it_m][i] != None:
        S_new[i] = s_list[it_m][i]
        P_new[i] = poles[it_m][i]
  
      else:
        i = len(s_list[it_m]) + 1

    res2 = dot_arr(S_new, P_new)
    temp = cross(res1,res2)
    res1 = temp
  
  G_val = 0.0  # To store numerical value at a given external frequency
  #if len(r_freq[0])!=0:
  for j in range (0, len(temp)):
    G_val = temp[j]*G_eval(r_freq[j], r_mnta[j], nu_ext, g_sym[0],k_init,t,mu,beta) + G_val
    #print 'G_val = ', G_val
    #print '?????'
  #else:
    #for j in range (0, len(temp)):
      #G_val = temp[j] + G_val 
    #print '?????'   
    #print 'G_val = ', G_val
  #print -G_val/4/4/4/pi/pi/pi/pi/pi/pi
  out.append((-1)**(M+F)*(U**M)*G_val.real/4**M/pi**(2*M))
  out.append((-1)**(M+F)*(U**M)*G_val.imag/4**M/pi**(2*M))
 
  #print out
  return out[flag_imag]
      
#End
###################################

###################################
#########Acceptance ratio##########
###################################
#Begin
def acp_rt(M_P,M_O,NUMBER_P,NUMBER_O,int_P,int_O):  # Return the acceptance ratio
  sign_O=np.sign(int_O)
  if M_P>=(highest_order+1) or M_P<=0:  # If the order is out of range
    val=0.0
    return val,sign_O,0
  sign_P=np.sign(int_P)
  if int_O==0:
    val=1.0
    return val,sign_O,sign_P
  p = num_terms(M_P)/num_terms(M_O)
  val=p*abs(int_P)*W(M_O)/abs(int_O)/W(M_P)
  return val,sign_O,sign_P
#End
###################################

###################################
##########Update function##########
###################################
#Begin
def new_state(int_O,M_O,NUMBER_O):  # Accept or reject the proposed state
  M_P = update_order(M_O)
  if M_P<=0 or M_P>=(highest_order+1):
    return M_O, NUMBER_O, np.sign(int_O), int_O
  NUMBER_P = update_NUMBER(M_P)
  int_P=Integrand(M_P,NUMBER_P)
  R_signs=[]
  R_signs.append(acp_rt(M_P,M_O,NUMBER_P,NUMBER_O,int_P,int_O))
  R=R_signs[0][0]
  r=random()
  if R<= r:
    M_N = M_O
    NUMBER_N = NUMBER_O
    sign_new=R_signs[0][1]
    int_N=int_O
    return M_N,NUMBER_N,sign_new,int_N
  else:  
    M_N = M_P
    NUMBER_N = NUMBER_P 
    sign_new=R_signs[0][2]
    int_N=int_P
    return M_N,NUMBER_N,sign_new,int_N
#End
###################################

###################################
############Update_func############
###################################
#Begin
def update(M_O,NUMBER_O):
  x1=[]  # To store Markov's chain
  x2=[]  # To store Markov's chain
  x3=[]  # To store Markov's chain
  x4=[]
  for i_T in range (0,T):
    N1=0.0
    N2=0.0
    N3=0.0
    N4=0.0
    next_state=[]
    int_O=Integrand(M_O,NUMBER_O)
    for i_N in range (0,N):  # To propose N states 
      #print M_O
      
      next_state=new_state(int_O,M_O,NUMBER_O)[:] # Update 
      M_O=next_state[0]
      #print 'sign = ', next_state[0]
      #print next_state[0][2]
      NUMBER_O=next_state[1]  
      int_O=next_state[3]
      #print M_O
      if i_T>=0:
        if M_O == 1: # To compute N1
          N1 = N1 + next_state[2]
        if M_O == 2: # To compute N2
          N2 = N2 + next_state[2]
        if M_O == 3: # To compute N3
          N3 = N3 + next_state[2]
        if M_O == 4:
          N4 = N4 + next_state[2]
      next_state=[]
    M_O=1
    NUMBER_O=1
    print 'a2 = ', N2/N1
    print 'a3 = ', N3/N1
    print 'a4 = ', N4/N1
    x1.append(N1)
    x2.append(N2)
    x3.append(N3)
    x4.append(N4)
    
  return x1,x2,x3,x4
#End
###################################

###################################
###############Error###############
###################################
#Begin
def measurement(x1,x2,x3,x4):  # xi: Markov's chain corresponding to Ni
  var1=0.0  
  var2=0.0
  var3=0.0
  var4=0.0
  avg_1=0.0
  avg_2=0.0
  avg_3=0.0
  avg_4=0.0
  error_N2_N1=0.0
  error_N3_N1=0.0
  error_N4_N1=0.0

  for i in range (0,len(x1)):
    avg_1 = avg_1 + x1[i]
  avg_1=avg_1/len(x1)
  for i in range (0,len(x2)):
    avg_2 = avg_2 + x2[i]
  avg_2=avg_2/len(x2)
  for i in range (0,len(x3)):
    avg_3 = avg_3 + x3[i]
  avg_3=avg_3/len(x3)
  for i in range (0,len(x4)):
    avg_4 = avg_4 + x4[i]
  avg_4=avg_4/len(x4)
  #print 'avg_0= ', avg_0
  #print 'avg_1= ', avg_1
  #print 'avg_2= ', avg_2

  #var0=np.std(x0)
  #var0 = var0/sqrt(MC_iter)
  #var1=np.std(x1)
  #var1 = var1/sqrt(MC_iter)
  #var2=np.std(x2)
  #var2 = var0/sqrt(MC_iter)
  for i in range (0,T):  # Calculate the variance
    var1 = var1 + (x1[i]-avg_1)**2
    var2 = var2 + (x2[i]-avg_2)**2
    var3 = var3 + (x3[i]-avg_3)**2
    var4 = var4 + (x4[i]-avg_4)**2
  var1 = sqrt(var1/(len(x1)-1))
  var1 = var1/sqrt(len(x1))
  #print 'error_0 = ', var1
  var2 = sqrt(var2/(len(x2)-1))
  var2 = var2/sqrt(len(x2))
  #print 'error_1 = ', var2
  var3 = sqrt(var3/(len(x3)-1))
  var3 = var3/sqrt(len(x3))
  #print 'error_2 = ', var3
  var4 = sqrt(var4/(len(x4)-1))
  var4 = var4/sqrt(len(x4))
  #print 'error_2 = ', var2


  a2 = avg_2/avg_1
  #print 'N1/N0 = ', a1
  if a2!=0:
    #error_N1_N0 = abs(a1) * sqrt( (var0/avg_0)**2 + (var1/avg_1)**2 )
    error_N2_N1 = abs(a2) * (var1/abs(avg_1) + var2/abs(avg_2) )
    #print 'error_N1/N0 = ', error_N1_N0 
  a3=avg_3/avg_1
  #print 'N_2/N_0 = ', a2
  if a3!=0:
    #error_N2_N0 = abs(a2) * sqrt( (var0/avg_0)**2 + (var2/avg_2)**2 )
    error_N3_N1 = abs(a3) * ( var1/abs(avg_1) + var3/abs(avg_3) )
    #print 'error_N2/N0 = ', error_N2_N0
  a4=avg_4/avg_1
  #print 'N_2/N_0 = ', a2
  if a4!=0:
    #error_N2_N0 = abs(a2) * sqrt( (var0/avg_0)**2 + (var2/avg_2)**2 )
    error_N4_N1 = abs(a4) * ( var1/abs(avg_1) + var4/abs(avg_4) )
    #print 'error_N2/N0 = ', error_N2_N0
  Result = a2+a3+a4
  #print 'Result = ', Result
  if avg_2!=0:
    #error = sqrt (error_N1_N0**2 + error_N2_N0**2)
    error = error_N2_N1 + error_N3_N1 + error_N4_N1
  else:
    error = error_N3_N1 + error_N4_N1
  #print 'error = ', error
  return a2, error_N2_N1, a3, error_N3_N1, a4, error_N4_N1, a2+a3+a4, error
#End
###################################

###################################
######Read external variables######
###################################
#Begin
data = np.loadtxt('ext_vars.dat')  # reading external variables from a file
count = len(open('ext_vars.dat').readlines())
if count==1:
  EXT_VARS = [None]*count
  EXT_VARS[0]=data
if count>1:
  EXT_VARS = [None]*count  
  for i in range (0, count):
    EXT_VARS[i] = data[i,:]  
#print EXT_VARS
t=EXT_VARS[1][0]
U=EXT_VARS[1][1]
beta=EXT_VARS[1][2]
mu=EXT_VARS[1][3]
nu_ext=(2*EXT_VARS[1][4]+1)*pi/beta
px_ext=EXT_VARS[1][5]
py_ext=EXT_VARS[1][6]
lat_const=EXT_VARS[1][7]
print
print "Evaluation for the following parameters:"
print 't = ', t
print 'U = ', U
print 'beta = ', beta
print 'mu = ', mu
print 'nu_ext = ', nu_ext
print 'px_ext = ', px_ext
print 'py_ext = ', py_ext
print 'a = ', lat_const
#End
###################################
highest_order=3
flag_imag=0   # 1 for imaginary part, and 0 for real part computation


old_order=1  # Start from the first order term
old_number=1  # For now it is always 1
#for i in range (0,100):
  #print new_state(old_order,old_number)

N=500000# Length of Markov chain
T=100# Monte Carlo iteration
#print 'a = ', a
#print 'b = ', b
print 'N = ', N
print 'T = ', T

#measure(old_order,old_number)


chains = update(old_order,old_number)  # Measurement
output = measurement(chains[0],chains[1],chains[2])  # Output

if flag_imag==0:
  print '****'
  print 'Results for the real part:'
else:
  print '****'
  print 'Results for the imaginary part:'

print 'a2 = ', output[0], '+-', output[1]
print 'a3 = ', output[2], '+-', output[3]
print 'a4 = ', output[4], '+-', output[5]
print 'Result = ', output[6], '+-', output[7]

flag_imag=1   # 1 for imaginary part, and 0 for real part computation


old_order=1  # Start from the first order term
old_number=1  # For now it is always 1
#for i in range (0,100):
  #print new_state(old_order,old_number)

'''
N=10000000# Length of Markov chain
T=20# Monte Carlo iteration
#print 'a = ', a
#print 'b = ', b
print 'N = ', N
print 'T = ', T

#measure(old_order,old_number)


chains = update(old_order,old_number)  # Measurement
output = measurement(chains[0],chains[1],chains[2])  # Output

if flag_imag==0:
  print '****'
  print 'Results for the real part:'
else:
  print '****'
  print 'Results for the imaginary part:'

print 'a2 = ', output[0], '+-', output[1]
print 'a3 = ', output[2], '+-', output[3]
print 'Result = ', output[4], '+-', output[5]
'''
'''
flag_imag=0
T=20
N=50000
m=3
print 'T = ', T
print 'N = ', N
avg_re=0
avg_im=0
ch_re=[]
ch_im=[]
for i_T in range(0,T):
  for i_N in range(0,N):
    result=Integrand(m,1)
    re=result
    #im=result[1]
    avg_re=avg_re+re
    #avg_im=avg_im+im
  print (2*pi)**(2*m)*avg_re/N
  ch_re.append(avg_re/N)
  ch_im.append(avg_im/N)
  #print 're = ', avg_re/N
  #print 'im = ', avg_im/N
  avg_re=0 
  #avg_im=0

re=0
im=0
for i_T in range(0,T):
  re=re+ch_re[i_T]
  #im=im+ch_im[i_T]
re_fin = (2*pi)**(2*m)*re/T
im_fin = (2*pi)**(2*m)*im/T
var_re=0
var_im=0
for i_T in range(0,T):
  var_re=var_re+(re_fin-(2*pi)**(2*m)*ch_re[i_T])**2
  var_im=var_im+(im_fin-(2*pi)**(2*m)*ch_im[i_T])**2
var_re = sqrt(var_re/(T-1.0))
var_im = sqrt(var_im/(T-1.0))
er_re = var_re/sqrt(T)
er_im = var_im/sqrt(T) 
print 'Real Part = ', re_fin, '+-', er_re
print 'Imaginary Part = ', im_fin, '+-', er_im 
  

#for i in range (0,50):
  #print new_state(2,1)

'''





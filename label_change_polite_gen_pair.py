'''
Polite AMI+diagMC for the
self-energy up to the second order

Last updates: MAY-8-2019
              MAY-13-2019
'''

###################################
##########Import libraries#########
###################################
#Begin
from random import *
import pickle
from math import *
import numpy as np
from copy import deepcopy
from igraph import *
from Label_self import *
#from Symbolic_multi_AMI_new import *
#End
###################################

###################################
#######Frequency sign chanage######
###################################
#Begin
def freq_sign_chng(g_freq,p): 
  val=deepcopy(g_freq)
  #print 'val_init = ', val
  for i_g in range (0,len(val)):
    val[i_g][p-1]=-val[i_g][p-1]
  #print 'init=',g_freq
  #print 'final=',val
  return val 
#End
###################################

###################################
#######Momentum sign chanage#######
###################################
#Begin
def mnta_sign_chng(g_mnta,p): 
  val=deepcopy(g_mnta)
  #print 'val_init = ', val
  for i_g in range (0,len(val)):
    val[i_g][p-1]=-val[i_g][p-1]
  #print 'init=',g_mnta
  #print 'final=',val
  return val 
#End
###################################

###################################
#######Momentum shift by pi########
###################################
#Begin
def mnta_shift(g_freq,signs,p):
  val=signs
  for i_g in range (0,len(val)):
    if g_freq[i_g][p-1]!=0:
      val[i_g]=-val[i_g]
  #print 'new_signs=', val
  return val
#End
###################################

###################################
#########Frequency swap############
###################################
#Begin
def freq_swap(g_freq,p1,p2):
  val_freq=deepcopy(g_freq)
  for i_g in range(0,len(val_freq)):
    temp=g_freq[i_g][p1-1]
    g_freq[i_g][p1-1]=g_freq[i_g][p2-1]
    g_freq[i_g][p2-1]=temp
  return g_freq
#End
###################################

###################################
#########M#omentum swap############
###################################
#Begin
def mnta_swap(g_mnta,p1,p2):
  val_mnta=deepcopy(g_mnta)
  for i_g in range(0,len(val_mnta)):
    temp=g_mnta[i_g][p1-1]
    g_mnta[i_g][p1-1]=g_mnta[i_g][p2-1]
    g_mnta[i_g][p2-1]=temp
  return g_mnta
#End
###################################


############################
###two arrays are equal?####
############################
#Begin
def two_arr_eq(r1,r2):
  val = True
  for i in range (0, len(r1)):
    if r1[i] != r2[i]:
      val = False
  return val
#End
############################

########################################################
##Comarison of two arrays considering the overall sign##
########################################################
#Begin      
def arr_comp_new(r1, r2):
  ovl_sgn = 1
  r2_new = [None]*len(r2)
  val = False
  if two_arr_eq(r1,r2):
    val = True
    ovl_sgn = 1
  for i in range(0,len(r2)):
    r2_new[i] = -r2[i] 
  if two_arr_eq(r1,r2_new):
    val = True
    ovl_sgn = ovl_sgn*(-1)
  
  return val,ovl_sgn 
#End
########################################################

###################################
###########Green comp##############
###################################
#Begin
def g_comp(g1,g2):
  #flag_g_comp=True
  ovl_sgn=1
  out=[]
  g2_copy=deepcopy(g2)
  rel_sgn=[]
  count=0
  #print arr_comp_new(g1[0],g2[1])
  
  for i in range (0,len(g1)):
    #print 'g1 = ', g1[i]
    
      for j in range (0,len(g2_copy)):
       
        if arr_comp_new(g1[i],g2_copy[j])[0]:
          count=count+1
          out.append(g2_copy[j])
          flag_g_comp=True
          ovl_sgn=ovl_sgn*arr_comp_new(g1[i],g2_copy[j])[1]
          
          rel_sgn.append(arr_comp_new(g1[i],g2_copy[j])[1])
          g2_copy[j]=[0]*len(g2_copy[0])
          break
        #else: 
          #flag_g_comp=False
  #if len(out)==len(g2_copy):
    #val=True
  #print 'count = ', count
  if count==len(g2_copy):
    val=True
  else:
    val=False
  return val,ovl_sgn,rel_sgn
#End
###################################

###################################
######Is transformation valid?#####
###################################
def valid(g_freq,p):
  count_p=0
  for i in range(0,len(g_freq)):
    if g_freq[i][p-1]!=0:
      count_p=count_p+1
  #print count_p
  return count_p
#End
###################################  


###################################
#########Decision func#############
###################################
#Begin
def dec_polite(int_freq,int_mnta,int_signs,NEW_freq,NEW_mnta,NEW_signs,int_num_loop,NEW_num_loop):
  polite=1
  flag=True
  overal_sign_freq=1
  comp=g_comp(NEW_freq, int_freq)
    #print comp
  if comp[0]==False:
    flag=False
  overal_sign_freq=comp[1]*overal_sign_freq*((-1)**(int_num_loop-NEW_num_loop))
  #print 'flag_freq= ', flag
  #print overal_sign_freq

  if flag:
    comp=g_comp(NEW_mnta, int_mnta)
    #print comp
    if comp==False:
      flag=False
  #print 'flag_freq_mnta= ', flag
  check_freq=deepcopy(NEW_freq)
  if flag:
    #comp=arr_comp_new(int_signs,NEW_signs)
    sign_final=1
    for i in range (0,len(NEW_signs)):
      sign_final=sign_final*NEW_signs[i]
      for j in range (0, len(NEW_freq[i])):
        check_freq[i][j]=NEW_signs[i]*check_freq[i][j]
    #print 'sign_final = ', sign_final
    #print 'overal_sign = ', overal_sign_freq
    #print 'check_freq = ', check_freq
    comp=g_comp(check_freq,int_freq)
    #print 'comp = ', comp
    if sign_final==overal_sign_freq==-1 and comp[0] and comp[1]==1 and comp[2]==[1]*len(comp[2]) :
      #print 'Impolite'
      polite=0
  return polite
#End
###################################

impolites=[]
diags=[1, 2, 4, 6, 10, 17, 18, 20, 21, 23, 24, 27, 30]
M=4
repeat=10000
diag_num=35
print diags
print len(diags)
for j in range (0,100):
  for NUMBER in diags:
    file_name='m_'+str(M)+'_num_'+str(NUMBER)+'.graphml'
    G=load(file_name)
        
    reset_g(G,M)
    label_abs_ran(G,M,20,20)
    ami_in = deepcopy(AMI_Input(G.es["Label"],M))
    #file_name_G_sym = 'G_sym_m_'+str(M)+'_num_'+str(NUMBER)
    #with open(file_name_G_sym,"rb") as fp:
      #g_sym=pickle.load(fp)
    file_name_f = 'f_m_'+str(M)+'_num_'+str(NUMBER)
    with open(file_name_f,"rb") as fp:
      int_f=pickle.load(fp)
    int_freq_labels = deepcopy(ami_in)
    energy_signs=[1]*len(int_freq_labels)
    int_mnta_labels=deepcopy(int_freq_labels)
    if NUMBER not in impolites:
      flag=True
    else:
      flag=False
  
    for NUMBER_comp in diags:
      if flag and NUMBER!=NUMBER_comp and NUMBER not in impolites:
        #print 
        print NUMBER, NUMBER_comp
        file_name='m_'+str(M)+'_num_'+str(NUMBER_comp)+'.graphml'
        Gp=load(file_name)
        #print G
        #print G.es["F_or_B"]
        #print G.es["INT_or_EXT"]
        #loops=all_loops(G,M)
        #print 'loops = ', loops
        #F=len(loops)-1
        #print 'F = ', F
        #ext_max=0
        for z in range (0,100):
          reset_g(Gp,M)
          label_abs_ran(Gp,M,20,20)
          ami_in_p = deepcopy(AMI_Input(Gp.es["Label"],M))
          #print ami_in
          #file_name_G_sym_p = 'G_sym_m_'+str(M)+'_num_'+str(NUMBER_comp)
          file_name_f_p = 'f_m_'+str(M)+'_num_'+str(NUMBER_comp)
          with open(file_name_f_p,"rb") as fp:
            int_f_p=pickle.load(fp)
          int_freq_labels_p = deepcopy(ami_in_p)
          energy_signs_p=[1]*len(int_freq_labels_p)
          int_mnta_labels_p=deepcopy(int_freq_labels_p)

          new_signs=deepcopy(energy_signs)
          new_freq=deepcopy(int_freq_labels)
          new_mnta=deepcopy(int_mnta_labels)

          for i in range(0,repeat):
            r=random()
            if r<0.5:
              coin1=randint(1,M)
              coin2=randint(1,M)
              if coin1!=coin2 and valid(new_freq,coin1)!=1 and valid(new_freq,coin2)!=1:
                new_freq=deepcopy(freq_swap(new_freq,coin1,coin2))
                new_mnta=deepcopy(mnta_swap(new_mnta,coin1,coin2))  
            dec=dec_polite(int_freq_labels_p,int_mnta_labels_p,energy_signs_p,new_freq,new_mnta,new_signs,int_f,int_f_p)
            if dec==0:
              #print 'new_freq = ', new_freq
              #print 'new_mnta = ', new_mnta
              #print 'new_signs = ', new_signs
              
              if NUMBER not in impolites:
                impolites.append(NUMBER)
                impolites.append(NUMBER_comp)
                print impolites
                print i
              flag=False
              #break
            r=random()
            if r<0.5:
              coin=randint(1,M)
              if valid(new_freq,coin)!=1:
                new_freq=deepcopy(freq_sign_chng(new_freq,coin)) 
                new_mnta=deepcopy(mnta_sign_chng(new_mnta,coin)) 
            dec=dec_polite(int_freq_labels_p,int_mnta_labels_p,energy_signs_p,new_freq,new_mnta,new_signs,int_f,int_f_p)
            if dec==0:
              #print 'new_freq = ', new_freq
              #print 'new_mnta = ', new_mnta
              #print 'new_signs = ', new_signs
              
              if NUMBER not in impolites:
                impolites.append(NUMBER)
                impolites.append(NUMBER_comp)
                print impolites
                print i
              flag=False
              #break
  
            if r>=0.5:
              coin=randint(1,M)
              if valid(new_freq,coin)!=1:
                new_signs=mnta_shift(new_freq,new_signs,coin)
            dec=dec_polite(int_freq_labels_p,int_mnta_labels_p,energy_signs_p,new_freq,new_mnta,new_signs,int_f,int_f_p)
            #print new_freq
            #print new_mnta
            #print new_signs
            if dec==0:
              #print 'new_freq = ', new_freq
              #print 'new_mnta = ', new_mnta
              #print 'new_signs = ', new_signs
              
              if NUMBER not in impolites:
                impolites.append(NUMBER)
                impolites.append(NUMBER_comp)
                print impolites
                print i
              flag=False
              #break  
print 'impolites = ', impolites
polites=[]

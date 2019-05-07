'''
Polite AMI+diagMC for the
self-energy up to the second order

Last updates: MAY-7-2019
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
#########Decision func#############
###################################
#Begin
def dec_polite(int_freq,int_mnta,int_signs,NEW_freq,NEW_mnta,NEW_signs):
  polite=1
  flag=True
  overal_sign_freq=1
  for i in range (0,len(NEW_freq)):
    comp=arr_comp_new(NEW_freq[i], int_freq[i])
    #print comp
    if comp[0]==False:
      flag=False
      break
    overal_sign_freq=comp[1]*overal_sign_freq
  #print 'flag_freq= ', flag
  #print overal_sign_freq

  if flag:
    for i in range (0,len(NEW_mnta)):
      comp=arr_comp_new(NEW_mnta[i], int_mnta[i])
      #print comp
      if comp[0]==False:
        flag=False
        break
  #print 'flag_freq_mnta= ', flag

  if flag:
    comp=arr_comp_new(int_signs,NEW_signs)
    #print comp
    if comp[1]==overal_sign_freq==-1:
      #print 'Impolite'
      polite=0
  return polite
#End
###################################
    

M=3
NUMBER=6
file_name_G_sym = 'G_sym_m_'+str(M)+'_num_'+str(NUMBER)
with open(file_name_G_sym,"rb") as fp:
  g_sym=pickle.load(fp)
int_freq_labels = g_sym[0]
#print 'int_freq_labels = ',int_freq_labels
#init_energy_labels = g_sym[1]
#print 'init_energy_labels = ', init_energy_labels
energy_signs=[1]*len(int_freq_labels)
#print 'initial_energy_signs=', energy_signs
#print

#

int_mnta_labels=deepcopy(int_freq_labels)

#print int_mnta_labels

#freq_sign_chng(int_freq_labels,2)
#print
#mnta_sign_chng(int_mnta_labels,1)
#print
#mnta_shift(int_freq_labels,energy_signs,3)


'''
###################################
############Simple test############
###################################
#Begin
new_freq=deepcopy(int_freq_labels)
print
print 'initial_freq = ', new_freq
for i in range (1,4):
  new_freq=deepcopy(freq_sign_chng(new_freq,i))
print 'final_freq = ', new_freq 
print
new_mnta=deepcopy(int_mnta_labels)
print 'initial_mnta = ', new_mnta
for i in range (1,4):
  new_mnta=deepcopy(mnta_sign_chng(new_mnta,i))
print 'final_mnta = ', new_mnta 
print
new_signs=deepcopy(energy_signs)
for i in range (0,len(energy_signs)):
  new_signs=mnta_shift(int_freq_labels,new_signs,i)
print 'new_sign=', new_signs
print

#print dec_polite(int_freq_labels,int_mnta_labels,energy_signs,new_freq,new_mnta,new_signs)

#flag=True
#overal_sign_freq=1
#for i in range (0,len(new_freq)):
  #comp=arr_comp_new(new_freq[i], int_freq_labels[i])
  #print comp
  #if comp[0]==False:
    #flag=False
    #break
  #overal_sign_freq=comp[1]*overal_sign_freq
#print 'flag_freq= ', flag
#print overal_sign_freq

#if flag:
  #for i in range (0,len(new_mnta)):
    #comp=arr_comp_new(new_mnta[i], int_mnta_labels[i])
    #print comp
    #if comp[0]==False:
      #flag=False
      #break
#print 'flag_freq_mnta= ', flag


#if flag:
  #comp=arr_comp_new(energy_signs,new_sign)
  #print comp
  #if comp[1]==overal_sign_freq==-1:
    #print 'Impolite'
#End
###################################
'''

###################################
############Another test###########
###################################
#Begin
flag_loop=False
new_signs=deepcopy(energy_signs)
new_freq=deepcopy(int_freq_labels)
new_mnta=deepcopy(int_mnta_labels)
print int_freq_labels
for i in range(0,10000):
  r=random()
  if r<0.5:
    coin=randint(1,M)
    new_freq=deepcopy(freq_sign_chng(new_freq,coin)) 
    new_mnta=deepcopy(mnta_sign_chng(new_mnta,coin)) 
  dec=dec_polite(int_freq_labels,int_mnta_labels,energy_signs,new_freq,new_mnta,new_signs)
  if dec==0:
    print 'new_freq = ', new_freq
    print 'new_signs = ', new_signs
    print i
    break
  if r>=0.5:
    coin=randint(1,M)
    new_signs=mnta_shift(int_freq_labels,new_signs,coin)
  dec=dec_polite(int_freq_labels,int_mnta_labels,energy_signs,new_freq,new_mnta,new_signs)
  if dec==0:
    print 'new_freq = ', new_freq
    print 'new_signs = ', new_signs
    print i
    break  
  
print i






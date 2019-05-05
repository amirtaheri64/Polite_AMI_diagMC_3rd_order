'''
Symbolic AMI algorithm considering multiple poles

Last updates March-25-2019
             March-26-2019
             March-27-2019
             March-31-2019

'''
from __future__ import division
from math import *
import math
import scipy.integrate 
import numpy
import time
import numpy as np
import pylab as plt
from numpy.linalg import inv
from numpy import linalg as LA
from copy import deepcopy
start_time = time.time()

##################################
###Definition of Fermi Function###
##################################
#Begin
def fermi(e, m, B):  # Definition of fermi function
  val = 1.0/(math.exp(B*(e-m))+1.0)
  return val
#End
##################################

#################################################
##Dispersion relation of 2D tight-binding model##
#################################################
#Begin
def E( px, py, t, a): # 2D tight-binding dispersion
  val = 2*t*(cos(a*px) + cos(a*py)) 
  return val
#End
##################################################

##################################################
##############Factorial function##################
##################################################
#Begin
def fact(n):
  if n==0 or n==1 :
    val = 1
  else:
    val = n*fact(n-1)
  return val
#End
##################################################

#################################
#####Symbolic Representation#####
#################################
#Begin
def sym_rep_g(g):  # symbolic representation of multiplication of Green's functions
  num_l=len(g)
  e=[[0]*(num_l)for i in range (0,num_l)]
  for j in range (0,num_l):
    e[j][j]=1
  val=[]
  print
  print 'ami_in = ', g 
  print
  
  print
  for i in range(0,num_l):
    val.append(e[i])
    val.append(g[i])
  #print 'g=',g
  #print 'e=',e
  #print 'val=',val
  for i in range (0,len(g)):
    for j in range (i+1,len(g)):
      if arr_comp_new(g[i], g[j])[0]:
        #print arr_comp_new(g[i], g[j])
        e[j][j]=0
        e[j][i]=e[i][i]
  print 'e = ', e    
  return g,e,val  # g:frequency part, e: momentum part, val: total representation
#End
#################################

##################################################
######Derivation of single Green's function#######
##################################################
#Begin
def der(g_freq, g_mnta, p):  # g_freq: frequecy part of the input Green's function, g_mnta: momentum part of the green's function, p: the index of nu_p as derivation variable
  der_g_freq = [[None]*len(g_freq) for i in range(0, 2)]  # derivative of frequency part
  der_g_mnta = [[None]*len(g_mnta) for i in range(0, 2)]  # derivative of momentum part
  for i in range (0, len(g_freq)):
    if g_freq[p]==0:  # If pth frequency is absent 
      #s = 0
      for j in range (0,2):
        for k in range (0, len(g_freq)):
          der_g_freq[j][k] =0 #s*g_freq[k]
        for k in range (0, len(g_mnta)):
          der_g_mnta[j][k] =0 #s*g_freq[k]
    else:
      s = -1*g_freq[p] # s is the sign attached to the result of derivation with respect to pth frequency
      for j in range (0,2):
        for k in range (0, len(g_freq)):
          der_g_freq[j][k] = g_freq[k]/(s**j)
        for k in range (0, len(g_mnta)):
          der_g_mnta[j][k] = g_mnta[k]/(s**j)
  return der_g_freq, der_g_mnta
#End
####################################################

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

########################################################
###Derivative of multiplication of Green's functions####
########################################################
#Begin
def der_mul(g_mul_freq, g_mul_mnta, p):
  val_freq=[]
  val_mnta=[]
  if len(g_mul_freq)>0:
    num_f = len(g_mul_freq[0])  
    num_m = len(g_mul_mnta[0])
    l=len(g_mul_freq)
    val_freq = np.zeros ((len(g_mul_freq), len(g_mul_freq)+1, num_f))
    val_mnta = np.zeros ((len(g_mul_mnta), len(g_mul_mnta)+1, num_m))
    k = 0
  
    for i in range (0, l):
      for j in range (0, l+1):
        if i==k:
          a=der(g_mul_freq[k],g_mul_mnta[k],p)
          val_freq[k][k] = a[0][0] 
          val_freq[k][k+1] = a[0][1]
          val_mnta[k][k] = a[1][0] 
          val_mnta[k][k+1] = a[1][1]
          for h in range(k+2,l+1):
            val_freq[k][h] = g_mul_freq[h-1]
            val_mnta[k][h] = g_mul_mnta[h-1]
          for h in range(0,k):
            val_freq[k][h] = g_mul_freq[h] 
            val_mnta[k][h] = g_mul_mnta[h]  
          k = k + 1
  return val_freq, val_mnta  
#End
########################################################

##############################################
##To determine the multiplicity of each pole##
##############################################
#Begin
def decision(r_freq,r_mnta):
  a = [] 
  num = []
  result_freq = []
  result_mnta = []
  ovl_sgn = 1
  for i in range(0,len(r_freq)):
    a.append(0)
    num.append(0)
    result_freq.append(0)
    result_mnta.append(0)
  l = 0 
  for i in range (0,len(r_freq)):
    k = 0
    
    if a[i] == 0:
      a[i] = 1
      
      result_freq[l] = r_freq[i]
      result_mnta[l] = r_mnta[i]  
      for j in range (i+1, len(r_freq)):
        comp_freq = arr_comp_new(r_freq[i],r_freq[j])
        comp_mnta = arr_comp_new(r_mnta[i],r_mnta[j])
        if comp_freq[0] and comp_mnta[0] and comp_freq[1]==comp_mnta[1]:
          ovl_sgn = ovl_sgn*comp_freq[1]
          a[j] = 1
          k = k+1
      num[l] = k+1
      #print result[l], num[l], l
      l = l + 1
     
  return result_freq, result_mnta, num, ovl_sgn
#End
##############################################

##############################################
###To do one summation over p-th frequency####
##############################################
#Begin
def one_sum(G_fr,G_mo,p): # G_fr: frequency part, G_mo: momentum part, p: frequency 0,1,...,n-1
  TEMP1_freq = []  # Auxiliary variable
  TEMP1_mnta = []  # Auxiliary variable
  #print TEMP1
  final_index = 0
  dec = decision(G_fr,G_mo)
  #print 'decision = ', dec
  #print dec[1][2]
  #print len(dec[1])
  l = 0
  for i in range (0, len(dec[2])):
    if dec[2][i]!=0:
      #print i, dec[1][i], dec[0][i]
      l = l + 1
  #print l
  dec_new = [[None]*l for i in range(0,3)]
  #print 'dec_new = ', dec_new
  for i in range(0,l):
    dec_new[0][i] = dec[0][i]
    dec_new[1][i] = dec[1][i] 
    dec_new[2][i] = dec[2][i]
  overal_sign = dec[3]
  #print 'dec_new = ', dec_new
  #print 'overal_sign = ', overal_sign
  g_freq = []  # To store freequency part
  g_mnta = []  # To store momentum part
  G_freq = deepcopy(G_fr)  # Frequency part
  G_mnta = deepcopy(G_mo)  # Momentum part
  #val_freq = [[1]*(len(G_freq[0])) for i in range(len(G_freq))]  # To store freequency part
  #val_mnta = [[1]*(len(G_mnta[0])) for i in range(len(G_mnta))]  # To store momentum part
  alpha_p = []
  z = []
  sign = []
  mul = []  # To store multiplicity
  result_freq = []
  result_mnta = []
  root_freq = []  # To store poles' frequency
  root_mnta = []  # To store poles' mnta
  new_z=[]
  for i in range(0, len(dec_new[0])):
    g_freq.append(0)
    g_mnta.append(0)
    #z.append(None)
    #root_freq.append(None)
    alpha_p.append(None)
    mul.append(0)
    result_freq.append(0)
    result_mnta.append(0)
    #sign.append(1)
  for i in range(0, len(dec_new[0])):  
    g_freq[i] = deepcopy(dec_new[0][i])
    g_mnta[i] = deepcopy(dec_new[1][i])
    mul[i] = deepcopy(dec_new[2][i])
    alpha_p[i] = dec_new[0][i][p]
  #for i in range(0, len(G_fr)):
    #sign.append(None)
    #z.append(None)
    #root_freq.append(None)
    #root_mnta.append(None)
  root_mnta=[[None]*(len(G_mo[0]))for i in range (0,len(G_mo))]  # Initialize root_mnta
  root_freq=[[None]*(len(G_fr[0]))for i in range (0,len(G_fr))]  # Initialize root_freq
  #print 'g_freq = ', g_freq
  #print 'g_mnta = ', g_mnta
  #print 'mul = ', mul
  #print 'alpha = ', alpha_p  
  new_index = 0
  index = 0
  for j in range (0,len(g_freq)):
    #for i in range(0, len(g)):  
      #g[i] = deepcopy(b[i])
    for i in range(0, len(g_freq)):
      alpha_p[i] = g_freq[i][p]
    #print 'g_freq = ', g_freq[j]
    #print 'g_mnta = ', g_mnta[j]
  
    #print 'sign[j] = ', sign[j]
    #print 'root_freq[j] = ', root_freq[j]
    #print 'root_mnta[j] = ', root_mnta[j]
    #print 'z[j] = ', z[j]
    if mul[j] == 1:
      #print 'mul[j] = ', mul[j]
      #print '*********j = ', j
      #print 'g_freq[j] = ', g_freq[j]
      #print 'g_mnta[j] = ', g_mnta[j]
      #print 'g[j][p] = ', g[j][p]
      if g_freq[j][p] != 0:
        sign.append(g_freq[j][p])
        #sign[final_index] = g_freq[j][p]
        #root_freq[final_index] = -alpha_p[j]*g_freq[j]
        for i in range(0,len(g_freq[j])):
          root_freq[final_index][i] = -alpha_p[j]*g_freq[j][i]  
        for i in range(0,len(g_mnta[j])): 
          root_mnta[final_index][i]=-alpha_p[j]*g_mnta[j][i]  # Find the momnetum of the pole
        z.append(-alpha_p[j]*g_freq[j][0])
        #z[final_index] = -alpha_p[j]*g_freq[j][0]
        root_freq[final_index][p] = 0
      #print 'sign[final_index] = ', sign[final_index]
      #print 'root_freq[final_index] = ', root_freq[final_index]
      #print 'root_mnta[final_index] = ', root_mnta[final_index]
      if g_freq[j][p]!=0:
        for k in range (0, len(g_freq)):
          if k!=j:
            index = k   
        
            g_freq[index][p] = 0
            g_freq[j][p] = 0
        
            #if alpha_p[j]!=0:
              #z[j] = -alpha_p[j]*g[j][0]
              #root[j] = -alpha_p[j]*g[j] 
            #if alpha_p[j]!=0:  # If pth frequency is present
              #z[j] = -alpha_p[j]*g_freq[j][0]  # Get the zeroth element as flag
              #root_freq[j]=-alpha_p[j]*g_freq[j] # Find frequency of the pole
              #for i in range(0,len(g_mnta[j])): 
                #root_mnta[j][i]=-alpha_p[j]*g_mnta[j][i]  # Find the momnetum of the pole
        
            for i in range(0, len(g_freq[j])):
              g_freq[index][i] = g_freq[index][i]-alpha_p[index]*alpha_p[j]*g_freq[j][i] # Find the residue: frequency part
            for i in range(0, len(g_mnta[j])):
              g_mnta[index][i] = g_mnta[index][i]-alpha_p[index]*alpha_p[j]*g_mnta[j][i] # Find the residue: momentum part
        
            #if alpha_p[j] == +1:
              #sign[j] = +1
            #if alpha_p[j] == -1:
              #sign[j] = -1
        tot_len = 0 
        for i in range (0, len(mul)):
          tot_len = tot_len + mul[i]
        #print tot_len   
        val_freq = [ [None]*len(g_freq[0]) for i in range (0, tot_len-1) ]
        val_mnta = [ [None]*len(g_mnta[0]) for i in range (0, tot_len-1) ]
        #print val
        val_index = 0
        for i in range (0, len(mul)):
          for h in range (0, mul[i]):
            if i!= j:
              val_freq[val_index] = g_freq[i]
              val_mnta[val_index] = g_mnta[i]
              val_index = val_index + 1
        #print 'val_mul[j] = ', mul[j], val
        for i in range(0, len(dec_new[0])):  
          g_freq[i] = deepcopy(dec_new[0][i])
          g_mnta[i] = deepcopy(dec_new[1][i])
          mul[i] = deepcopy(dec_new[2][i])
        TEMP1_freq.append(val_freq)
        TEMP1_mnta.append(val_mnta) 
        #TEMP1_freq[final_index] = val_freq
        #TEMP1_mnta[final_index] = val_mnta 
        #print val
        final_index = final_index + 1 
      #for i in range(0, len(dec_new[0])):
        #alpha_p[i] = dec_new[0][i][p]
        #print g
    count_new = 0
    for i in range(0, len(dec_new[0])):
      g_freq[i] = deepcopy(dec_new[0][i])
      g_mnta[i] = deepcopy(dec_new[1][i])
      mul[i] = deepcopy(dec_new[2][i])
    if mul[j]>1 and g_freq[j][p]!=0:
      #for i in range(0, len(dec_new[0])):  
        #g[i] = deepcopy(dec_new[0][i])
        #mul[i] = deepcopy(dec_new[1][i]) 
      #for i in range(0, len(dec_new[0])):
        #alpha_p[i] = dec_new[0][i][p]
      #print 'alpha = ', alpha_p
      #print 'g[j] = ', g[j]
      #print 'mul[j] = ', mul[j]
    
      for k in range(0,len(dec_new[0])):
        if k!=j:
          count_new = count_new + 1
      g_new_freq = [ [None]*len(G_fr[0]) for i in range (0, count_new) ]
      g_new_mnta = [ [None]*len(G_mo[0]) for i in range (0, count_new) ]
      mul_new = [ [None]*len(G_fr[0]) for i in range (0, count_new) ]
      count_new = 0
      #print 'g_new_freq = ', g_new_freq
      #print 'g_new_mnta = ', g_new_mnta
      #print 'mul_new = ', mul_new
      for k in range(0,len(dec_new[0])):
        if k!=j:
          g_new_freq[count_new] = dec_new[0][k]
          g_new_mnta[count_new] = dec_new[1][k]
          mul_new[count_new] = dec_new[2][k]
          count_new = count_new + 1  
      #print 'g_new_freq = ', g_new_freq
      #print 'g_new_mnta = ', g_new_mnta
      #print 'mul_new = ', mul_new
      #print 'g_new[j] = ', g_new[j]
      #print 'mul_new = ', mul_new
      tot_len = 0
      for i in range (0, len(mul_new)):
        tot_len = tot_len + mul_new[i]
      #print 'tot_len = ', tot_len  
      val_freq = [ [None]*len(G_fr[0]) for i in range (0, tot_len) ]
      val_mnta = [ [None]*len(G_mo[0]) for i in range (0, tot_len) ]
      #print val
      val_index = 0
      for i in range (0, len(mul_new)):
        for h in range (0, mul_new[i]):
          val_freq[val_index] = g_new_freq[i]
          val_mnta[val_index] = g_new_mnta[i]
          val_index = val_index + 1
      #print 'val_freq = ', val_freq
      #print 'val_mnta = ', val_mnta 
      if mul[j] == 2:
        a = der_mul(val_freq,val_mnta,p)
        if len(a)>0:
          TEMP_freq = a[0]
          TEMP_mnta = a[1]
          der_g_new_freq = []
          der_g_new_mnta = []
          for i in range (0,len(TEMP_freq)):
            #print '*****', TEMP_freq[i]
            flag_non_zero=True
            for check1 in range (0,len(TEMP_freq[i])):
              for check2 in range (0,len(TEMP_freq[i][check1])):
                if two_arr_eq(TEMP_freq[i][check1],[0]*len(TEMP_freq[i][check1])):
                  flag_non_zero=False
            if flag_non_zero:   
              der_g_new_freq.append(TEMP_freq[i])
              der_g_new_mnta.append(TEMP_mnta[i])
          #print 'der_g_new_freq = ', der_g_new_freq
          #print 'der_g_new_mnta = ', der_g_new_mnta 
      else:
        length = 0
        a = der_mul(val_freq,val_mnta,p)
        if len(a)>0:
          der_old_freq = a[0]
          der_old_mnta = a[1]
          for h in range (1,mul[j]-1):
            length = 0
            for u in range (0, len(der_old_freq)):
              a = der_mul(der_old_freq[u],der_old_mnta[u],p) 
              if len(a)>0: 
                TEMP_freq = a[0]
                TEMP_mnta = a[1]
                len_TEMP = len(TEMP_freq)
                length = length + len_TEMP
            der_g_new_freq = []
            der_g_new_mnta = []
            len_TEMP = 0
            new_counter = 0
            count_non_zero=0
            for u in range (0, len(der_old_freq)):
              a = der_mul(der_old_freq[u],der_old_mnta[u],p)
              if len(a)>0:
                TEMP_freq = a[0]
                TEMP_mnta = a[1]
                #print 'len(TEMP) = ', len(TEMP)
            
            
                for i in range (0,len(TEMP_freq)):
                  #print '*****', TEMP_freq[i]
                  flag_non_zero=True
                  for check1 in range (0,len(TEMP_freq[i])):
                    for check2 in range (0,len(TEMP_freq[i][check1])):
                      if two_arr_eq(TEMP_freq[i][check1],[0]*len(TEMP_freq[i][check1])):
                        flag_non_zero=False
                  if flag_non_zero:   
                    der_g_new_freq.append(TEMP_freq[i])
                    der_g_new_mnta.append(TEMP_mnta[i])
                #count_non_zero=count_non_zero+1
                #new_counter = new_counter+1
              #len_TEMP = len(TEMP_freq)
            #count_non_zero=0
            #new_counter=0
            der_old_freq = der_g_new_freq
            der_old_mnta = der_g_new_mnta
          #print 'count_non_zero = ', count_non_zero
          #print 'der_g_new_freq = ', der_g_new_freq
          #print 'der_g_new_mnta = ', der_g_new_mnta
          #new_index = 0 
          #for i in range(0,len(der_g_new_freq)):
            #result[new_index+i] = der_g_new[i]
            #new_index = new_index + 1
          #len_result = new_index
        
        #len_result=0
        result_freq=[]
        result_mnta=[]
        #for i in range (0,len(der_g_new_freq)):
          #if der_g_new_freq[i]!=None:
            #len_result=len_result+1
            #result_freq.append(der_g_new_freq[i])
            #result_mnta.append(der_g_new_mnta[i])
        #print 'len_result = ', len_result
        #print 'len_result =', len_result
        #print 'result_freq = ', result_freq
        #print 'result_mnta = ', result_mnta
        #new_index = 0
        #result_freq = [ [None]*len(G_fr[0]) for i in range (0, len_result) ]
        #result_mnta = [ [None]*len(G_mo[0]) for i in range (0, len_result) ]
        #print 'result_freq = ', result_freq
        #print 'result_mnta = ', result_mnta
        #for i in range(0,len_result):
          #if der_g_new_freq[i]!=None:
            #result_freq[i] = der_g_new_freq[i]
            #result_mnta[i] = der_g_new_mnta[i]
        #print 'result_freq = ', der_g_new_freq
        #print 'result_mnta = ', der_g_new_mnta
        
        #print 'g[j] = ', g[j]
      if len(der_g_new_freq)!=0:
        result_new_freq = np.zeros((len(der_g_new_freq), len(der_g_new_freq[0]), len(der_g_new_freq[0][0])))
        result_new_mnta = np.zeros((len(der_g_new_mnta), len(der_g_new_mnta[0]), len(der_g_new_mnta[0][0])))
      
        for i in range(0, len(der_g_new_freq)):
          for h in range(0, len(der_g_new_freq[i])):
            for r in range(0, len(der_g_new_freq[i][h])):
            #print '*****', r,  result[i][h][r]
            #print 'g[j][r] = ', r,  g[j][r]
            #if r!=p:
              result_new_freq[i][h][r] = der_g_new_freq[i][h][r]-der_g_new_freq[i][h][p]*g_freq[j][p]*g_freq[j][r]
          #g_mnta[index][i] = g_mnta[index][i]-alpha_p[index]*alpha_p[j]*g_mnta[j][i]
            for r in range(0, len(der_g_new_mnta[i][h])):
              result_new_mnta[i][h][r] = der_g_new_mnta[i][h][r]-der_g_new_freq[i][h][p]*g_freq[j][p]*g_mnta[j][r]
            
        #print 'result_new_freq = ', result_new_freq
        #print 'result_new_mnta = ', result_new_mnta
        len_result = len(der_g_new_freq)
        #print 'len_result = ', len_result
        for i in range (0, len(result_new_freq)):
          if g_freq[j][p] != 0:
            sign.append((g_freq[j][p]**mul[j])/fact(mul[j]-1))
            #sign[final_index] = g_freq[j][p]/fact(mul[j]-1)
            #root_freq[final_index] = -alpha_p[j]*g_freq[j]
            for count_freq in range(0,len(g_freq[j])):
              root_freq[final_index][count_freq] = -alpha_p[j]*g_freq[j][count_freq] 
            for count_mnta in range(0,len(g_mnta[j])): 
              root_mnta[final_index][count_mnta]=-alpha_p[j]*g_mnta[j][count_mnta]  # Find the momnetum of the pole
            z.append(-alpha_p[j]*g_freq[j][0])
            #z[final_index] = -alpha_p[j]*g_freq[j][0]
            root_freq[final_index][p] = 0
          TEMP1_freq.append(result_new_freq[i])
          TEMP1_mnta.append(result_new_mnta[i])
          #TEMP1_freq[final_index] = result_new_freq[i]
          #TEMP1_mnta[final_index] = result_new_mnta[i]
          final_index = final_index + 1 
  
  #print 'TEMP1_freq = ', TEMP1_freq
  #print 'TEMP1_mnta = ', TEMP1_mnta
  #print 'sign = ', sign
  #print 'z = ', z
  #print 'root_freq = ', root_freq
  #print 'root_mnta = ', root_mnta 
  #print 'final_index = ', final_index
  #print 'len_TEMP1_freq = ', len(TEMP1_freq)
  #fin_z = []
  #fin_sign = []
  fin_root_freq = []
  fin_root_mnta = []
  #fin_result_freq = []
  #fin_result_mnta = []   
  for i in range(0, final_index):
    #fin_z.append(z[i])
    fin_root_freq.append(root_freq[i])
    fin_root_mnta.append(root_mnta[i])
    #fin_sign.append(sign[i])
    #fin_result_freq.append(TEMP1_freq[i])  
    #fin_result_mnta.append(TEMP1_mnta[i])   
  #print 'fin_result_freq =', fin_result_freq
  #print 'fin_result_mnta =', fin_result_mnta
  #print 'fin_sign = ', fin_sign
  #print 'fin_z = ', fin_z
  #print 'fin_root_freq = ', fin_root_freq
  #print 'fin_root_mnta = ', fin_root_mnta
  if len(sign)!=0:
    measure=True
  else:
    measure=False
  return sign, fin_root_freq, fin_root_mnta, TEMP1_freq, TEMP1_mnta, overal_sign, measure
  #b = [ [0,1,0], [0,1,1], [1,1,0] ]
  #print one_sum_n_G(b, p)[2]
#End
##############################################

###################################
#########Initial energies##########
###################################
#Begin
# To evaluate the numerical values of dispersions
def Energy(g_freq_init,p,h,M): # p: momentum, h: hopping amplitude, M: chemical potential
  l = len(g_freq_init)
  e=[None] *l
  
  for i in range(0,l):
    val=[0]*2
    for j in range (0,len(g_freq_init[0])):
      #print 'p','j=', j, p[j]
      for r in range(0,len(p[0])):
        val[r] = p[j][r]*g_freq_init[i][j]+val[r]
    #print len(p[0])
    #print 'freq = ', g_freq_init[i]
    #print 'val[0] =', val[0]
    #print 'val[1] =', val[1]
    e[i]=-2*h*(cos(val[0])+cos(val[1]))-M 
    #print 'e[i] = ', e[i]
   
  return e
#End
###################################

###################################
#########Energy evaluation#########
###################################
#Begin
# To evaluate the numerical values of dispersions
def Energy_eval(e_array,g_freq_init,p,h,M):
  #print 'e_array = ', e_array
  e=Energy(g_freq_init,p,h,M)
  val=0
  for i in range(0,len(e_array)):
    val = val - e_array[i]*e[i]
    #print 'e_array[i] = ', e_array
    #print 'e[i] = ', e[i]
  #print 'e = ', e 
  #print g_freq_init 
  return val
#End
###################################

###################################
###########Poles merging###########
###################################
#Begin
#To merge frequency and momenta poles
def Poles(p_f, p_m, g_freq_init,p,h,M): # p_f: frequency array, p_m: momenta array of poles
 
  num_len = len(g_freq_init[0])+1
  p_f_copy = deepcopy(p_f)
  term = [None]*num_len
  for i in range (0,len(p_f)):
    for j in range (0,len(p_f[i])):
      if p_f_copy[i][j]!=None:
        #print 'p_f = ', p_f[i][j]
        #print 'p_m = ', p_m[i][j]
        for k in range (0,len(p_f[i][j])):
          #print 'p_f[i][j][k] = ', p_f[i][j][k]  
          #print 'p_m[i][j][k] = ', p_m[i][j][k]
          term[0] = Energy_eval(p_m[i][j][k],g_freq_init,p,h,M)
          #print 'term[0] = ', term[0] 
          for r in range (1, num_len):
            term[r] = p_f_copy[i][j][k][r-1]
          #print 'term = ', term
          p_f_copy[i][j][k] = deepcopy(term)
  return p_f_copy
#End
###################################

###################################
#########fermi/bose function#######
###################################
# Definition of f Function
#Begin
def f(z,B):
  lz = len(z)
  index = 0
  nu = [0]*lz
  for i in range (1, lz):
    if (z[i]!=0):
      index = index + 1
  sigma = (-1)**index
  #print index
  reg=0
  #print 'sigma = ', sigma
  #print 'z = ', z[0]
  if sigma*exp(B*z[0]-reg)+1.0 == 0:
    reg = 0.001
    print '????'
  val = 1.0/(sigma*exp(B*z[0]-reg)+1.0)
  return val
#End
###################################  

###################################
############dot_num################
###################################
#Begin
# Definition of dot operation between two array of numbers
def dot_num(a, b):
  la = len(a)
  lb = len(b)
  if la!=lb:
    print 'BEEP! Two arrays mush have the same length!'
  else:
    l = la
    val = [None]*l
    for i in range(0, l):
      val[i] = a[i]*b[i]
    #print 'A dot B = ', a, '.', b, '= ', val
    #return 'Done!' 
  return val
#End
###################################

###################################
############dot_arr################
###################################
#Begin
# Definition of dot operation between array of numbers and array of arrays of numbers
def dot_arr(c, d):
  lc = len(c)
  ld = len(d)
  count = 0
  if lc!=ld:
    print 'BEEP! Two arrays must have the same length!'
  if lc==ld:
    l = lc
    result = [None]*l
    for i in range(0,l):
      
      if len(c[i]) == len (d[i]):
        result[count] = dot_num(c[i],d[i])
        count = count + 1
      else:
        print 'BEEP! Two arrays elements must have the same length!'
        return None
  return result   
#End
###################################

###################################
##############cross################
###################################
#Begin
# Definition of cross operation  
def cross(a, c):
  la = len(a)
  lc = len(c)
  lval = 0
  count = 0
  if la!=lc:
    print 'BEEP! Two arrays must have the same length!'   
  if la == lc:
    l = la
    for i in range (0, l):
      lci = len(c[i])
      for j in range (0, lci):
        lval = lval + 1
    val = [None]*lval
    for i in range (0, l):
      lci = len(c[i])
      for j in range (0, lci):
        val[count] = a[i] * c[i][j]
        count = count + 1
    #print count
    #print lval
    #print val
    return val 
#End
###################################

###################################
############Evaluation#############
###################################
#Begin
def G_eval(g_arr_freq, g_arr_mnta, ext, g_freq_init,p,h,M,bet):
  val = 1
  delta = 0.1
  num = len(g_arr_freq[0])
  #print 'num = ',num
  for i in range (0, len(g_arr_freq)):
    e=Energy_eval(g_arr_mnta[i],g_freq_init,p,h,M)
    #print e
    #print abs(e + g_arr_freq[i][num-1]*1j*ext )

    if abs(e + g_arr_freq[i][num-1]*1j*ext) <=0.001:
      val = 0
      return val
    else:
      val = val/( e + g_arr_freq[i][num-1]*(ext*1j) )  
      #print '*****'
  return val
#End
###################################
def AMI_arrays_out(G_init,n):
  temp1_freq = [None]*20   # auxiliary variable
  temp1_mnta = [None]*20   # auxiliary variable
  temp2_freq = [None]*8   # auxiliary variable
  temp2_mnta = [None]*8   # auxiliary variable
  temp3 = [None]*8    # auxiliary variable
  S_list = [[None]*8 for i in range(0,n+1)]    # auxiliary variable
  P_list_freq = [[None]*8 for i in range(0,n+1)]    # auxiliary variable
  P_list_mnta = [[None]*8 for i in range(0,n+1)]    # auxiliary variable
  S_fin = [[None]*8 for i in range(0,n+1)]    # auxiliary variable
  P_fin_freq = [[None]*8 for i in range(0,n+1)]    # auxiliary variable
  P_fin_mnta = [[None]*8 for i in range(0,n+1)]    # auxiliary variable

  #bet = 5  # Temperature
  #mu=0    # Chemical potential
  #t=1       # Hopping amplitude
  
  
  #G_init = [[-1, 1, 0, 0, 1,],[-1, 1, 0, 0, 1],[1, 1, 0, 0, 0]]
  #G_init = data[:]
  G_sym = sym_rep_g(G_init)
  #print G_sym
  #print 'G_freq = ', G_sym[0]  # Find the frequency part of the multiplication of green's functions 
  #print 'G_mnta = ', G_sym[1]  # Find the momentum part of the multiplication of green's functions 
  #print '*********************'
  #one_sum_n_G(G_sym[0], G_sym[1], 1)

  #print Energy(G_sym[0],k_init,t,mu)

  #n=1

  if n==1:   # only one summation over the first frequency
    g=one_sum(G_sym[0], G_sym[1], 0)    # summation over the first frequency 
    ovl = g[5]
    #print 'overal sign = ', ovl
    temp2_freq = g[3]  # Store frequency part of the first summation
    temp2_mnta = g[4]  # Store momentum part of the first summation
    #print 'RESULT For n = ', n
    #print 'beta ^', n
    for i in range (0, len(g[0])):     # To ensure the correct sign 
      g[0][i] = g[0][i]*ovl
    #print g[0], '.', 'f', 'at frequency', g[1], 'and momenta', g[2]         # print signs and poles
    S_list[1][0] = g[0]
    P_list_freq[1][0] = g[1]
    P_list_mnta[1][0] = g[2]
    #print S_list
    #print P_list_freq
    #print P_list_mnta 
    #print '*'
    #print '['
    #for i in range (0, len(temp2_freq)): # print Green's function expression after the first summation
      #if temp2_freq[i] != None:
        #print 'R_{',i+1, '}' 
        #print temp2_freq[i]
        #print temp2_mnta[i]
        #if i<len(temp2_freq)-1:
          #print ','
    #print ']'     
    #print temp2              # print Green's function expression after the first summation
    #print len(temp2)
    #print '********'
  #print G[2][0]
  
        
  c_s_p = 0
  if n>=2:    # summation over n frequency with n>=2
    temp2 = one_sum(G_sym[0], G_sym[1], 0)
    temp2_freq = temp2[3]   # first summation: frequency part
    temp2_mnta = temp2[4]   # first summation: momentum part
    #print len(temp2)
    step = len(temp2_freq)   # number of terms after the first summation
    old_step = 0
    new_step = 0
    #print 'temp2_freq =', temp2[3], 'temp2_mnta =', temp2[4]
    count = 0
    #print '******'
    #print 'RESULT For n = ', n
    #print '******'
    #print 'beta ^', n
    g = one_sum(G_sym[0], G_sym[1], 0)         # do first summation and save the result in G
    ovl = g[5]
    #print 'overal sign = ', ovl
    for i in range (0, len(g[0])):     # To ensure the correct sign 
      g[0][i] = g[0][i]*ovl
    #print g[0], '.', 'frequency poles',g[1], 'mnta poles',g[2]              # signs and poles of the first summation
    S_list[1][c_s_p] = g[0]
    P_list_freq[1][c_s_p] = g[1]
    P_list_mnta[1][c_s_p] = g[2]
    #c_s_p = c_s_p + 1
    #print 'x'
    temp2_freq = g[3]                  # Green's fumctions expression after the first summation
    temp2_mnta = g[4]
    #print 'temp2_freq', temp2_freq
    #print 'temp2_mnta', temp2_mnta
    #print temp2_freq[0]
    #print temp2_mnta[0]
  
  
  
  
    for k in range (0,n-1):         # to do next summations
      c_s_p = 0  
    
      #print 'Second sum'
    
      for i in range (0, step):  # to do summations for each term separately
        g = one_sum(temp2_freq[i], temp2_mnta[i], k+1)   # summation over next frequency for i-th term
        ovl = g[5]
        #print 'overall_sign = ', ovl
        for i in range (0, len(g[0])):     # To ensure the correct sign 
          if g[0][i]!=None:
            g[0][i] = g[0][i]*ovl
        #print 'k =', k
        #print 'G[2] =', G[2]
        #print 
        #print 'len G[2] =', len(G[2])
        #print g[0]
        #print g[1]
        #print g[2]
        #print g[3]
        #print g[4]
        for j in range (0, len(g[3])):   # to save result of the k-th summation
          #print 'count =', count
          if k==n-3:
            old_step = old_step + 1
          if k==n-2:
            new_step = new_step + 1
          if j == 0:
          
            #print g[0], '.', 'frequency poles',g[1], 'mnta poles',g[2]              # signs and poles of the first summation
            S_list[k+2][c_s_p] = g[0]
  	    P_list_freq[k+2][c_s_p] = g[1]
            P_list_mnta[k+2][c_s_p] = g[2]
            c_s_p = c_s_p + 1
            #print
          
          temp1_freq[count] = (g[3][j])
          temp1_mnta[count] = (g[4][j])
        
          #print temp1[count]
          count = count + 1 
      #if k==n-2:    
        #print '*'
      #else:
        #print 'x'
      temp2_freq = deepcopy(temp1_freq)     # result of the k-th summation which will be used as input for (k+1)-th summation
      temp2_mnta = deepcopy(temp1_mnta)
      #print 'old_step = ', old_step
    
      #for i in range (0, new_step):  # to have a neat output
        #temp3[i] = temp2[i+old_step] 
      step = count
      count = 0
    #print new_step
    #print old_step
    parcham = new_step
    old_step = 0
    new_step = 0 
    length = 0
    #print temp3
    #print '['
    for i in range (0, len(temp2_freq)):
      if temp2_freq[i] != None:
        #print 'R_{',i+1, '}' 
        #print temp2_freq[i]
        #print temp2_mnta[i]
        #if i<parcham-1:
          #print ','
        length = length + 1  
      else:
        i = len(temp2_freq)+1  
    #print length
    #print ']'     
    #print temp3[0:length]  # print final Green's functions expression
    #print temp2 
    #print len(temp2)
    #print length
    #print("--- Construction_time ---", (time.time() - start_time))
    #print 'S_list = ', S_list
    #print 'P_list_freq = ', P_list_freq
    #print 'P_list_mnta = ', P_list_mnta
  

  #print (temp2[1])

  
  return S_list, P_list_freq, P_list_mnta, temp2_freq, temp2_mnta, G_sym

def AMI_out(S_list, P_list_freq, P_list_mnta, temp2_freq, temp2_mnta, G_sym, omg,k_init,n,bet,t,mu):
  
  #bet = 10  # Temperature
  #mu=0    # Chemical potential
  #t=1       # Hopping amplitude
  #e_ar = [0,1,0]
  #print Energy_eval(e_ar,G_sym[0],k,t,mu)
  poles = Poles(P_list_freq, P_list_mnta,G_sym[0],k_init,t,mu)  # construct joint poles
  #print 'joint poles = ', poles
  for i in range (0,len(poles)):
    for j in range (0,len(poles[i])):
      if poles[i][j]==None:
        j = len(poles[i])+1
      else:
        for k in range (0,len(poles[i][j])):
          #print poles[i][j][k]
          poles[i][j][k] = f(poles[i][j][k],bet)  # Find the numerical values of the poles 
        
  #print poles


  ########################################
  ################AMI Algebra#############
  ########################################

  res1 = dot_num(S_list[1][0], poles[1][0])
  temp=res1
  #print res1
  for m in range (2,n+1):
    o = 0
    #for k in range (2, n+1):
    for i in range(0, len(S_list[m])):
      if S_list[m][i] != None:
        o = o + 1
      else:
        i = len(poles[m]) + 1
    #print o
    S_new = [None]*(o)
    P_new = [None]*(o)

    for i in range(0, len(S_list[m])):
      if S_list[m][i] != None:
        S_new[i] = S_list[m][i]
        P_new[i] = poles[m][i]
  
      else:
        i = len(S_list[m]) + 1

    #print o
    #print S_new
    #print P_new  
    #print res1
    res2 = dot_arr(S_new, P_new)
    #print res2
    temp = cross(res1,res2)
    res1 = temp
    #print 'temp = ', res1
    #print len(res1)
    #print S_list
  #print("--- Construction_time ---", (time.time() - start_time))
  #start_time_new = time.time()

  #file = open('Out.dat', 'w')  # To store output

  #print 'length =', len(temp)

  #for i in range (-100, 100):
  nu = (2*omg+1)*pi/bet  # External frequency on matsubara axis
  #nu_10 = i*0.01
  #nu_5 = i*0.01
  G_val = 0.0  # To store numerical value at a given external frequency
  #print i
  for j in range (0, len(temp)):
      #G_val = G_eval(temp2[j], nu_5)
      #G_val = G_eval(temp2[j], nu_5) + G_val 
    G_val = temp[j]*G_eval(temp2_freq[j], temp2_mnta[j], nu, G_sym[0],k_init,t,mu,bet) + G_val 
    #txt = str(nu) + '\t' + str(G_val.real) + '\t' + str(G_val.imag) + '\n'
    #file.write(txt)
  #file.close()
  #print '**********************'
  #print("--- Evaluation_time ---", (time.time() - start_time_new))
  #print("--- Total_time ---", (time.time() - start_time))
  
  #print '**********************'
  #print("--- Total_time ---", (time.time() - start_time))
  '''
  print
  print 'Plotting..............'


  data1 = np.loadtxt('Out.dat') 
  plt.plot(data1[:,0], data1[:,1]) 
  #plt.title(r'$\vec k_1=(0,0), \vec k_2 =\vec k_3= \vec k_4 = \vec k_5 = \vec k_6 = \vec k_7 = \vec k_8 = \vec k_9 = \vec k_{10} = (\pi, \frac{\pi}{3}), \beta=10, \beta=10, \mu=1.1 $')
  #plt.legend(loc='best')
  plt.xlabel(r'$i\nu$')
  #plt.ylabel(r'Im{$I$}')
  plt.ylabel(r'$Re I^{(3)} $')
  plt.savefig('Re_3rd.eps', format='eps', dpi=1200, bbox_inches='tight')  # Save figure
  plt.show()
  plt.plot(data1[:,0], data1[:,2]) 
  #plt.title(r'$\vec k_1=(0,0), \vec k_2 =\vec k_3= \vec k_4 = \vec k_5 = \vec k_6 = \vec k_7 =  \vec k_8 = \vec k_9 = \vec k_{10} = (\pi, \frac{\pi}{3}), \beta=10, \beta=10, \mu=1.1 $')
  #plt.legend(loc='best')
  plt.xlabel(r'$i\nu$')
  #plt.ylabel(r'Im{$I$}')
  plt.ylabel(r'$Im I^{(3)} $')
  plt.savefig('Im_3rd.eps', format='eps', dpi=1200, bbox_inches='tight')  # Save figure
  plt.show() # Show figure
  ''' 
  return G_val.real,G_val.imag
'''
data = np.loadtxt('test.dat')  # reading data from data.txt
count = len(open('test.dat').readlines())
if count==1:
  G_init = [None]*count
  G_init[0]=data
if count>1:
  G_init = data[:]  # multiplication of n Green's function in the at the first stage of summation
  for i in range (0, count):
    G_init[i] = data[i,:]  # multiplication of n Green's function in the at the first stage of summation
print '*********************'
print '***Green_Functions***'
#print G      # print the multiplication of the given Green's functions
print '*********************' 
#G_init = [[-1, 1, 0, 0, 1,],[-1, 1, 0, 0, 1],[1, 1, 0, 0, 0]]
#G_init = data[:]
G_sym = sym_rep_g(G_init)
#print G_sym
print 'G_freq = ', G_sym[0]  # Find the frequency part of the multiplication of green's functions 
print 'G_mnta = ', G_sym[1]  # Find the momentum part of the multiplication of green's functions 
print '*********************'

#print der(G_sym[0][0],G_sym[1][0],0)
#print 'frequency part= ', der_mul(G_sym[0], G_sym[1], 0)[0]
#print 'mnta part= ', der_mul(G_sym[0], G_sym[1], 0)[1]
#dec = decision(G_sym[0],G_sym[1])
#print 'freq -->', dec[0]
#print 'mnta -->', dec[1]
#print 'multiplicity -->', dec[2]
#print 'overal sign -->', dec[3]

############################################
############One sum over frequency r########
############################################
#Begin
print 
print 'One sum (over nu1) result for one of the third order diagrams' 
print 'with an intrinsic double pole:'
print
r=0  # Choose frequency
output = one_sum(G_sym[0],G_sym[1],r)  # Construct S, P and R arrays
print 'sign --> ', output[0]
print 'roots frequency part --> ', output[1]
print 'roots momentum part --> ', output[2]
print 'R array frequency part --> ', output[3]
print 'R array momentum part --> ', output[4]
print 'Overal sign --> ', output[5]
print 'It is non-zero --> ', output[6]
print
############################################
'''
'''
ami_in=[[0, 1, 0, 0], [1, 0, 0, 0], [0, 1, -1, 1], [0, 1, 0, 0], [0, 0, 1, 0]]
symb=sym_rep_g(ami_in)
print symb[0]
print symb[1]
print symb[2]
b=AMI_arrays_out(ami_in,3)
s_list=deepcopy(b[0])
p_list_freq=deepcopy(b[1])
p_list_mnta=deepcopy(b[2])
r_freq=deepcopy(b[3])
r_mnta=deepcopy(b[4])
g_sym=deepcopy(b[5])
print 'S_list = ', s_list
print 'P_list_freq = ', p_list_freq
print 'P_list_mnta = ', p_list_mnta
print 'R_freq = ', r_freq
print 'R_mnta = ', r_mnta
print 'G_sym = ', g_sym
'''

# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 11:33:12 2016

@author: lxie2
"""

import numpy as np
from numpy import random
import os,pickle
from FunctionList import growth, random_mutate,comm_distribute
from scipy import integrate

# assign constants
u1max = 0.7
u2max = 0.3
beta1 = 0.1
Y1sN = 1e3
Y2sN = 1e3
Y1AN = 1e3
r2N = 3e-3
K2K1 = 1.

N0 = 150.
A0 = 0.
s0 = 10.
p0 = 0.
T0 = 102*0.15


N10 = 10.
N20 = N0-N10
comm_type_num = 10
comm_rep_num = 10
max_popul = 1e5
mut_rate = 1e-3
pcs = 1e-15
beta1max = 1
#SampleTP = 200;

t_bin = 0.15
t_binnum = 102

# initiates empty matrices. Each column will be manipulated seperately 
B1_L_all = np.zeros((max_popul,comm_type_num*comm_rep_num))
B2_L_all = np.zeros((max_popul,comm_type_num*comm_rep_num))
B1_beta1_all = np.zeros((max_popul,comm_type_num*comm_rep_num))
B1_t_all = np.zeros((t_binnum,comm_type_num*comm_rep_num))
B2_t_all = np.zeros((t_binnum,comm_type_num*comm_rep_num))
s_all = np.zeros((t_binnum,comm_type_num*comm_rep_num))
A_all = np.zeros((t_binnum,comm_type_num*comm_rep_num))
p_all = np.zeros((t_binnum,comm_type_num*comm_rep_num))
parentnum_all = np.zeros(comm_type_num*comm_rep_num)
rseed_all = np.uint32(np.zeros(comm_type_num*comm_rep_num))

B1_L_all[0:N10,:] = 1
B2_L_all[0:N20,:] = 1
B1_beta1_all[0:N10,:] = beta1
random.seed()
rseed_all = np.uint32(np.random.randint(2**31-1,size = comm_type_num*comm_rep_num))
    
for n in range(1):
    foldername = "C"+str(n+1)
    if not os.path.exists(foldername):
        os.makedirs(foldername)
    # manipulate each column 
    for rep in range(comm_type_num*comm_rep_num):
        nodeseed = rseed_all[rep]
        random.seed(nodeseed)
        B1_L = B1_L_all[:,rep]
        B2_L = B2_L_all[:,rep]
        B1_beta1 = B1_beta1_all[:,rep]
        B1_t = np.zeros(t_binnum)
        B2_t = np.zeros(t_binnum)
        s = np.zeros(t_binnum)
        A = np.zeros(t_binnum)
        p = np.zeros(t_binnum)
        
        growth_para = [u1max,u2max,Y1AN,Y1sN,Y2sN,r2N,K2K1,B1_L,B2_L,A0,s0,p0,t_bin,max_popul,pcs]
        [B1_L,B1_beta1,B2_L,A[0],s[0],p[0],B1_t[0],B2_t[0]] = growth(B1_L,B1_beta1,B2_L,growth_para)        
        mut_para = [mut_rate,max_popul,pcs]
        [B1_L,[B1_beta1]] = random_mutate(B1_L,[B1_beta1],mut_para)
        [B2_L,[]] = random_mutate(B2_L,[],mut_para)
            
        for dt in range(1,t_binnum):
            growth_para = [u1max,u2max,Y1AN,Y1sN,Y2sN,r2N,K2K1,B1_L,B2_L,A[dt-1],s[dt-1],p[dt-1],t_bin,max_popul,pcs]
            [B1_L,B1_beta1,B2_L,A[dt],s[dt],p[dt],B1_t[dt],B2_t[dt]] = growth(B1_L,B1_beta1,B2_L,growth_para)        
            mut_para = [mut_rate,max_popul,pcs]
            [B1_L,[B1_beta1]] = random_mutate(B1_L,[B1_beta1],mut_para)
            [B2_L,[]] = random_mutate(B2_L,[],mut_para)
            
        B1_L_all[:,rep] = B1_L
        B1_beta1_all[:,rep] = B1_beta1
        B2_L_all[:,rep] = B2_L
        B1_t_all[:,rep] = B1_t
        B2_t_all[:,rep] = B2_t
        s_all[:,rep] = s
        A_all[:,rep] = A
        p_all[:,rep] = p

    # save the data    
    filename = os.path.join( foldername, "comm_all" )
    np.savez(filename,B1_L_all = B1_L_all,B1_beta1_all = B1_beta1_all,B2_L_all = B2_L_all,\
    B1_t_all = B1_t_all,B2_t_all = B2_t_all,A_all = A_all,s_all = s_all,p_all = p_all,\
    parentnum_all = parentnum_all,rseed_all = rseed_all)

    rng = random.get_state()
    filename = os.path.join( foldername, "rng" )
    with open(filename, 'wb') as myfile:
        pickle.dump(rng, myfile, -1) #HIGHEST_PROTOCOL     
    # sort the data
    sort_idx = np.argsort(-p_all[-1,:])
    B1_L_sorted = B1_L_all[:,sort_idx]
    B1_beta1_sorted = B1_beta1_all[:,sort_idx]
    B2_L_sorted = B2_L_all[:,sort_idx]
    B1_t_sorted = B1_t_all[:,sort_idx]
    B2_t_sorted = B2_t_all[:,sort_idx]
    A_sorted = A_all[:,sort_idx]
    s_sorted = s_all[:,sort_idx]
    p_sorted = p_all[:,sort_idx]
    parentnum_sorted = parentnum_all[sort_idx]
    rseed_sorted = rseed_all[sort_idx]

    # generate matrices for next cycle of n
    rep_counter = 0
    gen_counter = 0
    for i in range(comm_type_num*comm_rep_num):
        if rep_counter >= comm_type_num*comm_rep_num:
            break
        dil_factor = np.int32(np.round((B1_t_sorted[-1,i]+B2_t_sorted[-1,i])/N0))
        if dil_factor==0:
            continue
        comm_all_idx = min(comm_type_num*comm_rep_num,rep_counter+min(dil_factor,comm_rep_num))
        dist_para = [comm_type_num,comm_rep_num,max_popul,dil_factor,rep_counter,pcs]
        [B1_L_all[:,rep_counter:comm_all_idx],
         B1_beta1_all[:,rep_counter:comm_all_idx],
         B2_L_all[:,rep_counter:comm_all_idx]]\
         = comm_distribute(B1_L_sorted[:,i],B1_beta1_sorted[:,i],B2_L_sorted[:,i],dist_para)
        parentnum_all[rep_counter:comm_all_idx] = i
        gen_counter = gen_counter+1
        rep_counter = rep_counter+min(dil_factor,comm_rep_num)
    rseed_all = np.uint32(random.randint(2**31-1,size = comm_type_num*comm_rep_num))    
    B1_L_gen = B1_L_sorted[:,:gen_counter]
    B1_beta1_gen = B1_beta1_sorted[:,:gen_counter]
    B2_L_gen = B2_L_sorted[:,:gen_counter]
    B1_t_gen = B1_t_sorted[:,gen_counter]
    B2_t_gen = B2_t_sorted[:,gen_counter]
    A_gen = A_sorted[:,gen_counter]
    s_gen = s_sorted[:,gen_counter]
    p_gen = p_sorted[:,gen_counter]
    parentnum_gen = parentnum_sorted[:gen_counter]
    rseed_gen = rseed_sorted[:gen_counter]
    
    filename = os.path.join( foldername, "comm_gen" )
    np.savez(filename,B1_L_gen = B1_L_gen,B1_beta1_gen = B1_beta1_gen,B2_L_gen = B2_L_gen,\
    B1_t_gen = B1_t_gen,B2_t_gen = B2_t_gen,A_gen = A_gen,s_gen = s_gen,p_gen = p_gen,\
    parentnum_gen = parentnum_gen,rseed_gen = rseed_gen) 

filename = os.path.join( foldername, "comm_all" )
np.savez(filename,B1_L_all = B1_L_all,B1_beta1_all = B1_beta1_all,B2_L_all = B2_L_all,\
B1_t_all = B1_t_all,B2_t_all = B2_t_all,A_all = A_all,s_all = s_all,p_all = p_all,\
parentnum_all = parentnum_all,rseed_all = rseed_all)

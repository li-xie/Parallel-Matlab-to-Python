# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 14:43:30 2016

@author: lxie2
"""
import numpy as np
from numpy import random
from scipy import integrate

def as_conc(t,y,para):

    AN = y[0];
    sN = y[1];

    u1max,u2max,Y1AN,Y1sN,Y2sN,r2N,K2K1,B1,B2 = para
    
    u1_coef = AN/(AN+sN)*sN/(sN+1)+sN/(AN+sN)*AN/(AN+1)
    u2_coef = sN/(sN+K2K1);
    dAN = -u1max/Y1AN*u1_coef*B1+u2max*r2N*u2_coef*B2
    dsN = -u1max/Y1sN*u1_coef*B1-u2max/Y2sN*u2_coef*B2

    dy = [dAN,dsN]
    return dy

def as_conc_jacobian(t,y,para):

    AN = y[0];
    sN = y[1];
    
    u1max,u2max,Y1AN,Y1sN,Y2sN,r2N,K2K1,B1,B2 = para
    
#    % u1_coef = AN/(AN+sN)*sN/(sN+1)+sN/(AN+sN)*AN/(AN+1);
#    % u2_coef = sN/(sN+K2K1);
#    % dAN = -u1max/Y1AN*u1_coef*B1+u2max*r2N*u2_coef*B2;
#    % dsN = -u1max/Y1sN*u1_coef*B1-u2max/Y2sN*u2_coef*B2;
    
    
    du1_coefdAN = sN**2/(AN+sN)**2/(sN+1)+(sN**2-sN*AN**2)/(AN+sN)**2/(AN+1)**2;
    du1_coefdsN = AN**2/(AN+sN)**2/(AN+1)+(AN**2-AN*sN**2)/(AN+sN)**2/(sN+1)**2;
    du2_coefdsN = K2K1/(sN+K2K1)**2;
    
    dANdAN = -u1max/Y1AN*du1_coefdAN*B1;
    dANdsN = -u1max/Y1AN*du1_coefdsN*B1+u2max*r2N*du2_coefdsN*B2;
    dsNdAN = -u1max/Y1sN*du1_coefdAN*B1;
    dsNdsN = -u1max/Y1sN*du1_coefdsN*B1-u2max/Y2sN*du2_coefdsN*B2;
    return [[dANdAN, dANdsN],[dsNdAN, dsNdsN]]
    
def growth(B1_L,B1_beta1,B2_L,para):
    u1max,u2max,Y1AN,Y1sN,Y2sN,r2N,K2K1,B1_L,B2_L,A0,s0,p0,t_bin,max_popul,pcs = para
    B1_L_out = np.zeros(max_popul)
    B1_beta1_out = np.zeros(max_popul)
    B2_L_out = np.zeros(max_popul)
    B1_counter = np.count_nonzero(B1_L>pcs)
    B2_counter = np.count_nonzero(B2_L>pcs)
    B1_L_temp = B1_L[:B1_counter]    
    B2_L_temp = B2_L[:B2_counter] 
    B1_beta1_temp = B1_beta1[:B1_counter]
    ode_para = [u1max,u2max,Y1AN,Y1sN,Y2sN,r2N,K2K1,np.sum(B1_L_temp),np.sum(B2_L_temp)]
    r = integrate.ode(as_conc,as_conc_jacobian).set_integrator('vode',rtol = 1e-6, method = 'bdf')
    r.set_initial_value([0.,s0],0.).set_f_params(ode_para).set_jac_params(ode_para)
    tx = [0.]
    y = [np.array([0.,10.])]
    while r.successful() and r.t<t_bin:
        r.integrate(t_bin,step = True)
        tx.append(r.t)
        y.append(r.y)
    y = np.array(y)
    tx = np.array(tx)
    y = y[np.nonzero(tx<=t_bin)]
    tx = tx[np.nonzero(tx<=t_bin)]
    tx = np.append(tx,np.array([t_bin]),axis = 0)
    y = np.append(y,[r.integrate(t_bin)],axis = 0)
    AN = y[:,0]
    sN = y[:,1]
    u1_coef = AN/(AN+sN)*sN/(sN+1)+sN/(AN+sN)*AN/(AN+1);
    u2_coef = sN/(sN+K2K1);
    u1 = np.trapz(u1_coef,tx)*u1max;
    u2 = np.trapz(u2_coef,tx)*u2max;
    B1_L_out[0:B1_counter] = np.exp((1-B1_beta1_temp)*u1)*B1_L_temp;
    B2_L_out[0:B2_counter] = np.exp(u2)*B2_L_temp;
    p_out = sum(B1_beta1_temp*(B1_L_out[0:B1_counter]-B1_L_temp)/(1-B1_beta1_temp));
    A_out = AN[len(tx)-1];
    s_out = sN[len(tx)-1];
    B1_beta1_out[0:B1_counter] = B1_beta1_temp
    B1_t = sum(B1_L);
    B2_t = sum(B2_L);
    return [B1_L_out,B1_beta1_out,B2_L_out,A_out,s_out,p_out,B1_t,B2_t]
    
def mu_factor(n):
    u = random.rand(n)
    y = np.zeros(n)
    idx = np.nonzero(u>1-1e-11)
    y[idx] = 2e9*(u[idx]-1)+0.14
    idx = np.nonzero(np.logical_and(u <= 1-1e-11,u>1-1e-11-3e-6))
    y[idx] = 2e4*(u[idx]-1+1e-11)+0.12
    idx = np.nonzero(u <= 1-1e-11-3e-6)
    y[idx] = -np.log(1-u[idx])/260
    return y
    
def mu_spontaneous(n):
    u = random.rand(n)
    y = np.zeros(n)
    idx = np.nonzero(u<0.001)
    if np.array(idx).size!=0:
        y[idx] = -1
    idx = np.nonzero(np.logical_and(u<=0.95,u>0.9))
    if np.array(idx).size!=0:
        y[idx] = mu_factor(np.array(idx).size)
    idx = np.nonzero(u>0.95)
    if np.array(idx).size!=0:
        y[idx] = -1.*mu_factor(np.array(idx).size)
    return y    
    
def random_mutate(B_L,B_pheno,mut_para):
    mut_rate,max_popul,pcs = mut_para
    B_counter = np.count_nonzero(B_L>pcs)
    pheno_count = len(B_pheno)
    
    if pheno_count>0:
        B_pheno_out = [np.zeros(max_popul)]*pheno_count    
        
        div_idx = np.nonzero(B_L >= 2)[0]
        div_len = len(div_idx)
        if div_len  >0:
            for i in range(pheno_count):
                pheno_temp = B_pheno[i]
                mut_marker = np.array(random.rand(div_len)<=mut_rate).astype(float)
                mut_multiplier = mu_spontaneous(div_len)*mut_marker+1
                pheno_temp[div_idx] =  pheno_temp[div_idx]*mut_multiplier
                mut_marker = np.array(random.rand(div_len)<=mut_rate).astype(float)
                mut_multiplier = mu_spontaneous(div_len)*mut_marker+1
                pheno_temp[B_counter:B_counter+div_len] = pheno_temp[div_idx]*mut_multiplier
                B_pheno_out[i] = pheno_temp
            B_L[B_counter:B_counter+div_len] =  B_L[div_idx]/2
            B_L[div_idx] =  B_L[div_idx]/2
    else:
        div_idx = np.nonzero(B_L>=2)[0]
        div_len = len(div_idx)
        if div_len>0:
            B_L[B_counter:B_counter+div_len] = B_L[div_idx]/2
            B_L[div_idx] = B_L[div_idx]/2
            B_counter = B_counter+div_len
        B_pheno_out = []
    return [B_L,B_pheno_out]
    

def comm_distribute(B1_L,B1_beta1,B2_L,para):
    comm_type_num = para[0]
    comm_rep_num = para[1]
    max_popul = para[2]
    dil_factor = para[3]
    rep_counter = para[4]
    pcs = para[5]    
    
    B1_L_live = B1_L[B1_L>pcs]
    B1_beta1_live = B1_beta1[B1_L>pcs]
    B2_L_live = B1_L[B2_L>pcs]    
    
    B1_counter = len(B1_L_live)
    B2_counter = len(B2_L_live)
    rand_temp = np.floor(random.rand(B1_counter+B2_counter)*dil_factor)
    B1_rand = rand_temp[:B1_counter]
    B2_rand = rand_temp[B1_counter:]
    if dil_factor >= comm_rep_num:
        B1_L_new = np.zeros((max_popul,comm_rep_num))
        B1_beta1_new = np.zeros((max_popul,comm_rep_num))
        B2_L_new = np.zeros((max_popul,comm_rep_num))
        for i in range(comm_rep_num):
            B1_num = np.count_nonzero(B1_rand==i)
            if B1_num >= 1:
                B1_L_new[:B1_num,i] = B1_L_live[B1_rand == i]
                B1_beta1_new[:B1_num,i] = B1_beta1_live[B1_rand == i]
            B2_num = np.count_nonzero(B2_rand == i)
            if B2_num >= 1:
                B2_L_new[:B2_num,i] = B2_L_live[B2_rand == i]
        if rep_counter+comm_rep_num>comm_type_num*comm_rep_num:
            B1_L_new = B1_L_new[:,:comm_type_num*comm_rep_num-rep_counter]
            B1_beta1_new = B1_beta1_new[:,:comm_type_num*comm_rep_num-rep_counter]
            B2_L_new = B2_L_new[:,:comm_type_num*comm_rep_num-rep_counter]
    else:
        B1_L_new = np.zeros((max_popul,dil_factor))
        B1_beta1_new = np.zeros((max_popul,dil_factor))
        B2_L_new = np.zeros((max_popul,dil_factor))
        for i in range(dil_factor):
            B1_num = np.count_nonzero(B1_rand==i)
            if B1_num >= 1:
                B1_L_new[:B1_num,i] = B1_L_live[B1_rand == i]
                B1_beta1_new[:B1_num,i] = B1_beta1_live[B1_rand == i]
            B2_num = np.count_nonzero(B2_rand == i)
            if B2_num >= 1:
                B2_L_new[:B2_num,i] = B2_L_live[B2_rand == i]
        if rep_counter+dil_factor > comm_type_num*comm_rep_num:
            B1_L_new = B1_L_new[:,:comm_type_num*comm_rep_num-rep_counter]
            B1_beta1_new = B1_beta1_new[:,:comm_type_num*comm_rep_num-rep_counter]
            B2_L_new = B2_L_new[:,:comm_type_num*comm_rep_num-rep_counter]
    
    return [B1_L_new,B1_beta1_new,B2_L_new]
        
    
    
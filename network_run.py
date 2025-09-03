##################################################################################
# network_run.py -- Uses the simulator from network_simulator.py
# and runs a simulation of a plastic recurrent network as in:
#
# Ref: Sadeh, Clopath and Rotter (PLOS Computational Biology, 2015).
# Emergence of Functional Specificity in Balanced Networks with Synaptic Plasticity.
#
# Author: Sadra Sadeh <s.sadeh@ucl.ac.uk> // Created: 2014-2015
##################################################################################

import numpy as np
import pylab as pl
import pickle as cPickle
from imp import reload
import params; reload(params); from params import *
import network_simulator as NS

#################################################################################
# -- generating the weight matrix
w0_exc = np.concatenate(( J*np.random.binomial(1, eps_ee, (ne,ne)), \
                          J*np.random.binomial(1, eps_ei, (ne,ni)) ), 1)
w0_inh = np.concatenate(( -g*J*np.random.binomial(1, eps_ie, (ni,ne)), \
                          -g*J*np.random.binomial(1, eps_ii, (ni,ni)) ), 1)
W0 = np.concatenate((w0_exc , w0_inh))

#################################################################################
# -- before learning
print('### before plasticity')

x_bp = []
sim_time_test = sim_time
x_len_test = sim_time_test/dt
for nn in range(n):
    if nn < ne:
        p_rate = b_rate*(1+m_exc*np.cos(2*(th - po_init[nn])))
    else:
        p_rate = b_rate*(1+m_inh*np.cos(2*(th - po_init[nn])))
    try:
        rate_ev = np.random.poisson(lam=int(p_rate * dt / 1000), size=int(x_len_test)).tolist()
    except TypeError as e:
        print(f"TypeError: {e}")
        print(f"p_rate * dt / 1000 evaluated to {p_rate * dt / 1000} which is cast to {int(p_rate * dt / 1000)}")
        print(f"Size parameter x_len_test is {x_len_test}, should be integer but was cast to {int(x_len_test)}")
    x_bp.append(rate_ev)
x_bp = np.array(x_bp)

x_ap = np.copy(x_bp)

y0 = np.zeros((1,n)) 

y_bp, s_bp, ym_plst_bp, yp_plst_bp, y_avg_bp, Wf_bp = \
        NS._net_sim_(A = A, y0 = y0, x = x_bp, vth = vth, W0 = W0, synapse='static')

spk_bp = np.where(s_bp[0:n,:] != 0)


#################################################################################
# -- during learning
print('### within plasticity')

W_blk = W0
spk_wp_tot = []
W_blk_tot = []
stim_rng_tot = []
for blk in range(block_no):
    print(blk)
    stim_rng = np.random.uniform(0, np.pi, int(stim_no))
    stim_rng_tot.append(stim_rng)
    t_stim = sim_time / len(stim_rng)
    x_wp = []
    for nn in range(n):
        rates = []
        for st in stim_rng:
            if nn < ne: p_rate = b_rate*(1+m_exc*np.cos(2*(st - po_init[nn])))
            else: p_rate = b_rate*(1+m_inh*np.cos(2*(st - po_init[nn])))
            try:
                rate_ev = np.random.poisson(lam=int(p_rate * dt / 1000), size=int(x_len / len(stim_rng))).tolist()
            except TypeError as e:
                print(f"TypeError: {e}")
        rates = np.array(rates)
        x_wp.append(rates)
    x_wp = np.array(x_wp)

    y0 = np.zeros((1,n)) 

    y, s_wp, ym_plst, yp_plst, y_avg, W_blk = \
       NS._net_sim_(A = A, y0 = y0, x = x_wp, vth = vth, W0 = W_blk, synapse='plastic')
    spk_wp = np.where(s_wp[0:n,:] != 0)
    spk_wp_tot.append(spk_wp)
    W_blk_tot.append(W_blk)

stim_rng_tot = np.array(stim_rng_tot)
W_blk_tot = np.array(W_blk_tot)

Wf = W_blk

#################################################################################
# -- after learning
print('### after plasticity')

y0 = np.zeros((1,n)) 

y_ap, s_ap, ym_plst_ap, yp_plst_ap, y_avg_ap, Wf_ap = \
        NS._net_sim_(A = A, y0 = y0, x = x_ap, vth = vth, W0 = Wf, synapse='static')

spk_ap = np.where(s_ap[0:n,:] != 0)


#################################################################################
# -- spontaneous activity
spont_act = 0 # set it to 1 to simulate spontaneous activity
if spont_act:
    
    print('### spontaneous activity')

    block_no_sp = 10 

    W_sp = Wf

    W_sp_tot = []
    spk_sp_tot = []
    ## plastic spontaneous
    for blk in range(block_no_sp):
        print(blk)
        stim_rng = np.arange(1, stim_no+1)
        #stim_rng_tot.append(stim_rng)
        t_stim = sim_time / len(stim_rng)
        x_wp = []
        #np.random.seed(1234)
        for nn in xrange(n):
            rates = []
            for st in stim_rng:
                p_rate = b_rate/2#*st/stim_no
                rates = rates + np.random.poisson(p_rate*dt/1000., x_len/len(stim_rng)).tolist()
            rates = np.array(rates)
            x_wp.append(rates)
        x_wp = np.array(x_wp)

        y0 = np.zeros((1,n)) 

        y, s_sp, ym_plst, yp_plst, y_avg, W_sp = \
           NS._net_sim_(A = A, y0 = y0, x = x_wp, vth = vth, W0 = W_sp, synapse='plastic', inh_ltd=1)
        spk_sp = np.where(s_sp[0:n,:] != 0)
        spk_sp_tot.append(spk_sp)
        W_sp_tot.append(W_sp)

#################################################################################
# -- over-representing cardinal orientations
card_act = 0 # 1 simulates an over-representation of cardinal stimulus orientations
if card_act:
    
    print('### cardinal orientations')

    block_no_cd = 20

    W_cd = W0
    W_cd_tot = []
    spk_cd_tot = []
    stim_rng_tot_cd = []
    ## plastic spontaneous
    for blk in range(block_no_cd):
        print(blk)
        stim_rng = np.concatenate( (np.random.uniform(0, np.pi, stim_no/2), np.ones(stim_no/4)*0., np.ones(stim_no/4)*np.pi/2) )
        np.random.shuffle(stim_rng)
        stim_rng_tot_cd.append(stim_rng)
        t_stim = sim_time / len(stim_rng)
        x_cd = []
        for nn in xrange(n):
            rates = []
            for st in stim_rng:
                if nn < ne: p_rate = b_rate*(1+m_exc*np.cos(2*(st - po_init[nn])))
                else: p_rate = b_rate*(1+m_inh*np.cos(2*(st - po_init[nn])))
                rates = rates + np.random.poisson(p_rate*dt/1000., x_len/len(stim_rng)).tolist()
            rates = np.array(rates)
            x_cd.append(rates)
        x_cd = np.array(x_cd)

        y0 = np.zeros((1,n)) 

        y, s_cd, ym_plst, yp_plst, y_avg, W_cd = \
           NS._net_sim_(A = A, y0 = y0, x = x_cd, vth = vth, W0 = W_cd, synapse='plastic')
        spk_cd = np.where(s_cd[0:n,:] != 0)
        spk_cd_tot.append(spk_cd)
        W_cd_tot.append(W_cd)

#################################################################################
# -- saving the results

res_save = 1
if res_save:
    results = {}
    results['W0'] = W0

    results['spk_bp'] = spk_bp
    results['spk_wp_tot'] = spk_wp_tot
    results['spk_ap'] = spk_ap
    results['W_blk_tot'] = W_blk_tot
    results['stim_rng_tot'] = stim_rng_tot
    
    if spont_act:
        results['spk_sp_tot'] = spk_sp_tot
        results['W_sp_tot'] = W_sp_tot
    else:
        results['spk_sp_tot'] = []
        results['W_sp_tot'] = []
        
    if card_act:
        results['spk_cd_tot'] = spk_cd_tot
        results['W_cd_tot'] = W_cd_tot
        results['stim_rng_tot_cd'] = stim_rng_tot_cd
    else:
        results['spk_cd_tot'] = []
        results['W_cd_tot'] = []
        results['stim_rng_tot_cd'] = []

    fl = open('results', 'wb')
    cPickle.dump(results, fl)
    fl.close()

# ------

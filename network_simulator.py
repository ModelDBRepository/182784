##################################################################################
# network_simulator.py -- Network simulator to simulate
# plastic recurrent networks studied in:
#
# Ref: Sadeh, Clopath and Rotter (PLOS Computational Biology, 2015).
# Emergence of Functional Specificity in Balanced Networks with Synaptic Plasticity.
#
# Author: Sadra Sadeh <s.sadeh@ucl.ac.uk> // Created: 2014-2015
#################################################################################

import numpy as np
import pylab as pl

from imp import reload
import params; reload(params); from params import *

def _rect_(xx): return xx*(xx>0)

# -----------------------
# --- network simulator function: _net_sim_()
# simulates a network of N integrate-and-fire neurons with plastic synapses using Exact Integration

# -----------------------
# - takes as the input:
# A: coefficient matrix of the subthreshold dynamics (needed for Exact Integration) []
# y0: vector of initial membrane potential of neurons
# x: sequence of input spike trains for all neurons 
# W0: initial weight matrix
# synapse: the type of recurrent synapses (static / plastic)
# -----------------------
# - returns as the output:
# y: vector of membrane potentials for all neurons across time
# s: vector of spikes (0: no spike, 1: spike) for all neurons across time
# ym_plst, yp_plst: low-pass filtered versions of membrane potential needed for the (voltage-dependent) plasticity rule 
# y_avg: mean depolarization of the post-synaptic neuron
# W: final weight matrix
# -----------------------

def _net_sim_(A, y0, x, vth, W0, synapse='static'):
    B = np.exp(A*dt)
    W = np.copy(W0)    
    y = np.zeros(np.shape(x)) 
    ym_plst, yp_plst = np.zeros(np.shape(x)), np.zeros(np.shape(x))
    y_avg = np.zeros(np.shape(x))
    dw_m, dw_p = np.zeros(np.shape(x)), np.zeros(np.shape(x))
    dw_m_tot, dw_p_tot = [], []
    s, x_plst = np.zeros(np.shape(x)), np.zeros(np.shape(x))      

    for i in range(np.size(x,1)):
        if i == 0:
           y[:,0] = B*y0 + x[:,0]
           ym_plst[:,0] = Bm_plst*y0 + y0/tm_plst
           yp_plst[:,0] = Bp_plst*y0 + y0/tp_plst
        else:
           y[:,i] = B*y[:,i-1] + x[:,i] + s[:,i-1]*np.matrix(W)
           ym_plst[:,i] = Bm_plst*ym_plst[:, i-1] + y[:,i-1]/tm_plst
           yp_plst[:,i] = Bp_plst*yp_plst[:, i-1] + y[:,i-1]/tp_plst
           x_plst[:,i] = Bx_plst*x_plst[:,i-1] + s[:,i-1]/tx_plst
        
        s[:,i] = (y[:,i] > vth)
        y[:,i] *= (y[:,i] <= vth)
        
        if synapse == 'plastic':
            # running avg of mem pot
            t_avg = .1 # (s)
            ind_avg = t_avg*1000/dt
            i_avg = i-ind_avg
            if i_avg < 0: i_avg = 0
            if i != 0: y_avg[:,i] = np.mean(y[:,i_avg:i], 1)

            #### depression term
            dw_m = -dt * np.matrix(s[:,i]).T * np.matrix((A_ltd * y_avg[:,i]**2/70.) * _rect_(ym_plst[:,i] - vth_m))
            #### potentiation term
            dw_p = dt * A_ltp * np.matrix(x_plst[:,i]).T * np.matrix(_rect_(y[:,i] - vth_p) * _rect_(yp_plst[:,i] - vth_m) )
            
            ## exc to all
            W[0:ne,0:n] = W[0:ne,0:n] + np.array(dw_m[0:ne,0:n] + dw_p[0:ne,0:n])* (W[0:ne,0:n] != 0)
                ## bounds						
            W[0:ne] = (W[0:ne] >= w_max) * w_max + (W[0:ne] < w_max)* W[0:ne]
            W[0:ne] *= (W[0:ne] > 0)
            
            ## inh to exc
            W[ne:, 0:ne] = W[ne:, 0:ne] - dw_p[ne:, 0:ne] - dw_m[ne:, 0:ne]
                ## bounds
            W[ne:] = (W[ne:] <= w_max_inh) * w_max_inh + (W[ne:] > w_max_inh)* W[ne:]
            W[ne:] *= (W[ne:] <= 0)

    return y, s , ym_plst, yp_plst, y_avg, W

# ----

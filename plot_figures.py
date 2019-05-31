##################################################################################
# plot_figures.py -- Reads and analyzes the results generated from network_run.py
# and plots Figures 1 and 3 in:
#
# Ref: Sadeh, Clopath and Rotter (PLOS Computational Biology, 2015).
# Emergence of Functional Specificity in Balanced Networks with Synaptic Plasticity.
#
# Author: Sadra Sadeh <s.sadeh@ucl.ac.uk> // Created: 2014-2015
##################################################################################

import numpy as np
import pylab as pl
from imp import reload
import params; reload(params); from params import *
import pickle as cPickle

from mpl_toolkits.axes_grid.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable

################################################################################
# --- read the results

fl = open('results', 'rb')
results = cPickle.load(fl)
fl.close()

W0 = results['W0'] 

spk_bp = results['spk_bp']
spk_wp_tot = results['spk_wp_tot'] 
spk_ap = results['spk_ap'] 
spk_sp_tot = results['spk_sp_tot']
spk_cd_tot = results['spk_cd_tot']

W_blk_tot = results['W_blk_tot']
stim_rng_tot = results['stim_rng_tot']
W_sp_tot = results['W_sp_tot']
W_cd_tot = results['W_cd_tot']
stim_rng_tot_cd = results['stim_rng_tot_cd']

Wf = W_blk_tot[-1]


################################################################################
### Figure 1
#########

mksz = 2.
Ts = sim_time/1000.

def _temp_plot_(spk, ax, stim=0, yy=0):

	exid = np.where(spk[0] < ne)[0]
	inid = np.where(spk[0] >= ne)[0]

	htex = pl.histogram(spk[1][exid], bins=sim_time/10, range=(0, sim_time))
	htin = pl.histogram(spk[1][inid], bins=sim_time/10, range=(0, sim_time))

	hr = pl.histogram(spk[0], bins=n, range=(0, n)) 

	ax.plot(spk[1][exid]*dt, spk[0][exid], 'r.', markersize=mksz, label='Exc: '+ str(np.round(len(exid)/ne)) )
	ax.plot(spk[1][inid]*dt, spk[0][inid], 'b.', markersize=mksz, label='Inh: '+str(np.round(len(inid)/ne)) )

	ax.set_yticks([0, 99, 199, 299, ne-1, n-1])
	ax.set_yticklabels([])
	ax.set_ylim(0-10, n+10)
	ax.set_xlim([0-10, sim_time+10])

	ax.set_xticklabels([])	

	divider = make_axes_locatable(ax)	

	axHisty = divider.append_axes("right", size=.5, pad=0.1)
	#adjust_spines(axHisty,['left', 'bottom'], outward=0)

	axHisty.plot(hr[0]/(Ts), hr[1][0:-1], color='k', lw=2)
	
	if stim == 0: 
	   pl.text(.85, .5, str(np.round(len(exid)/ne / Ts,1)) +' Hz', transform = axHisty.transAxes, color='r')
	   pl.text(.85, .85, str(np.round(len(inid)/ni /Ts,1))+' Hz', transform = axHisty.transAxes, color='b')
	else:
	    pl.text(.9, .5, str(np.round(len(exid)/ne /Ts,1)) +' Hz', transform = axHisty.transAxes, color='r')
	    pl.text(.9, .85, str(np.round(len(inid)/ni /Ts,1))+' Hz', transform = axHisty.transAxes, color='b')

	axHisty.set_yticks([0, 99, 199, 299, ne-1, n-1])
	axHisty.set_yticklabels([])
	axHisty.set_xticks([0, 10])

	axHisty.set_ylim(0-10, n+10)

	axHistx = divider.append_axes("bottom", 1.2, pad=0.3)
	#adjust_spines(axHistx,['left', 'bottom'], outward=0)

	axHistx.plot(htex[1][0:-1], htex[0], color='r', lw=2, label='Exc')
	axHistx.plot(htin[1][0:-1], htin[0], color='b', lw=2, label='Inh')	
	
	axHistx.set_yticks([0, 50, 100, 150])
	axHistx.set_yticklabels([])
	
	if yy == 1:	
	   axHistx.set_xlabel('Time (ms)')
	   axHistx.set_ylabel('Population spike count')
	   axHistx.set_yticklabels([0, 50, 100, 150])
	   pl.legend(loc=1, frameon=False, prop={'size':12.5})
	   
	   axHisty.set_xlabel('Firing rate \n (spikes/s)', size=10)


fig = pl.figure(figsize=(16,8))

mycl = pl.imshow(np.random.uniform(0, 1, (100, 100)), cmap='hsv', vmin=0, vmax=1)
pl.clf()

##
ax1 = pl.subplot(141)
pl.title('Before Plasticity')
_temp_plot_(spk=spk_bp, ax=ax1, yy=1)
ax1.set_ylabel('Neuron #')
#ax1.set_yticks([0, 99, 199, 299, ne, n])
ax1.set_yticklabels([1, 100, 200, 300, ne, n])

for i in range(int(stim_no)): 
    ax1.plot([i*trial_time, (i+1)*trial_time], [-5, -5], '-', color=pl.cm.hsv(th/np.pi), lw=10)

##
ax2 = pl.subplot(142)
pl.title('Beginning of Plasticity')
_temp_plot_(spk=spk_wp_tot[0], ax=ax2, stim=1)

for i in range(int(stim_no)):
    clbr = ax2.plot([i*trial_time, (i+1)*trial_time], [-5, -5], '-', color=pl.cm.hsv(stim_rng_tot[0][i]/np.pi), lw=10)

cax = fig.add_axes([.485, .25, .01, .1])
clbr = pl.colorbar(mycl, cax=cax, orientation='vertical')
clbr.set_ticks([0, .25, .5, .75, 1])
clbr.set_ticklabels([0, 45, 90, 135, 180])

##
ax3 = pl.subplot(143)
pl.title('End of Plasticity')

_temp_plot_(spk=spk_wp_tot[block_no-1], ax=ax3, stim=2)

for i in range(int(stim_no)):
    ax3.plot([i*trial_time, (i+1)*trial_time], [-5, -5], '-', color=pl.cm.hsv(stim_rng_tot[-1][i]/np.pi), lw=10)

##
ax4 = pl.subplot(144)
pl.title('After Plasticity')
_temp_plot_(spk=spk_ap, ax=ax4)

for i in range(int(stim_no)):
    ax4.plot([i*trial_time, (i+1)*trial_time], [-5, -5], '-', color=pl.cm.hsv(th/np.pi), lw=10)

ax4.text(.25, .1, 'Sparser activity', size=15, transform = ax4.transAxes)

pl.subplots_adjust(left=.05, right=.95, bottom=.075, top=.95, wspace=.25)

pl.savefig('Fig1')

################################################################################
### Figure 3 (A-C)
###########

pl.figure(figsize=(14,5))

pl.subplot(131)
pl.title('Initial Weights (W0)')
pl.imshow(W0)
clb = pl.colorbar(shrink=.75)
clb.set_ticks([-4, 0, .5])
clb.set_ticklabels([-4, 0, .5])

pl.xlabel('Post-synaptic #')
pl.ylabel('Pre-synaptic #')

pl.subplot(132)
pl.title('Final Weights (Wf)')
pl.imshow(Wf, vmin=-5, vmax=2)
clb=pl.colorbar(shrink=.75)
clb.set_ticks([-5, -4, -3, -2, -1, 0, 1, 2])
clb.set_ticklabels([-5, -4, -3, -2, -1, 0, 1, 2])

pl.subplot(133)
pl.title('Weight Changes (Wf - W0)')
pl.imshow(Wf - W0, vmin=-1.5, vmax=1.5)
clb = pl.colorbar(shrink=.75)
clb.set_ticks([-1.5, -1, -.5, 0, .5, 1, 1.5])
clb.set_ticklabels([-1.5, -1, -.5, 0, .5, 1, 1.5])

pl.savefig('Fig3A')

###

## bid and fs
dW = Wf - W0
dpo, dw_po, Wf_po, W0_po = [], [], [], []
for ii in range(ne):
    for jj in range(ne):
        Wf_po.append(Wf[ii, jj])
        W0_po.append(W0[ii, jj])   
        dpo.append(po_init[ii] - po_init[jj])

W0_po, Wf_po = np.array(W0_po), np.array(Wf_po)
dw_po = np.array(dw_po)
dpo = np.array(dpo)

dpo_id1 = np.where( (abs(dpo) < np.pi/6) + (abs(dpo) > np.pi-np.pi/6) == True)
dpo_id2 = np.where( ((abs(dpo) < 2*np.pi/6) * (abs(dpo) > np.pi/6)) +  ((abs(dpo) > np.pi-2*np.pi/6) * (abs(dpo) < np.pi-np.pi/6)) )
dpo_id3 = np.where( ((abs(dpo) < 3*np.pi/6) * (abs(dpo) > 2*np.pi/6)) +  ((abs(dpo) > np.pi-3*np.pi/6) * (abs(dpo) < np.pi-2
*np.pi/6)) )

################################################################################
### Figure 3 (D-G)
###########

pl.figure(figsize=(14,5))

#pl.title('aligend weights')
def _plot_alignw_(rng1=range(0,ne), rng2=range(ne,n)):	
    for i in rng1:
        ax.plot(-po_init[rng2] + po_init[i], Wf[i][rng2], 'k.', ms=mksz, alpha=.5)
    for i in rng1:
        if i == rng1[0]: ax.plot(-po_init[rng2] + po_init[i], W0[i][rng2], 'r.', ms=1, alpha=1.)
        else: ax.plot(-po_init[rng2] + po_init[i], W0[i][rng2], 'r.', ms=mksz, alpha=.5)


ax = pl.subplot(141)
pl.title('Exc to Exc')
_plot_alignw_(rng1=range(0,ne), rng2=range(0,ne))
ax.set_xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
ax.set_xticklabels(['', -90, 0, 90, ''])
ax.set_yticks([0, .5, 1, 1.5, 2])
ax.set_yticklabels([0, .5, 1, 1.5, 2])
ax.set_ylim([0-.1, 2+.1])
ax.set_ylabel('Final Weights (mV)')
ax.set_xlabel('Pre-syn. PO - Post-syn. PO (deg)')

ax = pl.subplot(142)
pl.title('Exc to Inh')
_plot_alignw_(rng1=range(0,ne), rng2=range(ne,n))
ax.set_xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
ax.set_xticklabels(['', -90, 0, 90, ''])
ax.set_yticks([0, .5, 1, 1.5, 2])
ax.set_yticklabels([0, .5, 1, 1.5, 2])
ax.set_ylim([0-.1, 2+.1])
ax.set_ylabel('Final Weights (mV)')

ax = pl.subplot(143)
pl.title('Inh to Exc')
_plot_alignw_(rng1=range(ne,n), rng2=range(0,ne))
ax.set_xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
ax.set_xticklabels(['', -90, 0, 90, ''])
ax.set_yticks([-6, -5.5, -5, -4.5, -4, -3.5, -3])
ax.set_yticklabels([-6, -5.5, -5, -4.5, -4, -3.5, -3])
ax.set_ylim([-5.5-.1, -3+.1])
ax.set_ylabel('Final Weights (mV)')

ax=pl.subplot(144)
#adjust_spines(ax,['left', 'bottom'], outward=0, s=0)
pl.title('Connection Specificity')

xrng = [1, 2, 3]

def _part_dpo_(www, cl='k', lbl=[]):
	ym1, ys1 = np.mean(www[dpo_id1]), np.std(www[dpo_id1])
	ym2, ys2 = np.mean(www[dpo_id2]), np.std(www[dpo_id2])
	ym3, ys3 = np.mean(www[dpo_id3]), np.std(www[dpo_id3])

	ymrng = np.array([ym1, ym2, ym3])
	ysrng = np.array([ys1, ys2, ys3])

	ax.plot(xrng, ymrng, '-o', lw=2, color=cl, label=lbl)

_part_dpo_(W0_po, cl='r', lbl='Initial')
_part_dpo_(Wf_po, lbl='Final')

pl.legend(title='Exc to Exc', frameon=False, numpoints=1)

ax.set_xlim(0, 4)
ax.set_xticks([1, 2, 3])
ax.set_xticklabels(['0-30', '30-60', '60-90'], rotation=0)

ax.set_yticks([0, .2, .4, .6, .8, 1])
ax.set_yticklabels([0, .2, .4, .6, .8, 1])

ax.set_ylim(0, 1)


pl.xlabel('dpo range (deg)')
pl.ylabel('Average Weight (mV)')

pl.subplots_adjust(left=.05, right=.97, bottom=.15, top=.925, wspace=.45)

pl.savefig('Fig3B')

pl.show()

# ----------


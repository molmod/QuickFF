#! /usr/bin/env python

from molmod.periodic import periodic as pt
from molmod.units import *
from pytool.fit import fitpar
import numpy as np, cPickle, matplotlib.pyplot as pp, os

from lib.system import System
from lib.evaluators import ic_evaluator, energy_evaluator

root_model = '/home/louis/mil53/model/metals/covalent/est'
qunit = 'deg'
eunit = 'kjmol'
kunit = 'kjmol/rad**2'
#kunit = '%s/%s**2' %(eunit,qunit)


####  Setup  ####
#icname = 'dihed/H_PH.C_PH.C_PH.H_PH/cos-m2-0/dihed/K'
icname = 'angle/C_PH.C_PH.C_PC/harm/angle/K'
sysname = 'Al/linker'
setups = [['UncoupFree', None, None], ['UncoupFix', None, 3], ['SymFree', 'symmetric', None], ['SymFix', 'symmetric', 3]] #list of [name, coupling, free_depth]
if   icname.split('/')[0]=='bond' : amplitude_factor = 0.1
elif icname.split('/')[0]=='angle': amplitude_factor = 0.05
elif icname.split('/')[0]=='dihed': amplitude_factor = 0.01


####  A plot method  ####
def add_plot(ax, curves, title='', xunit='au', yunit='au'):
    for x, y, style, label in curves:
        ax.plot(x/parse_unit(xunit), y/parse_unit(yunit), style, label=label)
    ax.set_xlabel('IC [%s]' %xunit)
    ax.set_ylabel('Energy [%s]' %yunit)
    #ax.set_ylim([-50,50])
    ax.grid()
    ax.set_title(title)
    ax.legend(loc='lower left') 


####  Initialization  ####
pp.clf()
types = icname.split('/')[1]
with open('%s/data/%s/system.pp' %(root_model,sysname), 'r') as f: system = cPickle.load(f)


####  Calculations  ####
vs = {}
for isetup, (name, coupling, free_depth) in enumerate(setups):
    vs[name] = system.geometry_perturbation(icname, coupling=coupling, free_depth=free_depth)
fig, axs = pp.subplots(max([len(vs[name]) for name in vs]), len(setups))
if len(vs[name])==1:
    axs = np.array([axs])

for isetup, (name, coupling, free_depth) in enumerate(setups):
    for iv, v in enumerate(vs[name]):
        print 'Processing %s - %s: %s-%i' %(sysname, icname, name, iv)
        fn_xyz = '%s/data/%s/traj-%s-%s-%i.xyz' %(root_model, sysname, types, name, iv)
        amplitude = system.ics[icname][iv].value(system.sample['coordinates'])*amplitude_factor
        values = system.analyze_trajectory(
            v, evaluators=[ic_evaluator(icname, iv), energy_evaluator('totharm'), energy_evaluator('eiharm')], 
            fn_xyz=fn_xyz, amplitude=amplitude, steps=101
        )
        qs = np.array(values[0])
        totharm = np.array(values[1])
        eiharm = np.array(values[2])
        if coupling is not None: #if coupling is not None, the energy should be devided by the number of ics to get the energy per ic
            totharm /= len(vs[name])
            eiharm /= len(vs[name])
        pars = fitpar(qs, totharm, rcond=1e-6)
        k = 2*pars[0]
        q0 = -pars[1]/k
        pars_cov = fitpar(qs, totharm-eiharm, rcond=1e-6)
        k_cov = 2*pars_cov[0]
        q0_cov = -pars_cov[1]/k_cov
        covfit = pars_cov[0]*qs**2 + pars_cov[1]*qs + pars_cov[2]
        rmsd = np.sqrt( ((totharm - eiharm - covfit)**2).sum()/len(totharm) )
        title = "%s perturbation %s:\n" \
              + "---------------------------------------------\n" \
              + "Total:   k  = %4.0f    q0 = %4.2f\n" \
              + "Coval:   k  = %4.0f    q0 = %4.2f\n" \
              + "---------------------------------------------\n" \
              + "RMS(Total-Cov-Elec) = %.3e"
        title = title %(
            name, '-'.join(['%s[%i]' %(system.sample['ffatypes'][atindex],atindex) for atindex in system.ics[icname][iv].indexes]),
            k/parse_unit(kunit), q0/parse_unit(qunit), k_cov/parse_unit(kunit), q0_cov/parse_unit(qunit), rmsd/parse_unit(eunit)
        )
        curves = [
            [qs, totharm, 'b-', 'Total'], [qs, eiharm+covfit, 'b--', '_nolabel_'], [qs, eiharm, 'g-', 'Electrostatic'],
            [qs, totharm-eiharm, 'r-', 'Covalent'], [qs, covfit, 'r--', '_nolabel_'],
        ]
        add_plot(axs[iv, isetup], curves, title=title, xunit=qunit, yunit=eunit)

fig.set_size_inches([8*len(setups),8*len(vs[name])])
fig.tight_layout()
pp.savefig('%s/data/%s/energy-%s.pdf' %(root_model, sysname, types))
os.system('qevince %s/data/%s/energy-%s.pdf' %(root_model, sysname, types)) 

#!/usr/bin/env python

from lib.system import System, FFTable
from pytool.html import header as htmlheader, footer as htmlfooter, print_table
import os, sys, cPickle, numpy as np


def format_table(systems, fcs, icnames, tags):
    sorter = lambda icname: len(icname.split('/')[1].split('.'))
    kind_dic = {'bond': '--#b#lBond lengths', 'angle': '--#b#lBending angles', 'dihed': '--#b#lDihedral angles'}
    table = [['_ ']]
    for system in systems:
        metal = system.name.split('/')[0]
        if metal not in table[0]:
            table[0].append(metal)
    ickind = None
    for icname in sorted(icnames, key=sorter):
        row = []
        ictype = ' - '.join(icname.split('/')[1].split('.'))
        if '%s' in ictype: ictype = ictype %( ('M',)*ictype.count('%') )
        if ickind!=icname.split('/')[0]:
            ickind = icname.split('/')[0]
            table.append([kind_dic[ickind]]+[' ',]*len(table[0][1:]))
        row.append('#r'+ictype)
        for metal in table[0][1:]:
            cell = '<pre>'
            for system, fc in zip(systems, fcs):
                assert system.name==fc.name
                if system.name.split('/')[0]==metal:
                    cell += fc.html(icname, tags, suffix=' [%1s]' %(system.name.split('/')[1][0].upper()) )
            cell = cell.rstrip('<br>')
            cell += '</pre>'
            row.append(cell)
        table.append(row)
    return table

#Reading data
from init import icnames, sysnames, root_model
kind = sys.argv[1]
systems = []
ff = []
for name in sysnames:
    with open('%s/data/%s/system.pp' %(root_model, name), 'r') as f: systems.append(cPickle.load(f))
    with open('%s/data/%s/ff_%s.pp' %(root_model, name, kind), 'r') as f: ff.append(cPickle.load(f))

#Writing data to html
f = open('tot_%s.html' %kind, 'w')
print >> f, htmlheader
print >> f, '<p>'
print >> f, '<h1 style=\"text-align: left;\">FF parameters for the TOTAL interactions estimated from the DFT Hessian</h1>'
print_table(format_table(systems, ff, icnames, ['k']), caption='Force constants', f=f)
print >> f, '</p><br><p>'
print_table(format_table(systems, ff, icnames, ['q']), caption= 'Rest values', f=f)
print >> f, '</p>'
print >> f, htmlfooter
f.close()

f = open('cov_%s.html' %kind, 'w')
print >> f, htmlheader
print >> f, '<p>'
print >> f, '<h1 style=\"text-align: left;\">FF parameters for the COVALENT interactions estimated from the DFT Hessian</h1>'
print_table(format_table(systems, ff, icnames, ['k_cov']), caption='Force constants', f=f)
print >> f, '</p><br><p>'
print_table(format_table(systems, ff, icnames, ['q_cov']), caption= 'Rest values', f=f)
print >> f, '</p>'
print >> f, htmlfooter
f.close()
sys.exit()

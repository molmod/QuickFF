#! /usr/bin/env python

from molmod.units import *
import numpy as np

from tools import statistics


__all__ = ['FFArray', 'FFTable']


class FFArray(object):
    def __init__(self, data):
        if isinstance(data, np.ndarray):
            self.data = data
        else:
            self.data = np.array(data)
        self.mean, self.std, self.num = statistics(self.data)
    
    def __len__(self):
        return len(self.data)
    
    def _html(self, fmt, unit='au'):
        if self.num==0: return ''
        return fmt %(self.mean/parse_unit(unit), self.std/parse_unit(unit), self.num)



class FFTable(object):    
    def __init__(self, icnames):
        self.icnames = icnames
        self.k = {}
        self.q = {}
        self.units = {}
    
    @classmethod
    def from_ffit2(cls, model):
        print 'FFTAB  FFIT2: constructing FFTable from FFit2 model'
        icnames = []
        units = {}
        kdata = {}
        qdata = {}
        for rule in model.rules:
            for par in rule.pars:
                icname = '/'.join(par.prefix.split('/')[:2])
                kind = par.name[0].lower()
                if not icname in icnames: icnames.append(icname)
                unit = units.get(icname, {'k': None, 'q': None})
                unit[kind] = par.unit.replace('^', '**')
                units[icname] = unit
                if kind=='k':
                    kdata[icname] = np.array([par.value])
                elif kind=='q':
                    qdata[icname] = np.array([par.value])
                else:
                    raise ValueError('Invalid value for kind, recieved %s' %kind)
        fftab = FFTable(icnames)
        for icname in icnames:
            if icname.startswith('dihed/'):
                qdata[icname] = np.array([0.0])
                units[icname]['q'] = 'deg'
            fftab.add(icname, kdata[icname], qdata[icname], unit=units[icname])
        return fftab
    
    def add(self, icname, kdata, qdata, unit={'q': 'au', 'k': 'kjmol/au**2'}):
        assert icname in self.icnames
        assert icname not in self.k.keys()
        assert icname not in self.q.keys()
        self.k[icname] = FFArray(kdata)
        self.q[icname] = FFArray(qdata)
        self.units[icname] = unit
    
    def __getitem__(self, key, return_std=False):
        k = self.k[key].mean
        k_std = self.k[key].std
        q = None
        q_err = None
        if key in self.q.keys():
            q = self.q[key].mean
            q_std = self.q[key].std
        if not return_std:
            return k, q
        else:
            return k, k_std, q, q_std
    
    def print_screen(self):
        print 'FFTAB  PRINT: printing force field parameters to screen'
        print
        for icname in sorted(self.icnames):
            k, k_std, q, q_std = self.__getitem__(icname, return_std=True)
            k = u'% 8.2f \u00B1 % 8.2f %15s' %(k/parse_unit(self.units[icname]['k']), k_std/parse_unit(self.units[icname]['k']), self.units[icname]['k'])
            if q is None: q = 'None'
            else: q = u'% 8.3f \u00B1 % 8.3f %3s' %(q/parse_unit(self.units[icname]['q']), q_std/parse_unit(self.units[icname]['q']), self.units[icname]['q'])
            print '                %40s   K=%s     q=%s' %(icname, k, q)
        print

    def dump_pars_ffit2(self, fn):
        f = open(fn, 'w')
        print >> f, '# Parameters'
        print >> f, '# ------------------------------------------------------------------#------'
        print >> f, '# longname                                 unit               value # fx/fr'
        print >> f, '# ------------------------------------------------------------------#------'
        for icname in self.icnames:
            kind = icname.split('/')[0]
            atypes = icname.split('/')[1]
            k = self.k[icname].mean
            q0 = self.q[icname].mean
            if kind in ['bond', 'dist']:
                print >>f, '  %40s %12s % 12.6f #  free' %(
                    'bond/%s/harm/dist/K' %atypes,
                    'kjmol/A^2', k/(kjmol/angstrom**2)
                )
                print >> f, '  %40s %12s % 12.6f #  free' %(
                    'bond/%s/harm/dist/q0' %atypes,
                    'A', q0/angstrom
                )
            elif kind in ['bend', 'angle']:
                print >> f, '  %40s %12s % 12.6f #  free' %(
                    'angle/%s/harm/angle/K' %atypes,
                    'kjmol/rad^2', k/(kjmol/rad**2)
                )
                print >> f, '  %40s %12s % 12.6f #  free' %(
                    'angle/%s/harm/angle/q0' %atypes,
                    'deg', q0/deg
                )
            elif kind in ['dihedral', 'dihed', 'torsion']:
                print >> f, '  %40s %12s % 12.6f #  free' %(
                    'dihed/%s/cos-m2-0/dihed/K' %atypes,
                    'kjmol', k/kjmol
                )
        print >> f, '# ------------------------------------------------------------------#------'
        f.close()
        print 'FFTAB  DUMP : dumped parameters in FFit2 format to %s' %fn

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
        if self.std==None: self.std = 0.0
        return fmt %(self.mean/parse_unit(unit), self.std/parse_unit(unit), self.num)



class FFTable(object):    
    def __init__(self, icnames, units):
        self.units = units
        self.icnames = icnames
        self.k = {}
        self.q = {}
    
    def add(self, icname, kdata, qdata):
        assert icname in self.icnames
        assert icname not in self.k.keys()
        assert icname not in self.q.keys()
        self.k[icname] = FFArray(kdata)
        self.q[icname] = FFArray(qdata)
    
    def __getitem__(self, key):
        if not key in self.icnames: raise KeyError('%s is not a valid icname' %key)
        k = self.k[key].mean
        q = None
        if key in self.q.keys():
            q = self.q[key].mean
        return k, q
    
    def print_html(self, icname, attrs, suffix=None):
        """
            Returns strings containing the statistical information of all attributes 
            in <attrs> for the given icname. Allowed attributes: q, k
        """
        fmts  = {
            'q'     :   ('q<sub> </sub> = ' ,'% 8.3f &#177 %7.3f (%i)'),
            'k'     :   ('k<sub> </sub> = ' ,'% 8.1f &#177 %5.1f (%i)'),
        }
        if '%s' in icname: icname = icname %( (self.name.split('/')[0].upper(),)*icname.count('%') )
        result = ''
        for attr in attrs:
            if len(attrs)>1: result += fmts[attr][0]
            result += getattr(self, attr)[icname]._html(fmts[attr][1], unit=self.units[icname][attr])
            if len(result)>0:
                if suffix is not None: result += suffix
                result += '<br>'
        return result
    
    def print_screen(self):
        for icname in self.icnames:
            k, q = self[icname]
            k = '%.6f' %(k/parse_unit(self.units[key]['k']))
            if q is None: q = 'None'
            else: q = '%.3f' %(q/parse_unit(self.units[key]['q']))
            print '%s:    K = %.6f  q = %.6f' %(icname, k, q)

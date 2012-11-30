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
    
    def __getitem__(self, key, return_std=False):
        if not key in self.icnames: raise KeyError('%s is not a valid icname' %key)
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
        print ' FFTAB PRINT: printing force field parameters to screen'
        print
        for icname in self.icnames:
            k, k_std, q, q_std = self.__getitem__(icname, return_std=True)
            k = u'%7.2f \u00B1 %7.2f %15s' %(k/parse_unit(self.units[icname]['k']), k_std/parse_unit(self.units[icname]['k']), self.units[icname]['k'])
            if q is None: q = 'None'
            else: q = u'%7.3f \u00B1 %7.3f %3s' %(q/parse_unit(self.units[icname]['q']), q_std/parse_unit(self.units[icname]['q']), self.units[icname]['q'])
            print '                %20s   K=%s     q=%s' %(icname, k, q)

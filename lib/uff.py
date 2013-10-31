from molmod.units import angstrom, kjmol, kcalmol
from molmod.periodic import periodic as pt
import numpy as np, os

from context import context

__all__=['get_uff_sigmas_epsilons']

_fn  = context.get_fn('uff.prm')

def _read_uff_file():
    pars = {}
    _f = open(_fn, 'r')
    for line in _f.readlines():
        if line.startswith('#key'):
            keys = line.lstrip('#').split()[2:]
        elif line.startswith('param'):
            words = line.split()
            dic = {}
            for key, value in zip(keys, [float(x) for x in words[2:]]):
                dic[key] = value
            pars[words[1]] = dic
    _f.close()
    return pars

def _match_ffatypes(numbers, ufftypes):
    ufftype_matches = {}
    for number in numbers:
        matches = []
        symbol = pt[number].symbol
        for ufftype in ufftypes:
            symbol_uff = ufftype[:2].rstrip('_')
            if symbol_uff.lower() == symbol.lower():
                matches.append(ufftype)
        if len(matches)==0:
            raise ValueError('No matching UFF atom type found for %s' %symbol)
        else:
            ufftype_matches[number] = matches
    return ufftype_matches

def get_uff_sigmas_epsilons(numbers):
    '''
       A method to get the sigma and epsilon LJ-parameters from the
       Universal Force Field.

       **Arguments:**

       numbers
            A list of atomic numbers for which the LJ-parameters need
            to be returned.
    '''
    pars = _read_uff_file()
    ufftype_matches = _match_ffatypes(numbers, pars.keys())
    sigmas = np.zeros(len(numbers), float)
    epsilons = np.zeros(len(numbers), float)
    for i, number in enumerate(numbers):
        sarray = []
        earray = []
        for ufftype in ufftype_matches[number]:
            values = pars[ufftype]
            sarray.append(values['x1']*angstrom/2.0**(1.0/6.0))
            earray.append(values['D1']*kcalmol)
        sarray = np.array(sarray)
        earray = np.array(earray)
        #Multiple values found, take average
        if len(sarray)>1:
            sarray_mean = sarray.mean()
            sarray_std  = np.sqrt( ((sarray-sarray.mean())**2 ).sum()/(len(sarray)-1))
            earray_mean = earray.mean()
            earray_std  = np.sqrt( ((earray-earray.mean())**2 ).sum()/(len(earray)-1))
            if sarray_std!=0.0 or earray_std!=0.0:
                print 'WARNING: parameters for %s are not unique (sigma std = %.3f A  ,  epsilon std = %.6f kJ/mol)' % (
                    ffatype, sarray_std/angstrom, earray_std/kjmol
                )
                print '         Mean value was used: sigma = %.3f A  ,  epsilon = %.6f kJ/mol' %(
                    sarray_mean/angstrom, epsilons_mean/kjmol
                )
            sigmas[i] = sarray_mean
            epsilons[i] = earray_mean
        #Only one value found
        elif len(sarray)==1:
            sigmas[i] = sarray[0]
            epsilons[i] = earray[0]
        else:
            raise ValueError('No UFF vdW parameters found atomic number %' % number)
    return sigmas, epsilons

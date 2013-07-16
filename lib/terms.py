#! /usr/bin/env python

import numpy as np


__all__ = ['BaseTerm', 'HarmonicTerm', 'CosineTerm']

class BaseTerm(object):
    'An abstract class for valence force field terms'
    def __init__(self, ic, coords, k, q0, A):
        self.ic = ic
        self.coords = coords
        self.k = k
        self.q0 = q0
        self.A = A

    def calc_energy(self, **kwargs):
        'Method for calculating the energy'
        raise NotImplementedError

    def calc_gradient(self, **kwargs):
        'Method for calculating the energy gradient'
        raise NotImplementedError

    def calc_hessian(self, **kwargs):
        'Method for calculating the energy hessian'
        raise NotImplementedError

    def parse_kwargs(self, **kwargs):
        '''
            A method for interpreting keyword arguments given to a calc_*
            method. Missing arguments are compleneted with the default
            value stored as attributes (self.coords, self.k, ...)
        '''
        if 'coords' in kwargs.keys():
            coords = kwargs['coords']
        else:
            assert self.coords is not None
            coords = self.coords
        if 'k' in kwargs.keys():
            k = kwargs['k']
        else:
            assert self.k is not None
            k = self.k
        if 'q0' in kwargs.keys():
            q0 = kwargs['q0']
        else:
            assert self.q0 is not None
            q0 = self.q0
        if self.A is None:
            return coords, k, q0
        else:
            if 'A' in kwargs.keys():
                A = kwargs['A']
            else:
                A = self.A
            return coords, k, q0, A


class HarmonicTerm(BaseTerm):
    '''
        A harmonic term of the form

            0.5*k*(q-q0)^2
    '''
    def __init__(self, ic, coords, k, q0):
        BaseTerm.__init__(self, ic, coords, k, q0, None)

    def calc_energy(self, **kwargs):
        'Method for calculating the energy'
        coords, k, q0 = self.parse_kwargs(**kwargs)
        return 0.5*k*(self.ic.value(coords) - q0)**2

    def calc_gradient(self, **kwargs):
        'Method for calculating the energy gradient'
        coords, k, q0 = self.parse_kwargs(**kwargs)
        return k*(self.ic.value(coords) - q0)*self.ic.grad(coords)

    def calc_hessian(self, **kwargs):
        'Method for calculating the energy hessian'
        coords, k, q0 = self.parse_kwargs(**kwargs)
        q = self.ic.value(coords)
        g = self.ic.grad(coords)
        h = self.ic.hess(coords)
        return k*(np.outer(g, g) + (q-q0)*h)


class CosineTerm(BaseTerm):
    '''
        A single cosine term of the form

            0.5*k*[1-cos(A*(q-q0))]
    '''
    def calc_energy(self, **kwargs):
        'Method for calculating the energy'
        coords, k, q0, A = self.parse_kwargs(**kwargs)
        q = self.ic.value(coords)
        return 0.5*k*( 1.0 - np.cos(A*(q-q0)) )

    def calc_gradient(self, **kwargs):
        'Method for calculating the energy gradient'
        coords, k, q0, A = self.parse_kwargs(**kwargs)
        q = self.ic.value(coords)
        g = self.ic.grad(coords)
        return 0.5*k*A*np.sin(A*(q-q0))*g

    def calc_hessian(self, **kwargs):
        'Method for calculating the energy hessian'
        coords, k, q0, A = self.parse_kwargs(**kwargs)
        q = self.ic.value(coords)
        g = self.ic.grad(coords)
        h = self.ic.hess(coords)
        return 0.5*k*A**2*np.cos(A*(q-q0))*np.outer(g, g) +\
               0.5*k*A*np.sin(A*(q-q0))*h

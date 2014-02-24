# -*- coding: utf-8 -*-
#QuickFF is a code to quickly derive accurate force fields from ab initio input.
#Copyright (C) 2012 - 2014 Louis Vanduyfhuys <Louis.Vanduyfhuys@UGent.be>
#Steven Vandenbrande <Steven.Vandenbrande@UGent.be>, 
#Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center for Molecular Modeling
#(CMM), Ghent University, Ghent, Belgium; all rights reserved unless otherwise
#stated.
#
#This file is part of QuickFF.
#
#QuickFF is free software; you can redistribute it and/or
#modify it under the terms of the GNU General Public License
#as published by the Free Software Foundation; either version 3
#of the License, or (at your option) any later version.
#
#QuickFF is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, see <http://www.gnu.org/licenses/>
#
#--

__all__ = ['eval_ic', 'eval_energy']

def eval_ic(ic):
    '''
        Return an evaluator for the value of the given internal coordinate, i.e.
        return a method that calculates the value of the given ic for a certain
        atomic configuration.

        **Arguments**

        ic
            an instance of the class :class:`quickff.ic.IC`
    '''
    def evaluator(model, coords):
        return ic.value(coords)
    return evaluator


def eval_energy(part_name):
    '''
        Return an evaluator for the value of the `part_name` energy, i.e.
        return a method that calculates the `part_name` energy for a certain
        atomic configuration.

        **Arguments**

        part_name
            a string defining an attribute of a model instance that is needed
            to evaluate the correct part of the energy.
    '''
    def evaluator(model, coords):
        part_model = getattr(model, part_name)
        return part_model.calc_energy(coords)
    return evaluator

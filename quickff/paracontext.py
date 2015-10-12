# -*- coding: utf-8 -*-
# QuickFF is a code to quickly derive accurate force fields from ab initio input.
# Copyright (C) 2012 - 2015 Louis Vanduyfhuys <Louis.Vanduyfhuys@UGent.be>
# Steven Vandenbrande <Steven.Vandenbrande@UGent.be>,
# Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center for Molecular Modeling
# (CMM), Ghent University, Ghent, Belgium; all rights reserved unless otherwise
# stated.
#
# This file is part of QuickFF.
#
# QuickFF is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# QuickFF is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
#--

'''Convenience functions to enable using scoop.
'''

__all__ = ['ParaContext', 'paracontext']

class FakeFuture(object):
    def __init__(self, fun, *args, **kargs):
        self.args = args
        self.kargs = kargs
        self._result = fun(*args, **kargs)

    def result(self):
        return self._result


class ParaContext(object):
    def __init__(self):
        # initialize with serial version of map and submit
        self.use_stub()

    def use_stub(self):
        def my_map(fn, l, **kwargs):
            return [fn(i, **kwargs) for i in l]
        def my_wait_first(fs):
            return fs[:1], fs[1:]
        def debug_log(*args):
            with open('debug.log', 'a') as f:
                print >> f, ' '.join(str(arg) for arg in args)
            return 0
        self.map = my_map
        self.wait_first = my_wait_first
        self.submit = FakeFuture
        self.debug_log = debug_log

    def use_scoop(self):
        from scoop import futures
        from scoop import WORKER_NAME
        def my_map(*args, **kwargs):
            return list(futures.map(*args, **kwargs))
        def my_wait_first(fs):
            return futures.wait(fs, return_when=futures.FIRST_COMPLETED)
        def debug_log(*args):
            with open('debug_%s.log' % WORKER_NAME, 'a') as f:
                print >> f, ' '.join(str(arg) for arg in args)
            return 0
        self.map = my_map
        self.wait_first = my_wait_first
        self.submit = futures.submit
        self.debug_log = debug_log


paracontext = ParaContext()

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

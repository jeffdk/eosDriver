#from eosDriver import eosDriver


class barytropicEOS(object):

    # If youdon't understand **kwargs, go here:
    # http://stackoverflow.com/questions/3394835/args-and-kwargs
    def __init__(self, **kwargs):
        raise NotImplementedError("Subclasses should implement this!")

    def rho_bFromP(self, P):
        raise NotImplementedError("Subclasses should implement this!")

    def rho_bFromEps(self, eps):
        raise NotImplementedError("Subclasses should implement this!")

    def PfromRho_b(self, rho_b):
        raise NotImplementedError("Subclasses should implement this!")

    def PfromEps(self, eps):
        raise NotImplementedError("Subclasses should implement this!")

    def epsFromP(self, P):
        raise NotImplementedError("Subclasses should implement this!")

    def epsFromRho_b(self, rho_b):
        raise NotImplementedError("Subclasses should implement this!")


#!/usr/bin/python3
import numpy as np
from cosmoTransitions import pathDeformation as pd
import datetime
import subprocess

class Potential:
    """
    A sample potential. The depth of the absolute minimum is controlled with
    the parameters `fx` and `fy`.

    This potential has no physical significance whatsoever.
    """
    def __init__(self, c1=0.2434, c2=0.5233, c3=0.34234, c4=0.4747, c5=0.234808, c6=0.57023, c7=0.138912, c8=0.51723, c9=0.65889):
        self.params = c1,c2,c3,c4,c5,c6,c7,c8,c9

    def V(self, X):
        """
        This is a two-dimensional potential, so the input should be some
        array with a *last* axis of length 2. That is, the final index in the
        array should always be the one that specifies the field (in this case
        *x* or *y*). This is the convention that CosmoTransitions uses
        throughout.
        """
        x1,x2,x3,x4,x5,x6,x7,x8 = X[...,0], X[...,1], X[...,2], X[...,3], X[...,4], X[...,5], X[...,6], X[...,7]
        c1, c2, c3, c4, c5, c6, c7, c8, c9 = self.params
        r1 = c1*(x1-1.)**2 +c2*(x2-1.)**2 +c3*(x3-1.)**2 +c4*(x4-1.)**2 +c5*(x5-1.)**2 +c6*(x6-1.)**2 +c7*(x7-1.)**2 +c8*(x8-1.)**2
        r2 = x1**2 + x2**2 + x3**2 + x4**2 + x5**2 + x6**2 + x7**2 + x8**2
        return (r1-c9)*r2

    def dV(self, X):
        """
        The output of the gradient should have the same shape as the input.
        The last index specifies the direction of the gradient.
        """
        x1,x2,x3,x4,x5,x6,x7,x8 = X[...,0], X[...,1], X[...,2], X[...,3], X[...,4], X[...,5], X[...,6], X[...,7]
        c1, c2, c3, c4, c5, c6, c7, c8, c9 = self.params
        r1 = c1*(x1-1.)**2 +c2*(x2-1.)**2 +c3*(x3-1.)**2 +c4*(x4-1.)**2 +c5*(x5-1.)**2 +c6*(x6-1.)**2 +c7*(x7-1.)**2 +c8*(x8-1.)**2
        r2 = x1**2 + x2**2 + x3**2 + x4**2 + x5**2 + x6**2 + x7**2 + x8**2
        rval = np.empty_like(X)
        rval[...,0] = 2.*c1*(x1-1.)*r2 + 2.*x1*(r1-c8)
        rval[...,1] = 2.*c2*(x2-1.)*r2 + 2.*x2*(r1-c8)
        rval[...,2] = 2.*c3*(x3-1.)*r2 + 2.*x3*(r1-c8)
        rval[...,3] = 2.*c4*(x4-1.)*r2 + 2.*x4*(r1-c8)
        rval[...,4] = 2.*c5*(x5-1.)*r2 + 2.*x5*(r1-c8)
        rval[...,5] = 2.*c6*(x6-1.)*r2 + 2.*x6*(r1-c8)
        rval[...,6] = 2.*c7*(x7-1.)*r2 + 2.*x7*(r1-c8)
        rval[...,7] = 2.*c8*(x8-1.)*r2 + 2.*x8*(r1-c8)
        return rval

    def plotContour(self):
        nx = 100
        X = np.linspace(-.2,1.2,nx)[:,None] * np.ones((1,nx))
        Y = np.linspace(-.2,1.2,nx)[None,:] * np.ones((nx,1))
        XY = np.rollaxis(np.array([X,Y]), 0, 3)
        Z = self.V(XY)
        plt.contour(X,Y,Z, np.linspace(np.min(Z), np.max(Z)*.3, 200),
                    linewidths=0.5)
#############################


m = Potential()
print(datetime.datetime.today())
Y = pd.fullTunneling([[1.,1.,1.,1.,1.,1.,1.,1.],[0.,0.,0.,0.,0.,0.,0.,0.]], m.V, m.dV, tunneling_init_params={'alpha':2}, deformation_deform_params={'fRatioConv':0.02})
print(datetime.datetime.today())
print('S_E = ', Y.action)


import math
from scipy.special import erfinv
from random import random
import numpy as np
import awkward as ak
import numba

class CrystallBall:

    def __init__(self, m=0, s=1, a=10, n=10):
        self.pi = 3.14159
        self.sqrtPiOver2 = math.sqrt(self.pi/2.0)
        self.sqrt2 = math.sqrt(2.0)
        self.m = m
        self.s = s
        self.a = a
        self.n = n
        self.fa = abs(self.a)
        self.ex = math.exp(-self.fa * self.fa/2)
        self.A  = pow(self.n/self.fa, self.n) * self.ex
        self.C1 = self.n/self.fa/(self.n-1) * self.ex
        self.D1 = 2 * self.sqrtPiOver2 * math.erf(self.fa/self.sqrt2)

        self.B = self.n/self.fa - self.fa
        self.C = (self.D1 + 2 * self.C1)/self.C1
        self.D = (self.D1 + 2 * self.C1)/2

        self.N = 1.0/self.s/(self.D1 + 2 * self.C1)
        self.k = 1.0/(self.n - 1)

        self.NA = self.N * self.A
        self.Ns = self.N * self.s
        self.NC = self.Ns * self.C1
        self.F = 1 - self.fa * self.fa/self.n
        self.G = self.s * self.n/self.fa
        self.cdfMa = self.cdf(self.m - self.a * self.s)
        self.cdfPa = self.cdf(self.m + self.a * self.s)

    def cdf(self, x):
        d = (x - self.m)/self.s
        if(d < -self.a):
            if(self.F - self.s * d/self.G > 0):
                return self.NC/pow(self.F - self.s * d/self.G, self.n - 1)
            else:
                return self.NC
        if(d > self.a):
            if(self.F + self.s * d/self.G > 0):
                return self.NC * (self.C - pow(self.F + self.s * d/self.G, 1 - self.n))
            else:
                return self.NC * self.C
        return self.Ns * (self.D - self.sqrtPiOver2 * math.erf(-d/self.sqrt2))


    def invcdf(self, u):
        if(u < self.cdfMa):
            if self.NC/u > 0:
                return self.m + self.G * (self.F - pow(self.NC/u, self.k))
            else:
                return self.m + self.G * self.F
        if(u > self.cdfPa):
            if(self.C - u/self.NC > 0):
                return self.m - self.G * (self.F - pow(self.C - u/self.NC, -self.k))
            else:
                return self.m - self.G * self.F
        return self.m - self.sqrt2 * self.s * erfinv((self.D - u/self.Ns)/self.sqrtPiOver2)


@numba.njit
def drawFromCB(mean, sigma, n, alpha, builder):

    for i in range(len(mean)):
        
        x = random()

        m = mean[i]
        s = sigma[i]
        a = alpha[i]
        nn = n[i] 
        if (nn == 0.0 or a == 0.0): 
            builder.append(1.0)
        else:
            pi = 3.14159
            sqrtPiOver2 = math.sqrt(pi/2.0)
            sqrt2 = math.sqrt(2.0)
            
            fa = abs(a)
            ex = math.exp(-fa * fa/2)
            A  = pow(nn/fa, nn) * ex
            C1 = nn/fa/(nn-1) * ex
            D1 = 2 * sqrtPiOver2 * math.erf(fa/sqrt2)
    
            B = nn/fa - fa
            C = (D1 + 2 * C1)/C1
            D = (D1 + 2 * C1)/2
    
            N = 1.0/s/(D1 + 2 * C1)
            k = 1.0/(nn - 1)
    
            NA = N * A
            Ns = N * s
            NC = Ns * C1
            F = 1 - fa * fa/nn
            G = s * nn/fa
    
            
            d = (x - m)/s
            if(d < - a):
                if(F - s * d/G > 0):
                    result = NC/pow(F - s * d/G, nn - 1)
                else:
                    result = NC
            elif(d > a):
                if(F + s * d/G > 0):
                    result = NC * (C - pow(F + s * d/G, 1 - nn))
                else:
                    result = NC * C
            else: result =  Ns * (D - sqrtPiOver2 * math.erf(-d/sqrt2))
    
            builder.append(result)

    return builder
    

def get_rndm(eta, nL, cset):
    # obtain parameters from correctionlib

    mean = cset.get("cb_params").evaluate(abs(eta), nL, 0)
    sigma = cset.get("cb_params").evaluate(abs(eta), nL, 1)
    n = cset.get("cb_params").evaluate(abs(eta), nL, 2)
    alpha = cset.get("cb_params").evaluate(abs(eta), nL, 3)

    # get random number following the CB
    builder = ak.ArrayBuilder()
    counts = ak.num(mean)
    rndm = drawFromCB(ak.flatten(mean),ak.flatten(sigma), ak.flatten(n), ak.flatten(alpha), builder).snapshot()
    return ak.unflatten(rndm,counts)

    
def get_std(pt, eta, nL, cset):
    # obtain parameters from correctionlib

    param0 = cset.get("poly_params").evaluate(abs(eta), nL, 0)
    param1 = cset.get("poly_params").evaluate(abs(eta), nL, 1)
    param2 = cset.get("poly_params").evaluate(abs(eta), nL, 2)

    # calculate value and return max(0, val)
    sigma = param0 + param1 * pt + param2 * pt*pt
    counts = ak.num(sigma)
    sigma = ak.flatten(sigma)
    np.asarray(sigma)[sigma < 0] = 0 
    
    return ak.unflatten(sigma, counts)


def get_k(eta, var, cset):
    # obtain parameters from correctionlib
    k_data = cset.get("k_data").evaluate(abs(eta), var)
    k_mc = cset.get("k_mc").evaluate(abs(eta), var)

    # calculate residual smearing factor 
    # return 0 if smearing in MC already larger than in data
    counts = ak.num(k_data)
    k_data = ak.flatten(k_data)
    k_mc = ak.flatten(k_mc)
    k=np.zeros(ak.count(k_data))
    k[k_mc < k_data] = (k_data[k_mc < k_data]**2 - k_mc[k_mc < k_data]**2)**.5
    k = ak.from_numpy(k)
    return ak.unflatten(k, counts)


    
@numba.njit
def replaceNaNs(pt_corr, pt, builder):

    for i in range(len(pt_corr)):
        if math.isnan(pt_corr[i]):
            builder.append(pt_corr[i])
        else:
            builder.append(pt[i])

    return builder

def pt_resol(pt, eta, nL, var, cset):
    """"
    Function for the calculation of the resolution correction
    Input: 
    pt - muon transverse momentum
    eta - muon pseudorapidity
    nL - muon number of tracker layers
    var - variation (standard is "nom")
    cset - correctionlib object

    This function should only be applied to reco muons in MC!
    """
    k = get_k(eta, var, cset)
    rndm = get_rndm(eta, nL, cset)
    std = get_std(pt, eta, nL, cset)

    pt_corr = pt * (1 + k * std * rndm)

    counts = ak.num(pt_corr)
    builder = ak.ArrayBuilder()
    pt_corr = replaceNaNs(ak.flatten(pt_corr), ak.flatten(pt),builder).snapshot()
    pt_corr = ak.unflatten(pt_corr, counts)
    
    return pt_corr


def pt_scale(is_data, pt, eta, phi, charge, var, cset):
    """
    Function for the calculation of the scale correction
    Input:
    is_data - flag that is True if dealing with data and False if MC
    pt - muon transverse momentum
    eta - muon pseudorapidity
    phi - muon angle
    charge - muon charge
    var - variation (standard is "nom")
    cset - correctionlib object
    
    This function should be applied to reco muons in data and MC
    """
    if is_data:
        dtmc = "data"
    else:
        dtmc = "mc"

    a = cset.get("a_"+dtmc).evaluate(eta, phi, var)
    m = cset.get("m_"+dtmc).evaluate(eta, phi, var)

    return 1. / (m/pt + charge * a)

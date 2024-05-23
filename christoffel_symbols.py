import numpy as np
import sympy as smp
from itertools import product

class ChristoffelSymbols():
    
  def get_symbols_from_metrics(self,g: smp.MutableDenseMatrix,x: np.ndarray) -> np.ndarray:
    """This function compute the christoffel symbols of a given metric in a given coordinate system

    Args:
        g (smp.MutableDenseMatrix): Metric tensor of the space-time
        x (np.ndarray): coordinate vector

    Returns:
        np.ndarray: np.ndarray of shape = (x.shape[0],x.shape[0],x.shape[0])
    """
  
    self.x = x   
    self.g= g
    g_inv = self.g.inv()    


    dim = self.g.shape[0]

    gamma_mu = np.zeros(shape = dim, dtype = 'object')
    for mu in range(dim):
      gamma_d = np.zeros((dim,dim), dtype='object')
      for (alpha,beta,l) in product(range(dim),repeat = 3):
        gamma_d[alpha,beta]+= smp.Rational(1/2)*(g_inv[mu,l]*(smp.diff(g[l,alpha],x[beta])+
                                                smp.diff(g[l,beta],x[alpha])-
                                                smp.diff(g[alpha,beta],x[l])
                                                ))
      gamma_mu[mu] = smp.simplify(gamma_d)

    return gamma_mu
 
import numpy as np
from itertools import product
import sympy as smp
from christoffel_symbols import ChristoffelSymbols


class Ricci(ChristoffelSymbols):
    def __init__(self):
        super().__init__()

    def RicciTensor(self,g: smp.MutableDenseMatrix,x:np.ndarray)  -> np.ndarray:
        """This method computes the Ricci tensor of a generic space-time provided a metric g and a coordinate system x_mu

        Args:
            g (smp.MutableDenseMatrix): Metric Tensor of the space-time
            x (np.ndarray): coordinates of the space-time

        Returns:
            np.ndarray: The Ricci Tensor of the space time
        """
    
        self.g = g
        self.x = x
        christoffel = self.get_symbols_from_metrics(self.g, self.x)
        dim = self.g.shape[0]
        derivatives = np.zeros(shape = (dim,dim),dtype = 'object')
        for (i,j) in product(range(dim),repeat =2):
            for k in range(dim):
                derivatives[i,j] += smp.diff(christoffel[k][i][j],x[k]) - smp.diff(christoffel[k][i][k],x[j])
        squares = np.zeros(shape = (dim,dim),dtype = 'object')
        for (i,j) in product(range(dim),repeat =2):
            for (k,m) in product(range(dim),repeat =2):
                squares[i,j] += christoffel[k][i][j]*christoffel[m][k][m]- christoffel[k][i][m]*christoffel[m][j][k]

        ricci_tensor = derivatives + squares
        return ricci_tensor

    
    def RicciScalar(self,g: smp.MutableDenseMatrix,x: np.ndarray) ->smp.Mul:
        """This method computes the Ricci Scalar Curvature of a generic space-time sith metric g and a coordinate system x_mu

        Args:
            g (smp.MutableDenseMatrix): metric tensor of the space-time
            x (np.ndarray): coordinate system of the space-time

        Returns:
            smp.Mul: Ricci Scalar of the Curvature
        """
        self.x = x
        self.g = g
        g_inv = self.g.inv()
        self.R = self.RicciTensor(self.g, self.x)
        tr_R = 0
        for (mu,nu) in product(range(self.g.shape[0]),repeat=2):
            tr_R += g_inv[mu,nu]*self.R[mu,nu]
        

        return tr_R
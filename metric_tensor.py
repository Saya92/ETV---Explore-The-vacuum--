import numpy as np
import sympy as smp
from itertools import product
from numpy.linalg import inv
from IPython.display import display
from sympy import det,sqrt,diag
from time import sleep

class DefineMetricTensor:
    def __init__(self,x_mu):

        self.x_mu = x_mu

    def __repr__(self) -> str:
        return f'MetricTensor object whoose coordinates are {self.x_mu}'

    def generate_covariant_metric_tensor(self) -> smp.MutableDenseMatrix:
        """ Given a coordinate vector, This function lets the user to create a metric tensor manually

        Raises:
            ValueError: Raises a ValuError if elements does not contain function of the coordinate vector. For Example, if the original coordinates
            were [tau,r] and the user types [tai,r] the function will generate the ValueError. Numbers are allowed 

        Returns:
            smp.MutableDenseMatrix: _description_
        """

        
        self.g = smp.Matrix.zeros(self.x_mu.shape[0],self.x_mu.shape[0], real = True, positive = True)
        for (mu,nu) in product(range(len(self.x_mu)),repeat = 2):
            if mu >= nu:
                self.g[mu,nu] = input("")           
                g_sym = smp.symbols(f"g_{self.x_mu[mu]}_{self.x_mu[nu]}")
                display(smp.Eq(g_sym,self.g[mu,nu]))
                print('\n Are you sure you inserted a function of the right coordinate? [Yes/No]')
                answer = input('')
                if answer == 'Yes':
                    g_T = self.g.T
                    diagonal = np.zeros(shape = (self.x_mu.shape[0],self.x_mu.shape[0]), dtype = 'object')
                    np.fill_diagonal(diagonal, np.diag(self.g))
                    self.g = self.g + g_T - diagonal
                else:
                    print('='*10)
                    print('\n You Typed The Wrong Coordinates: Exiting The Program in 5 seconds')
                    sleep(10.0)
                    break

            else:
                self.g[mu,nu] = 0
        g_T = self.g.T
        diagonal = np.zeros(shape = (self.x_mu.shape[0],self.x_mu.shape[0]), dtype = 'object')
        np.fill_diagonal(diagonal, np.diag(self.g))
        self.g = self.g + g_T - diagonal
        return self.g
    
    def generate_contravariant_metric_tensor(self,g: smp.MutableDenseMatrix) -> smp.MutableDenseMatrix:
        """This function computes the contravariant components of the metric tensor.

        Args:
            g (smp.MutableDenseMatrix): covariant component of the metric tensor

        Raises:
            ValueError: Raise a ValueError if the determinant of the metric tensor is zero.

        Returns:
            smp.MutableDenseMatrix: returns the contravariant components of the metric tensor.
        """
        self.g = g

        self.g_contravariant = smp.Matrix.zeros(shape = (self.g.shape[0],self.g.shape[0]), real = True, positive = True)
        
        if det(self.g)!=0:
            self.g_contravariant = self.g.inv()
        else: 
            raise ValueError('Metric tensor is not invertible')
        
        return self.g_contravariant
    

    
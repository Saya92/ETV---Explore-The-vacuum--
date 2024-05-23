import sympy as smp
from sympy import diag
from sympy import sin
import numpy as np


class Schwartzshild():

    def __init__(self):
        pass
        
    

    def coordinates(self) -> np.ndarray:
        """This method returns the Schwartzschild polar coordinates

        Returns:
            np.ndarray: Schwartzschild polar coordinates
        """
        t,r,theta,phi = smp.symbols('t r theta phi', real = True,positive = True)
        self.t,self.r,self.theta,self.phi = t,r,theta,phi
        coordinates = np.array([self.t,self.r,self.theta,self.phi], dtype='object')
        return coordinates


    def scwhartzschild_metric(self, geometrized: bool = True, signature: str = 'plus_minus', general:bool = False) -> smp.MutableDenseMatrix:
        """This method computes the covariant metric tensor of the Schwartzschild space-time

        Args:
            geometrized (bool, optional): When set to True dimensionful constant like G or c are set to one. Defaults to True.
            signature (str, optional): String parameter to decide the signature of the metric. Defaults to 'plus_minus'.
            general (bool, optional): When set to true the metric will have a more general form, with a non constant mass function and an exponential radial function 
            multiplying the temporal coordinate. Defaults to False.

        Returns:
            smp.MutableDenseMatrix: The covariant metric tensor of the Schwartzschild space-time
        """
        x = self.coordinates()
        self.M , self.c, self.G  = smp.symbols('M c G', real = True, positive = True)
        #self.t,self.r,self.theta, self.phi = smp.symbols('t r theta phi')
        m = smp.Function('M')('r')
        delta = smp.exp(2*smp.Function('delta', positive = True, real = True)('r'))
        

        
        if general == True:
           self.f_r = (1-(2*m*self.G)/(x[1]*self.c**2))
        else:
           self.f_r = (1-(2*self.M*self.G)/(x[1]*self.c**2))
        


        if geometrized == True:
            self.f_r = self.f_r.subs({self.c:1,self.G:1})
        
        self.g = smp.Matrix.zeros(s_t,s_t)
        
        if signature:
            if general == False:
                self.g = diag(self.f_r,-self.f_r**(-1),-x[1]**2,-x[1]**2*sin(x[2])**2)
            else:
                self.g = diag(self.f_r*delta,-self.f_r**(-1),-x[1]**2,-x[1]*sin(x[2])**2)

        
        elif signature == 'minus_plus':
            if general == False:
                self.g = diag(-self.f_r,self.f_r**(-1),x[1]**2,x[1]**2*sin(x[2])**2)
            else:
                self.g = diag(-self.f_r*delta,self.f_r**(-1),x[1]**2,x[1]**2*sin(x[2])**2)


        
        return self.g
        

    



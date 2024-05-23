import sympy as smp
import numpy as np
from sympy import det,sqrt
from itertools import product
from IPython.display import display

def to_numeric_symbol(A: smp.MutableDenseMatrix,sub:dict) -> smp.MutableDenseMatrix:
    """_summary_

    Args:
        A (smp.MutableDenseMatrix): _description_
        sub (dict): _description_

    Returns:
        smp.MutableDenseMatrix: _description_
    """
    num_A = A.subs(sub)
    return num_A

def to_numpy(A: smp.MutableDenseMatrix) -> np.ndarray:
    """Function which converts sympy MutableDense Matrix to np.ndarray

    Args:
        A (smp.MutableDenseMatrix): A generic Matrix

    Returns:
        np.ndarray: The input array as a numpy ndarray
    """
    return np.array(A)

def determinant(g: smp.MutableDenseMatrix) -> smp.Mul:
    """This function simply computes the determinant of the metric tensor
    Args:
        g (smp.MutableDenseMatrix): _description_

    Returns:
        smp.Mul: Determinant of the metric tensor
    """

    g = g
    d_g = det(g)
    return d_g

def sqrt_det(g):
    g = g
    return sqrt(-g)

def contravariant_metric_tensor(g: smp.MutableDenseMatrix) -> smp.MutableDenseMatrix:
    """Given the covariante metric tensor, 
       this function computes the contravariant metric tensor 
       by inverting the metric tensor itself

    Args:
        g (smp.MutableDenseMatrix): Metric tensor of the space time

    Raises:
        ValueError: If the determinant of the metric tensor is 0, the function will raise a ValueError due to the non invertibility of g
    Returns:
        smp.MutableDenseMatrix: Contravariant Metric Tensor
    """
    if determinant(g)==0:
        raise ValueError('The metric tensor is not invertible')
    g_inv = g.inv()
    return g_inv

def trace(g:smp.MutableDenseMatrix,t:smp.MutableDenseMatrix) -> smp.Mul:
    """This function computes the trace of type-2 tensor as the contraction of the (contravariant) metric tensor and the tensor itself

    Args:
        g (smp.MutableDenseMatrix): Metric tensor of the space-time
        t (smp.MutableDenseMatrix): type 2 Tensor we want to compute the trace

    Raises:
        ValueError: Raises a ValueError exception whwnever the shape of the tensor t is different from the shape of the metric tensor

    Returns:
        smp.Mul: The trace of the tensor t
    """

    if t.shape!=g.shape:
        raise ValueError('Input tensor and metric tensor must have the same dimensions')
    g_inv = contravariant_metric_tensor(g)
    
    trace = 0
    for (mu,nu) in product(range(g.shape[0]),repeat=2):
        trace += g_inv[mu,nu]*t[mu,nu]
    
    return trace

def generate_contravariant_type_two_tensor(g: smp.MutableDenseMatrix,t: smp.MutableDenseMatrix) -> smp.MutableDenseMatrix:
    """Given a covariant type 2 tensor, this function generate the contravariant component of the tensor itself

    Args:
        g (smp.MutableDenseMatrix): Covariant component of the Metric Tensor 
        t (smp.MutableDenseMatrix): Covariant component of a tensor t
    Returns:
        smp.MutableDenseMatrix: returns the controvariant components of a tensor t
    """
    g_inv = contravariant_metric_tensor(g)
    t_contravariant = smp.Matrix.zeros(t.shape[0])
    for mu,nu in product(range(g.shape[0]),repeat =2):
         for alpha,beta in product(range(g.shape[0]), repeat = 2):
            t_contravariant[mu,nu] += g_inv[alpha,mu]*g_inv[beta,nu]*t[alpha,beta]
         
    return t_contravariant

def display_type_two_tensor(x_mu: np.ndarray,t: smp.MutableDenseMatrix,covariant: bool = True) -> None:
    """This function allows to display type two tensor in a latex-fashion way

    Args:
        x_mu (np.ndarray): space-time indeces
        t (smp.MutableDenseMatrix): Arbitrary Two Type tensor to display
        covariant (bool, optional): When set to True, indeces are covariant (upper). Defaults to True.
    """


    t_contravariant = generate_contravariant_type_two_tensor(t)
    symbols_t = smp.Matrix.zeros(x_mu.shape[0],x_mu.shape[0], real = True, positive = True)
    for (mu,nu) in product(range(x_mu.shape[0]),repeat = 2):
        if covariant:
                symbols_t[mu,nu] = smp.symbols(f"T_{x_mu[mu]}_{x_mu[nu]}")
        else:
                symbols_t[mu,nu] = smp.symbols(f"T^{x_mu[mu]}^{x_mu[nu]}")

        if covariant:
            display(smp.Eq(symbols_t,t))
        else:
            display(smp.Eq(symbols_t,t_contravariant))

def generate_contravariant_vector(g: smp.MutableDenseMatrix,v:smp.MutableDenseNDimArray) -> smp.MutableDenseNDimArray:
    """This function computes the contravariant component of a 4-vector v

    Args:
        g (smp.MutableDenseMatrix): covariant component of the metric tensor
        v (smp.MutableDenseNDimArray): covariant component of a 4-vector v
    Returns:
        smp.MutableDenseNDimArray: contravariant components of a 4-vector v
    """
    g_inv = g.inv()
    v_contravariant = g_inv@v
    return v_contravariant.T

def display_vector(g: smp.MutableDenseMatrix,v: smp.MutableDenseNDimArray,covariant:bool = True) -> None:
    """This function allows to display a 4 vector in a latex-fashion way

    Args:
        g (smp.MutableDenseMatrix): Metric tensor of the space time
        v (smp.MutableDenseNDimArray): Covariant components of the 4-vector to display
        covariant (bool, optional): Boolean value to decide to display covariant or contravariant components. Defaults to True.
    """
    v_contravariant = generate_contravariant_vector(g,v)
    symbols_v = smp.Array.zeros(v.shape[0], real = True, positive = True)
    for mu in v.shape[0]:
        if covariant:
                symbols_v[mu]= smp.symbols(f"v_{mu}")
        else:
                symbols_v[mu] = smp.symbols(f"v^{mu}")

        if covariant:
            display(smp.Eq(symbols_v,v))
        else:
            display(smp.Eq(symbols_v,v_contravariant))


     
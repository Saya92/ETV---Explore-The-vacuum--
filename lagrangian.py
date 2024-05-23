import sympy as smp

def gravity_lagrangian(R: smp.Mul,g: smp.MutableDenseMatrix) -> smp.Mul:
    """This function Computes the gravitational Lagrangian provided the space-time metric tensor and the associated Ricci Curvature


    Args:
        R (smp.Mul): Ricci scalar curvature of the space-time 
        g (smp.MutableDenseMatrix): metric tensor of the space time

    Returns:
        (smp.Mul): Gravitational Lagrangian
    """
 

    
    det_g = smp.det(g)
    L = R*smp.sqrt(-det_g)
    return L
        

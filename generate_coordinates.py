import numpy as np
import sympy as smp
from IPython.display import display

def get_coordinates(st_dim:int ) -> np.ndarray :
  """ This function let a user to insert manually the coordinate of a space-time, given the space-time dimension

  Args:
      st_dim (int): the dimension of the space-time

  Raises:
      ValueError: Raises a ValueError if st_dim is less-equal than 0

  Returns:
      np.ndarray: returns a np.ndarray containing the coordinates
  """  """
  """


  coordinates = []

  if st_dim <= 0:
    raise ValueError('The space-time dimension must greater than zero')

  while len(coordinates) <= st_dim:
    coordinates.append(input(f" Type the coordinate "))

  x_mu = coordinates
  return np.array(x_mu, dtype='object')

def display_coordinate(x_mu: np.ndarray, covariant:bool = True) -> None:
  """This function allows to display in a latex fashion style the coordinate vector

  Args:
      x_mu (np.ndarray): coordinate vector
      covariant (bool, optional): When set to true, indeces are displayed lower. Defaults to True.

  Raises:
      ValueError: Raises a Value Error if the the shape of x_mu is less-equal than 0
  """


  if x_mu.shape[0] <= 0:
    raise ValueError('Enter a valid coordinate vector: Its dimension must be greater than zero')
  dim = 0
  while dim < x_mu.shape[0]:
    if covariant:
        display(smp.Eq(smp.symbols(f"x_{dim}"),smp.symbols(x_mu[dim])))
        dim+=1
    else:
        display(smp.Eq(smp.symbols(f"x^{dim}"),smp.symbols(x_mu[dim])))
        dim+=1



======================
Hamiltonian tutorial
======================

The idea of this tutorial is ilustrate how to create a matrix represation of a hamiltonian using this package.
As it's expected, the hamiltonian that we consider are the one of spin models (the operators are the product of Pauli matrices of a given values of spin).


The Model
---------
Before we start using the package let introduce a example model (we will use it in all the tutorial). Let consider a triangle, this one is modeled as follow.

.. math::
  \mathcal{H} = J_1S_0S_1 + J_2S_1S_2 + J_3S_0S_3

The operator is defined as :math:`S_{i} = (\sigma_{x,i}, \sigma_{y,i}, \sigma_{z,i})`, where the :math:`\sigma_{i,j}` is the pauli matrices of a given spin in the site :math:`i`.
The product of this operators is a dot product, and the product in the resulting sum is a multiplication between the pauli matrices on different site. For example:

.. math::
   \sigma_{x,0}\sigma_{x,1} = (\sigma_{x,0} \otimes \mathbf{I} \otimes \mathbf{I})(\mathbf{I} \otimes \sigma_{x,1} \otimes \mathbf{I})

Its a simple matrix multiplication, but the pauli matrices are in different places.


Computing the paulices matrices of any spin
--------------------------------------------
The next point is how we can create the matrices of any given spin...



Important elements
------------------
Now that we know how some to compute some things, let list the important elements. The first one, we need to have a enumeration of the sites of the system to 
indicate the value of the spin of that site, so, we need to create a list with the values of the spin. The second one is a list of tuples that indicate 
which site is conected with other one (it's important that the first index is lower or equal that the second one). And finally, a list with the values of the exchange. 

Taking the :math:`\mathcal{H}` model, the three element should look as [0.5, 0.5, 0.5], [(0,1), (1,2), (0,2)] and [:math:`J_1`, :math:`J_2`, :math:`J_3`].

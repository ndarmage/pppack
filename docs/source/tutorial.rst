Tutorial
========

This section shows how to create polynomial-like approximants to multi-dimensional real-valued functions. A few reminders to the theory follow hereafter; complete information about the methods is available in De Boor's literature :cite:`deBoor1978prguide,deBoor1979TOMS,deBoor2001prguide`.

Define a function
-----------------

PPPACK objects allow to represent a multivariate real-valued function belonging to a multi-dimensional Cartesian space :math:`F^{(d)}(V^{(d)})` as :math:`f^{(d)}: V^{(d)} \rightarrow \mathbb{R}`, with :math:`V^{(d)} = V_1 \times \ldots \times V_d`, :math:`V_i=[a_i,b_i] \in \mathbb{R},\, \forall i`.

:math:`F^{(d)}` is constructed by tensor products with :math:`d` one-dimensional polynomial-like function spaces, :math:`F^{(d)} = F_1(V_1) \times \ldots \times F_d(V_d)`, where :math:`F_i=F_i^{(1)}`. The superscript ``(d)`` drops to access the class names of the objects, see the section :ref:`class-definitions`, where for instance :math:`F^{(d)}` becomes simply *Fd*.

Two types of one-dimensional function spaces :math:`F(=F^{(1)})` are available:

* :py:class:`~pppack.PP`, the piece-wise polynomial space of (polynomial) order :math:`k`, breakpoints :math:`\vec{\xi}=[\xi_1,\ldots,\xi_{l+1}]` and continuities of derivativies :math:`\vec{\nu}=[\nu_2,\ldots,\nu_{l}]`, :math:`\mathbb{P}_{k,\vec{\xi}(,\vec{\nu})}`, where :math:`l` is the number of polynomial pieces;
* :py:class:`~pppack.S`, the space constructed by B-splines of (polynomial) order :math:`k` and knot sequence :math:`\vec{t}`, :math:`\$_{k,\vec{t}}`.

The two spaces are equivalent if the requirements of the theorem of Curry and Schoenberg (1966) are fulfilled, see :py:meth:`~pppack.PP.PP2S`. Simple polynomial spaces can be obtained as well with :math:`l=1`, being simply :math:`\vec{\xi} = [a,b]`.

.. In this case, the data points :math:`\tau` are expected in :math:`[a,b]`.

Provided a base of functions on each dimension, i.e. :math:`F_j = \operatorname{span} ( \{\psi_{i_j},\,i_j=1,\ldots, I_j\} ),\; \forall j`, it is possible to write the multivariate function :math:`f` as:

.. math::
   f^{(d)}&= \sum_{i_1} C_{i_1} \psi_{i_1}(x_{i_1}) \otimes
                                             \ldots \otimes
             \sum_{i_d} C_{i_d} \psi_{i_d}(x_{i_d}) \\
          &= (f_{1} \otimes \ldots \otimes f_{d})(x_{i_1},\ldots,x_{i_d}) \\
          &= \sum_{\vec{i}} C_{\vec{i}} \psi_{\vec{i}} (\vec{x}),

with :math:`\vec{i}=[i_1, \ldots, i_d]` and :math:`\psi_{\vec{i}} = \prod_j{ \psi_{i_j} }`.

The coefficients :math:`C_{\vec{i}}`, instead of the heavier notation :math:`C_{i_1}C_{i_2} \ldots C_{i_d}`, are stored as *d*-dimensional FORTRAN array, see :py:class:`~pppack.fd`.

Possible one-dimensional (polynomial) functions :math:`f \in F` are:

* :py:class:`~pppack.pp`:math:`\;\in \mathbb{P}_{k,\vec{\xi}}`, a piece-wise polynomial function;
* :py:class:`~pppack.s`:math:`\;\in \$_{k,\vec{t}}`, a function constructed by B-splines.

We use the convention that functions are in lower case letters, whereas function spaces follow in upper case.

The tensor product operator is implemented by overloading the special method ``mult`` to build multi-dimensional entities.

:Example:

>>> from pppack import *
>>> k,t=2,[1., 1., 2., 3., 4., 4.]
>>> S1 = S(k,t)
>>> xi = [1., 4.] 
>>> P1 = PP(k,xi)
>>> F1 = S1 * P1
>>> print(F1)
F^(2) = F0 x F1
 with,
F0 = S_{ k=2, t }, with
 t = [1.000,1.000,2.000,3.000,4.000,4.000]
F1 = PP_{ k=2, xi, nu }, with
 xi = [1.000,4.000],  nu = [-,,-]
>>> s0 = s.inS(S1) # = S(k,t) = s(k,t)
>>> print(s0)
s \in S(k=2,t=[1.0, 1.0, 2.0, 3.0, 4.0, 4.0])
coef (to be determined)
---

In addition to inheriting the ``init`` method from their space, functions have also another class method for clearer usage (see the example above): :py:meth:`.s.inS`, :py:meth:`.pp.inPP` and :py:meth:`.fd.inFd`.

It is possible to assign a name to any object at ``init``; default name is the class name plus an integer counter, increasing internally as far as new instances are created. The creation of multi-dimensional spaces with muliple occurrencies of the same elements is not allowed.


Compute the coefficients
------------------------

The computation of the fitting coefficients requires the introduction of new linear functionals, :math:`\lambda_{\vec{j}}`, to apply to a target function :math:`g` :cite:`deBoor1978prguide`. As well, the functionals apply also to the approximant :math:`f` of :math:`g`, defining the kind of approximation:

.. math:: \lambda_{\vec{j}}g(\vec{x}) = \lambda_{\vec{j}}f(\vec{x}) = \sum_{\vec{i}} C_{\vec{i}} \lambda_{\vec{j}}\psi_{\vec{i}}(\vec{x}),

by linearity. Interpolation requires as many conditions as unknown coefficients to determine. The functionals can be integral or derivative conditions, or divided differences. The last form covers osculatory interpolation too:

.. math:: \lambda_{\vec{j}}g = [x_{j_1}(,*)] \otimes \ldots \otimes [x_{j_d}(,*)]g.

:math:`(,*)` means optional additional points for the directional divided differences. Simple evaluations of :math:`g` at the data points :math:`x_{j_i}`, now called :math:`\tau_{j_i}` as in the references, are

.. math:: \lambda_{\vec{j}}g = [x_{j_1}] \otimes \ldots \otimes [x_{j_d}]g = g(\tau_{j_1}, \ldots, \tau_{j_d}).

Other kinds of functionals are deferred to later versions.

The inversion of the Gramian matrices building :math:`\lambda_{\vec{j}}\psi_{\vec{i}}` and their inner product with :math:`\lambda_{\vec{j}}g` (index by index) yields the coefficients :math:`C_{\vec{i}}`. About the former, inversion is treated separately on the different dimensions: :math:`\lambda_{\vec{j}}\psi_{\vec{i}} = \lambda_{j_1}\psi_{i_1} \otimes \ldots \otimes \lambda_{j_d}\psi_{i_d}`, see :py:meth:`.s.cmpcoef` and :py:meth:`.pp.cmpcoef`. The following inner products are implemented according to De Boor's intuition in :cite:`deBoor1979TOMS`, see :py:meth:`.fd.cmpcoef`, in order to achieve high computational performances.


Evaluate a function
-------------------

Once the coefficients are known, function evaluation is accomplished by the ``call`` methods, for all classes of functions. Again, De Boor's intuition allows to handle smartly the many inner products on successive indices, see :py:meth:`.fd.__call__`. One-dimensional functions call back the low-level optimized routines of the F90 extension library.

:Example:

>>> import numpy as np
>>> g2= lambda x,y: np.sin(x)/(x+.1)+x*np.exp(y)+.25*y
>>> P1.name="linear"
>>> V2=Fd(S1,P1, name="2D")
>>> h2=fd.inFd(V2)
>>> tau=[]
>>> tau.append( S1.cmptau() )
>>> tau.append([1., 4.])
>>> Ni=[len(t) for t in tau]
>>> N=np.prod(Ni)
>>> #
>>> # Use FORTRAN order for continuous elements in memory
>>> # when creating arrays. Double precision arithmetics
>>> # is used as default in the fortran libraries. If simple
>>> # precision is requested, specify dtype=np.float32. However,
>>> # this will cause copies of data due to element type
>>> # mismatch. 
>>> gtau=np.zeros((N,), order='F')
>>> for x,i in zip(itt.product(*tau), \
...            itt.product(*[range(nx) for nx in Ni])):
...         gtau[getidx(i,Ni)]=g2(*x)
>>> h2.cmpcoef(tau,gtau)
>>> # verify interpolation at data sites
>>> for x in itt.product(*tau):
...     fx,gx=h2(x),g2(*x)
...     print("%+13.6e %+13.6e %+7.2f"%(fx,gx,(1.-fx/gx)*100.))
+3.733255e+00 +3.733255e+00   +0.00
+5.636312e+01 +5.636312e+01   +0.00
+6.119562e+00 +6.119562e+00   +0.00
+1.106293e+02 +1.106293e+02   +0.00
+8.450368e+00 +8.450368e+00   +0.00
+1.648400e+02 +1.648400e+02   +0.00
+1.093854e+01 +1.093854e+01   +0.00
+2.192080e+02 +2.192080e+02   +0.00
>>> print(h2)
fd \in Fd(S(k=2,t=[1.0, 1.0, 2.0, 3.0, 4.0, 4.0]), PP(k=2,xi=[1.0, 4.0],nu=[]))
coef = 
[  3.73325545   6.11956243   8.45036807  10.93854134  17.5432894
  34.8365788   52.1298682   69.42315761]
 tau = [array([ 1.,  2.,  3.,  4.]), [1.0, 4.0]]
gtau = [   3.73325545    6.11956243    8.45036807   10.93854134   56.36312366
  110.62929884  164.83997268  219.20801416]
---
>>> print("%13.6e"%h2([1.89,3.05]))
 7.337242e+01


Array management
----------------

FORTRAN passes conventionally dummy arguments by reference, that is by passing a pointer to the data. This avoids internal copies and intermediate conversion steps before addressing data, so providing efficient memory management.

This feature is largely used in this package to exchange the multi-dimensional arrays: arrays may appear with indices of different sizes w.r.t. to the calling routine. For instance, a matrix mat(M,N) can be transfered as mat(N,M) to an auxiliary function or subroutine. The ordering of mat changes from mat(m + M(n-1)) to mat(n + N(m-1)) in the called module.

In order to keep the pass-by-reference behavior when C-wrapping, it is necessary to use f2py_ directives in the FORTRAN code. Here is an example of the chosen implementation strategy for the subroutine :f:subr:`~pppack/spli2d`.

.. code-block:: fortran
  
  subroutine spli2d ( tau, gtau, t, n, k, m, work, q, bcoef, iflag )
  
  ! ...
  
  !f2py integer, depend(gtau) :: n = shape(gtau,0)
  !f2py integer, check(n*m == shape(gtau,0)*shape(gtau,1)), depend(gtau,n) :: m = shape(gtau,1)
  integer ( kind = 4 ), intent(in) :: m, n
  
  !f2py real ( kind = 8 ), intent(in), dimension(*,*) :: gtau
  real ( kind = 8 ), intent(in) :: gtau(n,m)
  
  ! ...
  
  end subroutine

With the directives above there is basic check of array dimensions, and an exception is raised in case of mismatch of `gtau` element types or `gtau` is not contiguous. Pearu Peterson, the main f2py_ developer, is greatly acknowledged for suggesting this implementation scheme.

Please, remind that copies are avoided only if the dummy argument is FORTRAN contiguous and the element types match exactly. In all other cases, data is copied. So use :py:class:`numpy.ndarray` with parameter `order='F'` and `dtype=np.float64` for floats and `dtype=np.int32` for integers, according to the FORTRAN libraries of pppack (`real` in double precision with `kind = 8`, and `integer` of `kind = 4`). However, internal type cast is implemented wherever possible to help the users. When array copies are produced by f2py, a message is prompted as a warning.

.. deprecated:: 1.0.0
  
  A less efficient alternative is available in :py:meth:`.fd.cmpcoef` and :py:meth:`.fd.__call__` with explicit calls to array vectorization and matricization (by :py:meth:`~pppack.vectorize` and :py:meth:`~pppack.matricize` respectively) to implement correctly the calculation of the tensor coefficients and the evaluation of the multi-dimensional functions. Code changes are kept commented as reminder of previous implementation (in pure Python) and internal data organization.

.. _f2py: https://docs.scipy.org/doc/numpy-dev/f2py/


Topical tests
-------------

Here are a few examples to understand the usage of the package.

Esteban, I'm all ears with the tests.

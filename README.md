# face-lattice

Simple python script to enumerate the set of k-faces, based on the algorithm described by Loechner & Wilde, 1997:
https://link.springer.com/article/10.1023/A:1025117523902

### Uses
* islpy - https://pypi.org/project/islpy/
* Polylib - http://www.irisa.fr/polylib/

### Setup steps
* install the PolyLib libraries from the link above
* compile the polylib/test.c file to produce the constraints2rays binary (just use the polylib/Makefile)
* place the polylib/constraints2rays and polylib/constraints2rays_wrapper.sh in your environment PATH
* install the islpy python package (use the requirements.txt to optionally create a virtual environment)

### Example

If everything is set up correctly, then you should be able to do the following from a python session:
```
>>> from face_lattice import *
>>>
>>> face_lattice('[N] -> { [i,j] : 0<=i,j and i+j<=N and N>=3 }')
[N] -> { [i,j] : 0<=i,j and i+j<=N and N>=3 }

C_hat (constraints):
[[ 1  0  0  0]		[N] -> { [i, j] : i >= 0 }
 [ 0  1  0  0]		[N] -> { [i, j] : j >= 0 }
 [-1 -1  1  0]		[N] -> { [i, j] : N - i - j >= 0 }
 [ 0  0  1 -3]		[N] -> { [i, j] : -3 + N >= 0 }
 [ 0  0  0  1]]		

R_hat (rays/vertices):
[[0 0 1 3 0 0]
 [0 1 0 0 0 3]
 [1 1 1 3 3 3]
 [0 0 0 1 1 1]]

Dimension:
3

k-faces:
k=3
(array([1., 1., 1., 1., 1., 1.]), set())

k=2
(array([1., 1., 0., 0., 1., 1.]), {0})
(array([1., 0., 1., 1., 1., 0.]), {1})
(array([0., 1., 1., 1., 0., 1.]), {2})
(array([0., 0., 0., 1., 1., 1.]), {3})

k=1
(array([1., 0., 0., 0., 1., 0.]), {0, 1})
(array([0., 1., 0., 0., 0., 1.]), {0, 2})
(array([0., 0., 0., 0., 1., 1.]), {0, 3})
(array([0., 0., 1., 1., 0., 0.]), {1, 2})
(array([0., 0., 0., 1., 1., 0.]), {1, 3})
(array([0., 0., 0., 1., 0., 1.]), {2, 3})

k=0
(array([0., 0., 0., 0., 1., 0.]), {0, 1, 3})
(array([0., 0., 0., 0., 0., 1.]), {0, 2, 3})
(array([0., 0., 0., 1., 0., 0.]), {1, 2, 3})
```

To summarize, the `face_lattice(s)` function does the following:
* constructs an ISL set based on the input string `s`
* extracts the constraints matrix `C_hat` from ISL
* passes the constraints matrix to PolyLib to produce the rays/verticies representation `R_hat`. Each column in `R_hat` represents a ray or vertex with the last row indicating which one it is (0 for rays, 1 for verticies). In this example, there are 3 rays and 3 verticies.
* constructs the incidence matrix from `C_hat` and `R_hat`
* make several calls to ScanFaces (based on Loechner & Wilde's recursive ScanFaces algorithm) to compute the various k-faces using the incidence matrix 

The example above has a 3-dimensional combined index and parameter space. There are four 2-faces (i.e., planes):
```
k=2
(array([1., 1., 0., 0., 1., 1.]), {0})
(array([1., 0., 1., 1., 1., 0.]), {1})
(array([0., 1., 1., 1., 0., 1.]), {2})
(array([0., 0., 0., 1., 1., 1.]), {3})
```

The first reported 2-face is described by the 0th, 1st, 4th, and 5th columns in `R_hat`, i.e., at the points/vectors (i,j,N):
* (0,0,1)
* (0,1,1)
* (0,0,3)
* (0,3,3)

and comes about when saturating the 0th constraint (i.e., i=0).

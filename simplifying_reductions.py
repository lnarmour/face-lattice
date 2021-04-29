from face_lattice import *


def recursive_simplify(k, op=None, fp=None, fd=None, node=None, lattice=None):
    labels = []

    for face in lattice.graph.neighbors(node):
        print(face)
    



def simplify(op=None, fp=None, s=None, fd=None):
    C_hat, lattice, bset, dim_P = face_lattice(s)

    k = dim_P
    root = lattice.get_root()
    recursive_simplify(k, op=op, fp=fp, fd=fd, node=root, lattice=lattice)


def main():
    simplify(op='+',
             fp='{[i,j,k]->[i]}',
              s='{[i,j,k] : 1<=i<=100 and 1<=j<=i-1 and 1<=k<=i-j }',
             fd='{[i,j,k]->[k]}')


if __name__ == '__main__':
    main()
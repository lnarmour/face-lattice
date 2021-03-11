from islpy import BasicSet
from islpy import dim_type
from subprocess import Popen
from subprocess import PIPE
import numpy as np


def intersect(V, I_i):
    Vp = V.copy()
    where_0 = I_i != V
    Vp[where_0] = 0
    return Vp


def card(V):
    return np.count_nonzero(V == 1)


def ScanFaces(k=None, i=None, needed=None, F=None, V=None, I=None):
    f = len(I) - 1
    if needed == 0:
        for j in range(f):
            if j not in F and np.array_equal(intersect(V, I[j]), V):
                return
        print((V, F))
        return
    if i + needed > f:
        return
    Vp = intersect(V, I[i])
    if card(Vp) >= k + needed:
        Fp = F.copy()
        Fp.add(i)
        ScanFaces(k, i+1, needed-1, Fp, Vp, I)
    if i + needed >= f:
        return
    ScanFaces(k, i+1, needed, F, V, I)


def make_incidence_matrix(C, R):
    I = np.matmul(C, R)
    where_0 = np.where(I == 0)
    where_not_0 = np.where(I != 0)
    I[where_0] = 1
    I[where_not_0] = 0
    return I


def make_constraints_matrix(bset):
    C = []
    for c in bset.get_constraints():
        vars = c.get_var_dict()
        index_vals = [int(c.get_coefficient_val(v[0], v[1]).to_str()) for k, v in vars.items() if
                      v[0] != dim_type.param]
        param_vals = [int(c.get_coefficient_val(v[0], v[1]).to_str()) for k, v in vars.items() if
                      v[0] == dim_type.param]
        const_val = [int(c.get_constant_val().to_str())]
        polylib_prefix = [0] if c.is_equality() else [1]
        C.append(polylib_prefix + index_vals + param_vals + const_val)

    C.append([0] * (1 + len(index_vals) + len(param_vals)) + [1])
    C = np.array(C)

    where_inequality = C[:,0] == 1
    index_param_space = C[where_inequality][:,1:-1]
    dim_P = np.linalg.matrix_rank(index_param_space)

    return C, dim_P


def make_rays_matrix(C):
    polylib_C = C[:-1, :]

    M, N = polylib_C.shape
    ret = '{} {}\n'.format(M, N)
    for i in range(M):
        for j in range(N):
            ret += '{} '.format(polylib_C[i, j])
        ret += '\n'
    tmp_file = '/tmp/.polylib-main'
    with open(tmp_file, 'w') as f:
        f.write(ret)

    # this is a call to two PolyLib functions: Constraints2Polyhedron, then Polyhedron2Rays
    pipe = Popen(['constraints2rays_wrapper.sh'], stdout=PIPE)
    result_bytes = pipe.stdout.read()

    R = []
    for line in result_bytes.decode('utf-8').split('\n')[1:]:
        if not line:
            continue
        ray = [int(r) for r in line.split(' ')]
        R.append(ray)

    R = np.array(R)
    return R


def pretty_print_constraints(bset, C_hat):
    ret = '\n'.join([str(c) for c in bset.get_constraints()]) + '\n'
    print('C_hat (constraints):')
    for row, constraint in zip(str(C_hat).split('\n'), ret.split('\n')):
        print('{}\t\t{}'.format(row, constraint))


def initialize_facets(C):
    F = set()
    # put all equalities from C_hat in F
    # equalities in C_hat have a value of 0 in the first column, excluding the last row
    where_equality = np.where(C[:-1, 0] == 0)
    for k in where_equality[0]:
        F.add(k)
    return F


def face_lattice(s):
    if not s:
        return
    print(s)
    print()
    bset = BasicSet(s)
    C, dim_P = make_constraints_matrix(bset)
    R = make_rays_matrix(C)

    C_hat = C[:, 1:]
    R_hat = R[:, 1:].transpose()

    pretty_print_constraints(bset, C_hat)
    print('\nR_hat (rays/vertices):\n{}\n'.format(R_hat))

    I = make_incidence_matrix(C_hat, R_hat)

    F = initialize_facets(C)
    V = np.ones(I.shape[1]) # all at first

    print('Dimension:\n{}\n'.format(dim_P))
    print('k-faces:')

    # find all k-faces of P
    for k in range(dim_P, -1, -1):
        print('k={}'.format(k))
        ScanFaces(k=k, i=0, needed=dim_P-k, F=F, V=V, I=I)
        print()

if __name__ == '__main__':
    face_lattice('{ [i,j,k] : 0<=i,j,k and k+j>=i }')

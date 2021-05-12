from face_lattice import *
from islpy import *
from enum import Enum
from itertools import combinations
from functools import reduce


def ex1():
    # ex1
    simplify(op='max',
             fp='{[i,j,k]->[i]}',
              s='{[i,j,k] : 1<=i<=100 and 1<=j<=i-1 and 1<=k<=i-j }',
             fd='{[i,j,k]->[k]}')


def ex3():
    # ex3 - needs reduction decomposition with not-so-obvious factorization of fp
    simplify(op='max',
             #fp='{[i,j,k]->[i,j+k]}',
             fp='{[i,j,k]->[i]}',
             s='{[i,j,k] : j>=i and 2i>=j and k>=i and 3i>=j+k and j>=1 and k>=1 and 0>=j-100 and 0>=k-100 and 0>=i-100 and i>=1}',
             fd='{[i,j,k]->[j,k]}')


def inverse(op):
    store = {
        '+': '-',
        'max': None
    }
    if op not in store:
        raise Exception('Operator "{}" not supported yet.'.format(op))
    return store[op]


# creates the matrix representation of f
# if f = '{[i,j,k]->[i]}', then this returns [[1,0,0]]
# if f = '{[i,j,k]->[i+k,j]}', then this returns [[1,0,1],[0,1,0]]
def map_to_matrix(f):
    mat = []
    for c in f.get_constraints():
        vars = c.get_var_dict()
        index_vals = [-1 * int(c.get_coefficient_val(v[0], v[1]).to_str()) for k, v in vars.items() if
                      v[0] != dim_type.param]
        mat.append(index_vals)
    return np.array(mat)


def mat_to_set(mat, condition='='):
    indices = ','.join(['i{}'.format(i) for i in range(len(mat[0,:]))])
    constraints = []
    for c in mat:
        constraints.append('+'.join(['{}{}'.format(a[0],a[1]) for a in zip(c, indices.split(','))]) + '{}0'.format(condition))
    s = '{{[{}] : {}}}'.format(indices, ' and '.join(constraints))
    return BasicSet(s)


def ker(f):
    if type(f) == str:
        f = BasicMap(f)
    mat = []
    for c in f.get_constraints():
        vars = c.get_var_dict()
        index_vals = [-1 * int(c.get_coefficient_val(v[0], v[1]).to_str()) for k, v in vars.items() if
                      v[0] != dim_type.param]
        mat.append(index_vals)
    mat = np.array(mat)
    indices = ','.join(['i{}'.format(i) for i in range(len(mat[0, :]))])
    constraints = []
    for c in mat:
        constraints.append(
            '+'.join(['{}{}'.format(a[0], a[1]) for a in zip(c, indices.split(','))]) + '=0')
    s = '{{[{}] : {}}}'.format(indices, ' and '.join(constraints))
    return BasicSet(s)


def is_strict(face, fp, C):
    c_mat = C[np.array(list(face)),:-1]
    fp_mat = map_to_matrix(fp)
    ker_c = mat_to_set(c_mat)
    ker_fp = mat_to_set(fp_mat)
    return not ker_c.intersect(ker_fp).is_empty() and ker_fp.is_subset(ker_c)


LABEL = {
    0: 'ADD',
    1: 'INV',
    2: 'SUB'
}


def d_ary(d, n):
    if n == 0:
        return '0'
    nums = []
    while n:
        n, r = divmod(n, d)
        nums.append(str(r))
    return ''.join(reversed(nums))


def rho_from_labels(faces, C, labels, fd):
    # check if this combo is possible given the feasible_rho
    # c1     c2     c3
    # 'ADD', 'ADD', 'INV'
    # c1*rho>=0 and c2*rho>=0 and c3*rho=0  and in ker(fp) and it saturates current node
    # {[i,j,k] : j>=0 and i-j-k>=0 and k=0  and i=0   }  <- is this empty?  if yes, then prune out this combo
    # if no, then check that it's "in ker(fp)
    map = {'ADD': '>', 'INV': '=', 'SUB': '<'}
    bsets = list()

    # must be in ker(fd)
    bsets.append(mat_to_set(map_to_matrix(fd)))

    for face, label in zip(faces, labels):
        # select constraints representing this face, and create the set according to its label
        bset = mat_to_set(C[np.array(list(face)),:-1], map[label])
        bsets.append(bset)

    # check that it's not trivially just the origin
    origin = BasicSet('{{[{}]}}'.format(('0,'*(C.shape[1]-1))[:-1]))
    result = reduce(lambda s0,s1: s0.intersect(s1), bsets)

    feasible_rho = result - origin

    if feasible_rho.is_empty():
        return None

    rho = list()
    # either lexmin or lexmax
    funcs = [feasible_rho.lexmin, feasible_rho.lexmax]
    for func in funcs:
        try:
            func().foreach_point(rho.append)
        except Exception:
            continue
        finally:
            if rho:
                break

    return rho[0]


def recursive_simplify(k=None, fp_str=None, fd_str=None, node=None, lattice=None, C=None, legal_labels=None):
    fp = BasicMap(fp_str)
    fd = BasicMap(fd_str)

    print('STEP 3.0 - recursive call with:\n')

    print('node = {}'.format('{}'))
    print('fp = {}'.format(fp_str))
    print('fd = {}'.format(fd_str))
    print()

    print('STEP 3.1 - for each child face in node, label as "weak" or "strict" given fp\n')
    faces = list(lattice.graph.neighbors(node))
    for face in faces:
        face = set(face)
        label = 'strict' if is_strict(face, fp, C) else 'weak'
        print(face, label)
    print()

    candidate_faces = [face for face in faces if not is_strict(face, fp, C)]
    #candidate_faces = faces

    print('STEP 3.2 - determine all possible combos\n')

    label_combos = list()
    header = ['{}'.format(set(f)) for f in candidate_faces]
    print(header)
    print('-' * len(str(header)))
    for i in range(len(legal_labels) ** len(candidate_faces)):
        labels = [LABEL[int(d_ary(len(legal_labels), i).zfill(len(candidate_faces))[j])] for j in range(len(candidate_faces))]
        rho = rho_from_labels(candidate_faces, C, labels, fd)
        if rho:
            label_combos.append(labels + [rho])
            print('{}  possible -> rho = {}'.format(labels, rho))
        else:
            print('{}'.format(labels))
    print()

    print('STEP 3.3 - construct complete possible labeling w/ any strict faces\n')



def simplify(op=None, fp=None, s=None, fd=None):
    print('-'*80)

    print('STEP 1 - construct the face lattice from context domain of reduction body\n')
    C, lattice, bset, dim_P = face_lattice(s)

    print('STEP 2 - set up recursion\n')

    # todo add counter for num remaining dims of reuse available
    # not necessarily a need to decompose the dependence function

    legal_labels = ['ADD', 'INV']
    if inverse(op):
        legal_labels.append('SUB')

    print('op = {}'.format(op))
    print('labels = {}'.format(legal_labels))
    print()

    root = lattice.get_root()

    recursive_simplify(k=dim_P, fp_str=fp, fd_str=fd, node=root, lattice=lattice, C=C, legal_labels=legal_labels)


def ex_manual():
    # ex2 - needs reduction decomposition of fp first
    simplify(op='max',
             fp='{[i,j,k]->[i,k]}',
             s='{[i,j,k] : j>=1 and i>=j and k>=1 and 0>=i+k-100 and 0>=k-100 and 0>=j-100 and 0>=i-99 and i>=1}',
             fd='{[i,j,k]->[j,k]}')


if __name__ == '__main__':
    ex_manual()













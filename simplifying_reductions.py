from face_lattice import *
from islpy import *
from enum import Enum
from itertools import combinations
from functools import reduce



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

    return feasible_rho


def recursive_simplify(k=None, op=None, fp=None, fd=None, node=None, lattice=None, C=None):

    print('STEP 3 - for each child face of the current node, label as "weak" or "strict"\n')
    faces = list(lattice.graph.neighbors(node))
    for face in faces:
        face = set(face)
        label = 'strict' if is_strict(face, fp, C) else 'weak'
        print(face, label)
    print()

    # ignore strict faces
    candidate_faces = [face for face in faces if not is_strict(face, fp, C)]

    print('STEP 4 - determine the feasible space of legal reuse vectors\n')
    #   - if the operator admits an inverse, then this is given by the kernel of fd
    #   - else, construct the feasible reuse space s.t. no subtraction faces can be induced
    # this is a set of inequalities, all non-boundary faces must be INV or ADD (i.e., rho dot c >= 0)
    if inverse(op):
        legal_labels = ['ADD', 'SUB', 'INV']
    else:
        legal_labels = ['ADD', 'INV']

    print('STEP 5 - construct set of all possible label combos & corresponding rho\n')

    num_faces = len(candidate_faces)
    label_combos = list()
    num_labels = len(legal_labels)
    for i in range(num_labels ** len(candidate_faces)):
        labels = [LABEL[int(d_ary(num_labels, i).zfill(num_faces)[j])] for j in range(num_faces)]
        feasible_rho = rho_from_labels(candidate_faces, C, labels, fd)
        if feasible_rho.is_empty():
            print(labels)
            continue
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

        rho = rho[0]
        print('{}  possible ->   rho = {}'.format(labels, rho))
        label_combos.append(labels + [rho])
    print()

    print('STEP 6 - remove redundant combinations\n')
    for combo in label_combos:
        print('{}  possible ->   rho = {}'.format(combo[:-1], combo[-1]))
    print()
    # 1 'ADD', 'ADD', 'INV'  <-- recursive calls to the first ADD, then the second ADD are successful
    # 2 'INV', 'INV', 'ADD'

    # STEP 7 - for a given combo, recurse into each ADD face
    #   each recursive call must return a successful simplification
    #   need some "global" data structure (i.e. a table at each node to keep track of which paths have already been explored)

    # double check this
    # STEP 7.5, before each recursive call, do decomposition to construct new alphaz expression
    #   construct new projection function
    #   unroll constant dims, domain must have equalities







    # select minimal possible ones





def simplify(op=None, fp=None, s=None, fd=None):
    print('-'*80)

    print('STEP 1 - construct the face lattice\n')
    C, lattice, bset, dim_P = face_lattice(s)
    print()

    print('STEP 2 - get the root node in lattice and start recursion\n')
    root = lattice.get_root()
    print('node = {}'.format('{}'))
    print('fp = {}'.format(fp))
    print('fd = {}'.format(fd))
    print()

    # todo add counter for num remaining dims of reuse available
    # not necessarily a need to decompose the dependence function

    recursive_simplify(k=dim_P, op=op, fp=BasicMap(fp), fd=BasicMap(fd), node=root, lattice=lattice, C=C)


def main():

    # ex1
    simplify(op='max',
             fp='{[i,j,k]->[i]}',
              s='{[i,j,k] : 1<=i<=100 and 1<=j<=i-1 and 1<=k<=i-j }',
             fd='{[i,j,k]->[k]}')

    # ex2
    simplify(op='+',
             fp='{[i,j,k]->[i]}',
             s='{[i,j,k] : k>=1 and 0>=i+k-100 and j>=1 and i>=j and 0>=k-100 and 0>=j-100 and i>=1 and 0>=i-99}',
             fd='{[i,j,k]->[j,k]}')


if __name__ == '__main__':
    main()












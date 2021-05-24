from face_lattice import *
from islpy import *
from copy import deepcopy
from functools import reduce
import builtins
import json


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


def ker_from_map(f):
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


def ker_from_facet_normal(facet, parent, C):
    np_C = np.array(C)
    mat = np_C[np.array(list(facet-parent)),:-1]
    indices = ','.join(['i{}'.format(i) for i in range(len(mat[0, :]))])
    constraints = []
    for c in mat:
        constraints.append('+'.join(['{}{}'.format(a[0], a[1]) for a in zip(c, indices.split(','))]) + '=0')
    s = '{{[{}] : {}}}'.format(indices, ' and '.join(constraints))
    return BasicSet(s)


def is_strict(facet, parent, fp, C):
    ker_c = ker_from_facet_normal(facet, parent, C)
    ker_fp = ker_from_map(fp)
    return not ker_c.intersect(ker_fp).is_empty() and ker_fp.is_subset(ker_c)


# convert base-10 number n to base-d
def d_ary(d, n):
    if n == 0:
        return '0'
    nums = []
    while n:
        n, r = divmod(n, d)
        nums.append(str(r))
    return ''.join(reversed(nums))

LABEL = {
    0: 'ADD',
    1: 'INV',
    2: 'SUB'
}


def enumerate_labels(facets, labels):
    num_labels = len(labels)
    num_facets = len(facets)
    for i in range(num_labels ** num_facets):
        labels = [LABEL[int(d_ary(num_labels, i).zfill(num_facets)[j])] for j in range(num_facets)]
        yield labels


def rho_from_labels(faces, parent, Lp, C, labels, fd):
    # check if this combo is possible given the feasible_rho
    # c1     c2     c3
    # 'ADD', 'ADD', 'INV'
    # c1*rho>=0 and c2*rho>=0 and c3*rho=0  and in ker(fp) and it saturates current node
    # {[i,j,k] : j>=0 and i-j-k>=0 and k=0  and i=0   }  <- is this empty?  if yes, then prune out this combo
    # if no, then check that it's "in ker(fp)
    map = {'ADD': '>', 'INV': '=', 'SUB': '<'}
    bsets = list()
    bsets.append(Lp)
    # must be in ker(fd)
    bsets.append(mat_to_set(map_to_matrix(fd)))

    for face, label in zip(faces, labels):
        # select constraints representing this face, and create the set according to its label
        bset = mat_to_set(C[np.array(list(face-parent)),:-1], map[label])
        bsets.append(bset)

    # check that it's not trivially just the origin
    origin = BasicSet('{{[{}]}}'.format(('0,'*(C.shape[1]-1))[:-1]))
    result = reduce(lambda s0,s1: s0.intersect(s1), bsets)

    feasible_rho = result - origin

    if not feasible_rho:
        return None

    rho = list()
    # either lexmin or lexmax, not sure yet how to determine upfront which one to use, so for now just try both
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


def prune_combos(label_combos):
    combos = [[1 if l == 'ADD' else 0 for l in c[:-1]] for c in label_combos]
    combos.sort(key=lambda c: np.sum(c), reverse=1)
    # sort by combos with most ADD faces first
    # discard a combo if its ADD faces can
    ret = []
    for i,combo in enumerate(combos):
        _combo = np.array(combo)
        for other in combos[i+1:]:
            other = np.array(other)
            _combo = _combo - other
        if 1 in _combo:
            ret.append(combo)
    # TODO - return only label_combos rows that match ret
    unique_combos = []
    for lc in label_combos:
        labels = [1 if l == 'ADD' else 0 for l in lc[:-1]]
        if labels in ret:
            unique_combos.append(lc)
    return unique_combos


# take an isl_point obj and parse it as a vector
# given {Point}{ [1, 1, 0] }
# return [1, 1, 0]
def point_to_vec(isl_point):
    N = len(isl_point.get_var_dict())
    vd = isl_point.get_var_dict()
    vec = [int(isl_point.get_coordinate_val(vd[v][0], vd[v][1]).to_str()) for v in vd]
    return vec


def boundary_label(facet, parent, C, rho):
    c_mat = C[np.array(list(facet-parent)), :-1]
    rho_vec = point_to_vec(rho)
    orientations = np.matmul(c_mat, rho_vec)
    zero_vec = np.zeros(len(c_mat))
    return 'inward' if orientations >= zero_vec else 'outward'


# theorem 2 - Lp is the intersection of all of the effectively saturated constraints
def get_Lp(node, C):
    # start with universe
    Lp = BasicSet('{{[{}]}}'.format(','.join(['i{}'.format(i) for i in range(C.shape[1] - 1)])))
    np_C = np.array(C)
    for c in list(node):
        mat = np_C[[c], :-1]
        indices = ','.join(['i{}'.format(i) for i in range(len(mat[0, :]))])
        constraints = []
        for c in mat:
            constraints.append('+'.join(['{}{}'.format(a[0], a[1]) for a in zip(c, indices.split(','))]) + '=0')
        s = '{{[{}] : {}}}'.format(indices, ' and '.join(constraints))
        Lp = Lp.intersect(BasicSet(s))
    return Lp


def simplify(k=None, fp_str=None, fd_str=None, node=None, lattice=None, C=None, legal_labels=None, rho=None):
    def pprint(*args, **kwargs):
        print('@{} '.format(set(node) if node else '{}'), end='')
        print(*args, **kwargs)

    pprint('STEP A.1 - if k==0 return else continue. k={}'.format(k))
    if k == 0:
        pprint()
        pprint('Success - reached bottom, no more available dimensions of reuse.')
        pprint()
        return True
    pprint()

    fp = BasicMap(fp_str)
    fd = BasicMap(fd_str)
    pprint('node = {}'.format(set(node)))
    pprint('fp = {}'.format(fp_str))
    pprint('fd = {}'.format(fd_str))
    pprint()

    pprint('STEP A.2 - identify "strict" boundaries given fp')
    pprint()
    facets = list(lattice.graph.neighbors(node))
    for facet in facets:
        label = 'boundary' if is_strict(facet, node, fp, C) else ''
        pprint(set(facet), label)
    pprint()

    pprint('STEP A.3 - construct list of candidate facets (i.e., non-boundary facets)')
    pprint()
    candidate_facets = [facet for facet in facets if not is_strict(facet, node, fp, C)]
    pprint('candidate_facets = {}'.format([set(cf) for cf in candidate_facets]))
    pprint()

    pprint('STEP A.4 - determine all possible combos')
    pprint()
    label_combos = list()
    header = ['{}'.format(set(f)) for f in candidate_facets]
    pprint(header)
    pprint('-' * len(str(header)))
    # theorem 2
    Lp = get_Lp(node, C)
    for labels in enumerate_labels(candidate_facets, legal_labels):
        rho = rho_from_labels(candidate_facets, node, Lp, C, labels, fd)
        if rho:
            label_combos.append(labels + [rho])
            pprint('{}  possible -> rho = {}'.format(labels, rho))
        else:
            pprint('{}  impossible'.format(labels))
    pprint()
    if not label_combos:
        pprint('FAILURE - no possible combos')
        return False

    pprint('STEP A.5 - prune out redundant possible combos')
    pprint()
    pprint(header)
    pprint('-' * len(str(header)))
    unique_label_combos = prune_combos(label_combos)
    for combo in unique_label_combos:
        labels,rho = combo[:-1],combo[-1]
        pprint('{}  -> rho = {}'.format(labels, rho))
    pprint()

    pprint('STEP A.6 - incorporate boundary facets')
    pprint()
    boundary_facets = [f for f in facets if f not in candidate_facets]
    header = '{} {}'.format(header, [str(set(bf)) for bf in boundary_facets])
    pprint(header)
    pprint('-' * len(str(header)))
    full_label_combos = []
    for combo in unique_label_combos:
        labels,rho = combo[:-1],combo[-1]
        boundary_labels = []
        full_label_combo = deepcopy(labels)
        for boundary_facet in boundary_facets:
            bl = boundary_label(boundary_facet, node, C, rho)
            boundary_labels.append(bl)
            full_label_combo.append(bl)
            pprint('{} {}  -> rho = {}'.format(labels, boundary_labels, rho))
        full_label_combo.append(rho)
        full_label_combos.append(full_label_combo)
    pprint()

    pprint('STEP A.7 - recurse into "ADD" and "inward" boundary facets')
    pprint()
    successful_combos = []
    for combo in full_label_combos:
        labels,rho = combo[:-1],combo[-1]
        abort = None
        for label,facet in zip(labels,candidate_facets + boundary_facets):
            if label != 'ADD' and label != 'inward':
                continue
            pprint('recursing into {} facet'.format(set(facet)))
            pprint()
            ret = simplify(k=k-1, fp_str=fp_str, fd_str=fd_str, node=facet, lattice=lattice, C=C, legal_labels=legal_labels)
            abort = not ret
            if abort:
                break
            # TODO - figure out how to propagate the success back up the recusion
            # TODO - likely, attach to the face lattice accordingly
        if not abort:
            successful_combos.append(combo)
        pprint()
    if len(successful_combos) == 0:
        pprint('FAILURE - no successful combos')
    pprint()

    pprint('STEP B - todo...')
    pprint()
    # TODO - implement step 6b from Algorithm 2 in 2006 paper

    pprint('STEP C - reduction decomposition - todo...')
    pprint()
    # TODO - figure out how to do the decomposition
    # TODO - i.e., how to navigate the infinite space of decompositions

    return successful_combos


def ex_manual():
    # ex2 - needs reduction decomposition of fp first
    op = 'max'
    fp = '{[i,j,k]->[i]}'
    s = '{[i,j,k] : j>=1 and i>=j and k>=1 and 0>=i+k-100 and 0>=k-100 and 0>=j-100 and 0>=i-99 and i>=1}'
    fd = '{[i,j,k]->[j,k]}'

    # ex1
    #op = 'max'
    #fp = '{[i,j,k]->[i]}'
    #s = '{[i,j,k] : j>=1 and i>=j+1 and k>=1 and i>=j+k and 0>=k-100 and i>=2 and 0>=i-100}'
    #fd = '{[i,j,k]->[k]}'


    legal_labels = ['ADD', 'INV']
    if inverse(op):
        legal_labels.append('SUB')

    C, lattice, bset, dim_P = face_lattice(s)
    root = lattice.get_root()
    node = root

    simplify(k=2, fp_str=fp, fd_str=fd, node=node, lattice=lattice, C=C, legal_labels=legal_labels)


if __name__ == '__main__':
    ex_manual()













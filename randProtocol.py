from scipy.optimize import linprog


def to_add(x):
    """ takes an integer x as input,
    finds the position of least significant 0,
    returns an integer with only set bit at that position """
    c = 1
    while x & 1 == 1:
        c = c << 1
        x = x >> 1
    return c


def get_extensions(x):
    """ takes an integer x,
    returns a list of size x+1, at each index of which is a list of all the extensions of the corresponding
    rowset/columnset """
    extensions = [[]]
    for i in range(x - 1, 0, -1):
        temp = []
        j = to_add(i)
        y = i + j
        temp.append(y)
        for p in extensions[x - y]:
            temp.append(p)
            temp.append(p - j)
        extensions.append(temp)
    extensions.append([])
    extensions.reverse()
    return extensions


def extend(rectangle, row_extension, col_extension):
    """
    inputs:
    rectangle: a tuple with rowset at first index and columnset at second
    row_extension/col_extension: a list of list of all the extensions of each rowset/columnset

    returns  two lists
    1. a list of all the rectangles (as a tuple) that is an extension of the input rectangle
    2. a boolean list that tells if its a row extension or not
    """
    extensions = []
    is_row_ext = []
    row, col = rectangle[0], rectangle[1]
    for i in row_extension[row]:
        extensions.append((i, col))
        is_row_ext.append(True)
    for j in col_extension[col]:
        extensions.append((row, j))
        is_row_ext.append(False)
    return extensions, is_row_ext


def all_rectangles(matrix):
    """
    takes the matrix as input
    returns a list of all the rectangles (as a tuple of rowset and columnset) contained in the matrix
    """
    m, n = 1 << len(matrix), 1 << len(matrix[0])
    rectangles = []
    for i in range(1, m):
        for j in range(1, n):
            rectangles.append((i, j))
    return rectangles


def get_partitions(x):
    """
    takes an integer x,
    returns a list of size x+1,
    at each index of which is a list of all the partitions of the corresponding rowset/columnset
    """
    partitions = [[]]
    for i in range(1, x+1):
        temp = []
        least = i & -i
        if least == i:
            temp.append(0)
        else:
            y = i - least
            for p in partitions[y]:
                temp.append(p)
                temp.append(p + least)
        partitions.append(temp)
    for i in range(1, x+1):
            partitions[i].remove(0)
    return partitions


def all_outer_rectangles(rectangle, row_extension, col_extension):
    """
    inputs:
    rectangle: a tuple with rowset at first index and columnset at second
    row_extension/col_extension: a list of list of all the extensions of each rowset/columnset

    returns a list of all rectangles that contains the input rectangle including itself
    """
    outer_rectangles = []
    row, col = rectangle[0], rectangle[1]
    for i in row_extension[row] + [row]:
        for j in col_extension[col] + [col]:
            outer_rectangles.append((i, j))
    return outer_rectangles


def maximize(matrix, steps, one_sided=False):
    """
    inputs:
    matrix of output values,
    number of steps,
    a boolean value that tells if one-sided error is to be calculated

    returns the success probability with which the output can be communicated in given steps at min
    """

    eq_constraints = []
    eq_rhs = []
    a, b = len(matrix), len(matrix[0])
    m, n = (1 << a) - 1, (1 << b) - 1
    all_partitions = get_partitions(max(m, n))
    indices = {}
    index = 1
    for k in range(steps + 1):
        for r in all_rectangles(matrix):
            partitions = []
            if k < steps:
                for partition in all_partitions[r[0]]:
                    partitions.append((partition, r[1]))
                for partition in all_partitions[r[1]]:
                    partitions.append((r[0], partition))
            partitions.append(0)
            partitions.append(1)
            for q in partitions:
                indices[(k, r, q)] = index
                index += 1
    row_ext = get_extensions(m)
    col_ext = get_extensions(n)
    for k in range(steps + 1):
        for r in all_rectangles(matrix):
            extensions, is_row_ext = extend(r, row_ext, col_ext)
            partitions = []
            if k < steps:
                for partition in all_partitions[r[0]]:
                    partitions.append((partition, r[1]))
                for partition in all_partitions[r[1]]:
                    partitions.append((r[0], partition))
            partitions.append(0)
            partitions.append(1)
            constraint = [0]*index
            for q in partitions:
                constraint[indices[(k, r, q)]] = 1
            rhs = 0
            if k == 0:
                if not extensions:
                    rhs = 1  # RHS is 1 when k = 0 and R is a full rectangle
            else:
                if extensions:
                    for i, ext in enumerate(extensions):
                        if (k-1, ext, r) in indices:
                            constraint[indices[(k - 1, ext, r)]] = -1
                        else:  # if r is not a part in the partition list, look for ext\r
                            if is_row_ext[i]:
                                r1 = (ext[0] - r[0], r[1])
                            else:
                                r1 = (r[0], ext[1] - r[1])
                            constraint[indices[(k - 1, ext, r1)]] = -1
            eq_rhs.append(rhs)
            eq_constraints.append(constraint)
    ineq_constraints = []
    ineq_rhs = []

    for i in range(a):
        for j in range(b):
            q = matrix[i][j]
            if q is not None:
                inp = (1 << i, 1 << j)
                constraint = [0] * index
                constraint[0] = 1
                rhs = 0
                if one_sided:
                    if q == 0:
                        constraint[0] = 0
                        rhs = -1  #success probability is 1 for 0-inputs
                for r in all_outer_rectangles(inp, row_ext, col_ext):
                    for k in range(steps + 1):
                        constraint[indices[(k, r, q)]] = -1
                ineq_constraints.append(constraint)
                ineq_rhs.append(rhs)
    bound = [(0, 1)]*index
    bound[0] = (None, None)
    obj = [0]*index
    obj[0] = -1

    res = linprog(obj, A_ub=ineq_constraints, b_ub=ineq_rhs, A_eq=eq_constraints, b_eq=eq_rhs, bounds=bound)
    return -res.fun


def binary_search(matrix, min, max, p, one_sided=False):
    """
    returns the minimum bits to communicate with which the output can be calculated
    with probability 'p' for given matrix
    within 'min' and 'max' steps
    """
    if max == min:
        return max
    mid = int((min + max)/2)
    if maximize(matrix, mid, one_sided) < p:
        return binary_search(matrix, mid + 1, max, p)
    else:
        return binary_search(matrix, min, mid, p)


def trivial_upper_bound(matrix):
    """
    given a matrix, returns the trivial upper bound on the communication complexity based on the size of the matrix(not
    the actual matrix)
    """
    a = min(len(matrix), len(matrix[0]))
    c = 0
    result = 0
    while a != 1:
        if a & 1 == 1:
            c += 1
        result += 1
        a = a >> 1
    if c == 0:
        return result + 1
    else:
        return result + 2


def min_bits(matrix, p, one_sided=False):
    """
    returns the minimum bits to communicate with which the output is calculated with probability 'p' for given matrix
    """
    p = p - 0.01
    m = trivial_upper_bound(matrix)
    return binary_search(matrix, 0, m, p, one_sided)


matrix = [[0, 1, 0],
          [1, 0, 1],
          [1, None, None]]

m1 = [[0, 1, 1, 1],
      [1, None, 1, 0],
      [0, 0, None, 1],
      [1, 0, None, 0]]

m2 = [[0, 1], [None, 0]]

eq = [[1, 0, 0, 0],
      [0, 1, 0, 0],
      [0, 0, 1, 0],
      [0, 0, 0, 1]]

disj = [[1, 1, 1, 1],
        [1, 0, 1, 0],
        [1, 1, 0, 0],
        [1, 0, 0, 0]] 



#print(min_bits(disj, 0.5, True))

#print(max_bit(5))

print(maximize(disj, 2, True))














import math
import node


def to_binary(x, s):
    """
    given integers x and s, returns a binary string for integer x of length s
    """
    binary = ''
    for i in range(s-1, -1, -1):
        binary = binary + str((x >> i) & 1)
    return binary


def make_2d(value, k, l):
    """
    initializes a 2d array of size m * n with ever element as value
    """
    return [[value for i in range(l)] for j in range(k)]


def get_index(x):
    """
    returns log base 2 of x (in our case, x is always a power of 2)
    """
    result = 0
    while x != 1:
        result += 1
        x = x >> 1
    return result


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


class Protocol:
    def find_opt_protocol(matrix):
        """
        input: a matrix of output values of size a * b

        returns following 2d matrices of size 2^a * 2^b
        cost(integer):
        value at index i,j is the optimal cost to communicate on rectangle (i,j)

        speaker(boolean):
        value at index i,j who initiates communication on rectangle (i,j) for optimal result

        choice(integer):
        value at index i,j gives the partition of rectangle (i, j) for optimal communication

        output(integer):
        value at index i,j tell if (i,j) is a monochromatic rectangle and if it is gives the color
        -1 if non-monochromatic, else the color (0, 1 or None)
        """

        m = (1 << len(matrix))

        n = (1 << len(matrix[0]))

        partitions = get_partitions(max(m, n))

        cost = make_2d(math.inf, m, n)
        output = make_2d(-1, m, n)
        speaker = make_2d(None, m, n)
        choice = make_2d(0, m, n)

        for x in range(1, m):
            row_partitions = partitions[x]
            for y in range(1, n):
                col_partitions = partitions[y]
                if not row_partitions:
                    if not col_partitions:
                        cost[x][y] = 0
                        row_index = get_index(x)
                        col_index = get_index(y)
                        output[x][y] = matrix[row_index][col_index]
                        continue
                    else:
                        y1 = col_partitions[0]
                        y2 = y - y1
                        if cost[x][y1] == 0 and cost[x][y2] == 0:
                            if output[x][y1] == output[x][y2] or output[x][y2] is None:
                                cost[x][y] = 0
                                output[x][y] = output[x][y1]
                                continue
                            if output[x][y1] is None:
                                cost[x][y] = 0
                                output[x][y] = output[x][y2]
                                continue
                else:
                    x1 = row_partitions[0]
                    x2 = x - x1
                    if cost[x1][y] == 0 and cost[x2][y] == 0:
                        if output[x1][y] == output[x2][y] or output[x2][y] is None:
                            cost[x][y] = 0
                            output[x][y] = output[x1][y]
                            continue
                        if output[x1][y] is None:
                            cost[x][y] = 0
                            output[x][y] = output[x2][y]
                            continue
                for x1 in row_partitions:
                    x2 = x - x1
                    temp = 1 + max(cost[x1][y], cost[x2][y])
                    if cost[x][y] > temp:
                        cost[x][y] = temp
                        speaker[x][y] = True
                        choice[x][y] = x1

                for y1 in col_partitions:
                    y2 = y - y1
                    temp = 1 + max(cost[x][y1], cost[x][y2])
                    if cost[x][y] > temp:
                        cost[x][y] = temp
                        speaker[x][y] = False
                        choice[x][y] = y1
        return cost, speaker, choice, output

    def get_root(self, row_set, col_set):
        """
        given a rectangle (row_set, col_set), returns the root of the protocol with all its descendants built
        recursively
        """
        root = node.Node()
        root.row_set, root.col_set = row_set, col_set
        root.speaker = self.speaker[row_set][col_set]
        root.cost = self.cost[row_set][col_set]
        root.choice = self.choice[row_set][col_set]
        root.rectangle = to_binary(row_set, self.row_length) + ',' + to_binary(col_set, self.col_length)
        if root.cost != 0:
            if root.speaker:
                left_row = root.choice
                right_row = row_set - root.choice
                left_col = col_set
                right_col = col_set
            else:
                left_row = row_set
                right_row = row_set
                left_col = root.choice
                right_col = col_set - root.choice
            root.left = self.get_root(left_row, left_col)
            root.right = self.get_root(right_row, right_col)
        else:
            root.leaf = True
            root.output = self.output[row_set][col_set]
        return root

    def __init__(self, matrix, row_set=None, col_set=None):
        self.matrix = matrix
        self.row_length = len(matrix)
        self.col_length = len(matrix[0])
        if row_set is None:
            row_set = (1 << self.row_length) - 1
        if col_set is None:
            col_set = (1 << self.col_length) - 1
        self.cost, self.speaker, self.choice, self.output = Protocol.find_opt_protocol(matrix)
        self.root = self.get_root(row_set, col_set)

    def show_protocol(self):
        """ prints the protocol """
        self.root.print_node()


neq = [[0, 1, 1, 1],
      [1, 0, 1, 1],
      [1, 1, 0, 1],
      [1, 1, 1, 0]]

disj = [[1, 1, 0, 1],
        [1, 0, 1, 0],
        [1, None, 0, 0],
        [1, 0, None, 0]]

root = Protocol(disj)

root.show_protocol()
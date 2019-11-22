class Node:

    def __init__(self):
        self.row_set = None
        self.col_set = None
        self.speaker = None
        self.cost = None
        self.choice = None
        self.left = None
        self.right = None
        self.leaf = False
        self.output = None
        self.rectangle = None

    def print_node(self, level=0):
        """ prints a node at given level with all its descendants """
        print('\t' * level + 'Rectangle: (' + self.rectangle, end=') ')
        if self.leaf:
            print('Output: ' + str(self.output))
        else:
            if self.speaker:
                speaker = 'Alice'
            else:
                speaker = 'Bob'
            print('Speaker: ' + speaker + ' Cost: ' + str(self.cost))
            self.left.print_node(level + 1)
            self.right.print_node(level + 1)






        




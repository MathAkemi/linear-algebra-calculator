import copy

### This is the general class containing all of the methods of the calculator.
### calculations_performed, int: the number of calculations that have been performed
### Methods: matrix addition, matrix subtraction,, add constant, multiply constant,
### transpose, multiply matrices, invert, invertible, determinant
class Matrix:
    calculations_performed = 0
    
    # Initialize a new Matrix object with the matrix potentially specified.
    def __init__(self, matrix = []):
        self.calculations_performed = 0
        self.matrix = matrix
        self.r = len(matrix)
        if self.r > 0:
            self.c = len(matrix[0])

    # Get user input for the matrix; makes use of define_matrix to create it.
    def get_matrix(self):
        r = float(input("Number of rows: "))
        c = float(input("Number of columns: "))
        self.r = int(r)
        self.c = int(c)
        print("Please enter the entires of each row one at a time.")
        self.matrix = self.define_matrix()
        return self.matrix

    # Algorithmically add elements to a matrix.
    def define_matrix(self):
        matrix = []
        
        for i in range(self.r):
            row = []
            for j in range(self.c):
                row.append(float(input("Entry " + "(" + str(i) + ", " + str(j) + "): ")))
            matrix.append(row)
        
        return matrix

    # Add two matrices element-wise.
    def addition(self, matrix2):
        if isinstance(matrix2, Matrix) == False:
            raise TypeError("That's not a matrix!")
        if self.r != matrix2.r or self.c != matrix2.c:
            raise ValueError("Addition cannot be performed.")
        
        row = [None] * self.c
        returnMatrix = [row] * self.r
        for i in range(self.r):
            for j in range(self.c):
                self.matrix[i][j] = self.matrix[i][j] + matrix2.matrix[i][j]
        
        self.calculations_performed += 1
        return self

    # Subtract two matrices element-wise.
    def subtraction(self, matrix2):
        if isinstance(matrix2, Matrix) == False:
            raise TypeError("That's not a matrix!")
        if self.r != matrix2.r or self.c != matrix2.c:
            raise ValueError("Addition cannot be performed.")
        
        row = [None] * self.c
        returnMatrix = [row] * self.r
        for i in range(self.r):
            for j in range(self.c):
                self.matrix[i][j] = self.matrix[i][j] - matrix2.matrix[i][j]
        
        self.calculations_performed += 1
        return finalMatrix

    # Add a constant to a mtrix.
    def add_constant(self, c):
        if not isinstance(c, float) and not isinstance(c, int):
            raise TypeError("That's not a number!")
        for i in range(self.r):
            for j in range(self.c):
                self.matrix[i][j] = self.matrix[i][j] + c

        self.calculations_performed += 1
        return self

    # Multiply a matrix by a constant.
    def mult_constant(self, c):
        for i in range(self.r):
            for j in range(self.c):
                self.matrix[i][j] = self.matrix[i][j] * c

        self.calculations_performed += 1
        return self

    # Find the transpose of a matrix.
    def transpose(self):
        returnMatrix = []
        for i in range(self.c):
            row = []
            for j in range(self.r):
                row.append(self.matrix[j][i])
            returnMatrix.append(row)

        self.calculations_performed += 1
        transposeMatrix = Matrix(returnMatrix)
        return transposeMatrix

    # Invert a matrix if it is inertible.
    # Not yet finished!!
    #
    #def invert(self):
    #    if self.invertible() == False:
    #        raise ValueError("The matrix is not invertible.")
    #
    #    self.calculations_performed += 1
    #    return self.matrix

    # Determine if a matrix is invertible, using the determinant to check.
    def invertible(self):
        if len(self.matrix) != len(self.matrix[0]):
            raise ValueError("Matrix must be nxn.")
        if self.determinant() == 0:
            return False
        else:
            return True
        
    # Compute the determinate of a matrix.
    def determinant(self, mul=1):
        # Ensures that the matrix is square.
        if self.r != self.c:
            return ValueError("Determinant cannot be computed.")

        w = self.c

        # Takes care of the base case.
        if w == 1:
            return mul * self.matrix[0][0]
        
        # Cofactor expansion expansion for higher-dimensional matrices.
        else:
            sign = -1
            det = 0
            for i in range(w):
                m = []
                for j in range(1, w):
                    b = []
                    for k in range(w):
                        if k != i:
                            b.append(self.matrix[j][k])
                    m.append(b)
                sign *= -1
                matrix = Matrix(m)
                det = det + mul * matrix.determinant(sign * self.matrix[0][i])
        return det

    # Use Cramer's method to solve Ax = b, where A = self.matrix and b is a vector.
    def cramers_method(self, b):

        # Check necessary preliminaries to make sure a solution exists.
        if len(b) != self.r:
            raise ValueError("b must have the same length as the matrix!")
        if self.invertible() == False:
            raise ValueError("There is no unique solution to this equation.")

        # Find the dterminant of the matrix for each column replaced with b.
        detMatrix = self.determinant()
        x_values = [None] * len(b)
        for i in range(self.r):
            tempMatrix = copy.deepcopy(self.matrix)
            for j in range(self.c):
                tempMatrix[j][i] = b[j]
            cramerMatrix = Matrix(tempMatrix)
            x_values[i] = cramerMatrix.determinant() / detMatrix
        return x_values

            
        cramerDet.append(det())                

    # Prints the matrix itself when called as a string.
    def __str__(self):
        return self.matrix

    # Prints the matrix itself when called.
    def __repr__(self):
       	return str(self.matrix)

    # Magic method to support addition
    def __add__(self, other):
        if isinstance(other, Matrix):
            return self.addition(other)
        elif isinstance(other, int):
            return self.add_constant(other)

    # Magic method to support subtraction
    def __sub__(self, other):
        if isinstance(other, Matrix):
            return self.subtraction(other)
        elif isinstance(other, int):
            return self.add_constant(-other)

    # Magic method to support constant multiplication
    def __mul__(self, other):
        if isinstance(other, int):
            return self.mul_constant(other)

    __rmul__ = __mul__
    __lmul__ = __mul__

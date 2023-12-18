from sympy import symbols, Matrix, solve, Poly, Sum , Array, permutedims, BlockMatrix

# Creating a class for Lie Algebra
class LieAlg:
    def __init__(self, basis, bracket_dict):
        self.basis=basis
        self.dict=bracket_dict

    # Define structure constant function
    def SC(self):
        basis=self.basis
        n = len(basis)
        # initialize structure constatns to zeros
        C = [[[0 for k in range(n)] for j in range(n)] for i in range(n)]
        for i in range(n):
            for j in range(n):
                # compute [e_i, e_j]
                expr = self.dict[(basis[i],basis[j])]
                for k in range(n):
                    # find the coeff of e_k in [e_i, e_j]
                    C[i][j][k]=expr.coeff(basis[k])
        return C
    # Define Ad matrix
    def Ad(self):
        basis=self.basis
        n= len(basis)
        # define standard basis vector
        # e = [Matrix([int(i == j) for j in range(n)]) for i in range(n)]
        
        # # Solve for T_{e_i}
        # # Define the matrix T_{e1} symbolically
        # A = Matrix(n, n, symbols('a:{}'.format(n * n)))
        # C = SC(basis, bracket_dict)
        # # express [e_i,e_j] as a sum
        # summation = [[Matrix.zeros(n,1) for i in range(n)] for j in range(n)]
        # for i in range(n):
        #     for j in range(n):
        #         for k in range(n):
        #             summation[i][j] +=C[i][j][k]*e[k]
        # # Define the equations A * e_i = [e1, ei] for each i
        # equations = [A * e[i] -  summation[0][i] for i in range(n)]
        # return solve(equations, A)

        # Make a list of matrices T_{ek} from SC(basis, bracket_dict)
        C = permutedims(self.SC(), (2,0,1))
        matrices = [Matrix.zeros(n,n) for k in range(n)]
        for i in range(n):
            matrices[i]= Matrix(C[i])
        return Matrix(BlockMatrix(matrices))
    # Define trivial d1 matrix
    def trivial(self):
        n = len(self.basis)
        matrices = [Matrix.zeros(n,n) for k in range(n)]
        for i  in range(n):
            matrices[i] = Matrix(C[i])
        return Matrix(BlockMatrix(matrices))
    
    # Define Center [g,g]
    def center(self):
        return self.Ad().nullspace()

    # Check if the input satisfies jacobi
    def check_jacobi_identity(self):
        n = len(self.basis)
        basis=self.basis
        C = self.SC()
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    for m in range(n):  # Added innermost loop for 'm'
                        sum_jacobi = 0
                        for l in range(n):  # Added loop for 'l'
                            sum_jacobi += C[i][j][l] * C[l][k][m] + C[j][k][l] * C[l][i][m] + C[k][i][l] * C[l][j][m]
                        if sum_jacobi != 0:
                            return False
        return True
    
    pass
# Consider Ad: g \to End(g)
#         x \mapsto Ad_x:= [x , \cdot ]
#     We have Ad(e_k)=k-th column of Ad, and note that Ad(x) is a n^2 vector which can be converted to a matrix T_x
#     Hence center(g)=null(Ad)={x \in g, Ad_x=0}
#------------------------------------------------------------------

def translate(e, basis, bracket_dict, a, b):
    if bracket_dict[(basis[a],basis[b])] == symbols('0'):
        return 0
    else:
        expr = bracket_dict[(basis[a],basis[b])]
        return expr.coeff(bracket_dict[(basis[a],basis[b])])*e[basis.index(bracket_dict[(basis[a],basis[b])])]
        



# test case
e1,e2 = symbols('e1 e2')
zero = symbols('0')
basis = [e1, e2]
n= len(basis)
e = [Matrix([int(i == j) for j in range(n)]) for i in range(n)]
dictionary = {
    (e1,e1): zero,
    (e1,e2): e2,
    (e2,e1): -1*(e2),
    (e2,e2): zero,
   
}
A=LieAlg(basis, dictionary)

C = A.SC()
# C1=permutedims(C, (2,0,1))
# matrices = [Matrix.zeros(n,n) for k in range(n)]
# for i in range(n):
#     matrices[i]= Matrix(C1[i])
# M = Matrix(BlockMatrix(matrices))
print(A.Ad)
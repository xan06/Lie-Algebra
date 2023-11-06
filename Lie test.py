from sympy import symbols, Matrix, solve, Poly


# Define structure constant function
def SC(basis, bracket_dict):
    n = len(basis)
    # initialize structure constatns to zeros
    C = [[[0 for k in range(n)] for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(n):
            # compute [e_i, e_j]
            expr = bracket_dict[(basis[i],basis[j])]
            for k in range(n):
                # find the coeff of e_k in [e_i, e_j]
                C[i][j][k]=expr.coeff(basis[k])
    return C

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
        
# Define Ad matrix
def Ad(basis, bracket_dict):
    n= len(basis)
    # define standard basis vector
    e = [Matrix([int(i == j) for j in range(n)]) for i in range(n)]
    
    # Solve for T_{e_i}
    # Define the matrix T_{e1} symbolically
    A = Matrix(n, n, symbols('a:{}'.format(n * n)))
    # Define the equations A * e_i = [e1, ei] for each i
    equations = [A * e[i] - translate(e, basis, bracket_dict, 0, i) for i in range(n)]
    return solve(equations, A)


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
a=translate(e,basis,dictionary,1,0)
print(a)

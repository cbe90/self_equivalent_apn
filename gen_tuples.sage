# this is to generate all the LE-automorphisms to consider for n = 6,7,8, listed in Corollary 1

# returns the dimension of Ord(A,i)
def check_order_space(A,i):
	n = 0
	V = VectorSpace(GF(2),A.nrows())
	for v in V:
		if (A^i*v == v):
			n = n+1
	return log(n,2)

# converts the matrix A to a look-up table
def matrix_to_sbox(A):
	T = []
	V = VectorSpace(GF(2),A.nrows())
	for v in V:
		T.append(ZZ(list(A*v),base=2))
	return T

# checks if a matrix tuple (A,B) can be a self invariance for a permutation, i.e., check whether, for all i, Ord(A,i) has the same dimension as Ord(B,i)
def is_permutation_tuple(A,B):
	for i in range(2**(A.nrows())):
		if (check_order_space(A,i) != check_order_space(B,i)):
			return False
	return True

# returns all invertible companion matrices of dimension d over GF(2)
def blocks_for_rcf(d):
	Q = PolynomialRing(GF(2),'X')
	R = []
	V = VectorSpace(GF(2),d+1)
	for v in V:
		if (v[0] == 1):
			if (v[d] == 1):
				R.append(companion_matrix(Q(list(v))))
	return R

# returns all rational canonical forms of GF(2)-matrices of dimension n
def get_rcfs(n):
	R = []
	for p in Partitions(n):
		V = []
		for i in p:
			V.append(blocks_for_rcf(i))
		CV = cartesian_product(V)
		for c in CV:
			app = True
			poly = c[0].minimal_polynomial()
			for i in range(len(c)):
				if (not c[i].minimal_polynomial().divides(poly)):
					app = False
				poly = c[i].minimal_polynomial()
			if (app == True):	
				R.append(block_diagonal_matrix(c[::-1]))
	return R


# checks if A=[A[0],A[1]] is power similar to B=[B[0],B[1]]
def is_power_similar(A,B):
	if ((A[0].multiplicative_order() == B[0].multiplicative_order()) and (A[1].multiplicative_order() == B[1].multiplicative_order())):
		for i in range(max(A[0].multiplicative_order(),A[1].multiplicative_order())):
			if ((A[0]**i).is_similar(B[0]) and (A[1]**i).is_similar(B[1])):
				return True
	return False

# checks if A=[A[0],A[1]] is extended power similar to B=[B[0],B[1]]
def is_extended_power_similar(A,B):
	if (is_power_similar(A,B)):
		return True
	if (is_power_similar([A[0].inverse(),A[1].inverse()],[B[1],B[0]])):
		return True
	return False

# generates all possible matrix tuples that need to be considered for self equivalence (i.e., 17 for n=6, 27 for n=7, and 32 for n=8)
def gen_permutation_tuples(n):
	T = get_rcfs(n)
	RC = []
	for t in T:                 
        	if (t.multiplicative_order().is_prime()):                       
			RC.append(t)
	T = tuples(RC,2)
	G = []
	for t in T:                
		if (t[0].multiplicative_order() == t[1].multiplicative_order()):
			G.append(list(t)+[1])
	GG = []
	ctr = 0                    
	for i in range(len(G)):                                           
		for j in range(len(G)):                         
			if not(i==j):              
				if(G[j][2]==1):  
					if(is_extended_power_similar(G[j],G[i])):                  
						G[i][2]=0                                     
						ctr = ctr+1  
						break      
		print(i,ctr)
	for g in G:
		if (g[2]==1):
			GG.append(g[0:2])
	GGG = []
	for i in range(len(GG)):  
		t = GG[i]                                                     
		print(i,len(GGG))                                         
		if (is_permutation_tuple(t[0],t[1])):
			GGG.append(t)
	return GGG

# generates the file AllTuples.h for the matrix tuples in list T
def generate_tuples_file(T):
	L = []
	for i in range(len(T)):
		L.append(matrix_to_sbox(T[i][0]) + matrix_to_sbox(T[i][1]))
	file = open("AllTuples.h", "w")
	file.write('#define N_TUPLES ' + repr(len(L)) + '\n\n')
	file.write('int AllTuples[N_TUPLES][' + repr(len(L[0])) + ']={\n')
	for ind in range(len(L)-1):
		file.write('{')
		for i in range(len(L[ind])-1):
			file.write(repr(L[ind][i]) + ', ')
		file.write(repr(L[ind][len(L[ind])-1]) + '},\n' )
	file.write('{')
	for i in range(len(L[len(L)-1])-1):
		file.write(repr(L[len(L)-1][i]) + ', ')
	file.write(repr(L[len(L)-1][len(L[len(L)-1])-1]) + '}\n' )
	file.write('};\n')
	file.close()

T6 = gen_permutation_tuples(6)
T7 = gen_permutation_tuples(7)
T8 = gen_permutation_tuples(8)


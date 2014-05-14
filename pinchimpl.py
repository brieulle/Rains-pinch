# Every review, remarks and corrections are welcomed  
# and to be sent to l.brieulle(at)gmail(dot)com


from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

def	pinch_method(p, n, m, f = None, g = None):
    '''
    Function that given three integers p, n and m such as :

    - p is prime,
    - m is dividing n and there exists a primitive mth root that spans F_{p^n} 
    and none of its subfields,

    returns the image, in a finite field G, of the polynomial generator of a 
    finite field F with the same cardinality as G.
    '''
    c, w = cputime(), walltime()
    R.<x> = PolynomialRing(GF(p))

    # If no polynomials are given, we compute both of them randomly.
    if f is None:
        f = R.irreducible_element(n, algorithm='random')
    if g is None:
        g = R.irreducible_element(n, algorithm='random')
    while f == g:
        g = R.irreductible_element(n, algorithm='random')

    # We compute two fields of cardinality p^n and two primitive m-rooth
    rootmf, rootmg, F, G = find_mroots_and_fields(p, n, m, f, g)

    # The matrixes will contain the coefficients of rootmf and rootmg in the 
    # basis x^i and y^i respectively
    A = matrix(GF(p), n, n)
    B = matrix(GF(p), n, n)

    for i in range(n):
    	A[i,:] = (rootmf**i).vector()

    # Failsafe
    try:
    	Ainv = A.inverse()
    except ZeroDivisionError:
    	print 'erreur'
        return A

    # We will try to find the power s such as phi(rootmf) = rootmg^s, since 
    # rootmg and rootmf are both primitive mrooth it is bound to happen if 
    # m is prime and statifies the conditions we gave above.
    s = 1	


    while s <= m :	
    	for i in range(n):
			B[i,:] = ((rootmg**s)**i).vector()


        # This will be the isomorphism's matrix
    	C = Ainv*B	

    	v = C[1,:]	  # The second line correponds to the image of x
        res = G(v[0])

        # I realized that you could try to find if the image of rootmf is also
        # a zero of the minimal polynomial of rootmg but it would force us to 
        # to compute yet another minimal polynomials. Instead, if you find that
        # the image is a root of the minimal polynomial of x, then you win!
        if f(res) == 0:
	    	print 'CPU %s, Wall %s' % (cputime(c), walltime(w))	
	    	return (res, C, rootmf, rootmg, s, f, F, G)

    	s = s + 1
	
    print 'No isomorphism found, check your m.'
    return 1
	
def find_mroots_and_fields(p, n, m, f, g):
    '''
    Computes explicitly two finite fields of cardinality of p^n and two 
    primitive m-rooths in both of them.
    '''    
    F.<a> = GF(p^n, modulus = f)
    G.<b> = GF(p^n, modulus = g)

    fact = m.factor()
    cofact = F.cardinality() // m


    # We look for m-rooth of order exactly m.
    rootmf = 1
    while any(rootmf**(m // k[0]) == 1 for k in fact):
        rootmf = F.random_element()**cofact
    rootmg = 1
    while any(rootmg**(m // k[0]) == 1 for k in fact):
        rootmg = G.random_element()**cofact

    return (rootmf, rootmg, F, G)

def calcul_img(elem, img_gen, F, G, mat = None): 
    '''
    Function that computes the image by an isomorphism of 
    an elem in F from the image of a generator or the matrix of said isomorphim.

    F, G : finite fields of the same cardinality,
    elem : the element of which we wish to find the image,
    img_gen : the image of F.gen()
    mat : the matrix of the isomorphism
    '''
    c, w = cputime(), walltime()
    p = F.characteristic()
    n = F.degree()
    res = 0
    
    if mat is None:
        for i in range(n):
            if elem[i] != 0:
                res = res + elem[i]*(img_gen)**i
    else:
        res = sum(G(C[i]) for i in range(n))

    print 'CPU %s, Wall %s' % (cputime(c), walltime(w))
    return res




    
    
    
    
    



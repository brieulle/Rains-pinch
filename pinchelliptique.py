'''
R.<X> = PolynomialRing(GF(5))

E = EllipticCurve(GF(5), [0,0,0,4,3])
f = X**11 + 3*X + 3
g = X**11 + X**5 + 4*X**3 + 4*X +2
m=23


F.<a> = GF(5**11, modulus = f)
G.<b> = GF(5**11, modulus = g)

d = a**3 + 4*a + 3
e = 2 + b + 4*b**2 + b**3

A = matrix(GF(5), 11, 11)
B = matrix(GF(5), 11, 11)

EFT = E.change_ring(F).quadratic_twist(d)
EGT = E.change_ring(G).quadratic_twist(e)
'''

def routine(p,n,m, E, f = None, g = None):
    c, w = cputime(), walltime()
    R = PolynomialRing(GF(p), 'X') 

    if f is None:
        f = R.irreducible_element(n, algorithm='random')
    if g is None:
        g = R.irreducible_element(n, algorithm='random')
    while f == g:
        g = R.irreducible_element(n, algorithm='random')
    
    F = GF(p**n, name='x', modulus = f)
    G = GF(p**n, name='y',modulus = g)

    EFT = E.change_ring(F)
    EGT = E.change_ring(G)

    OF = EFT.point([0,1,0])
    OG = EGT.point([0,1,0])

    cofact = EFT.cardinality()//m

    A = OF
    while any(i*A == OF for i in range(1,m)):
        A = cofact*EFT.random_point()

    B = OG
    while any(i*B == OG for i in range(1,m)):
        B = cofact*EGT.random_point() 

    M = matrix(GF(p), n, n)
    N = matrix(GF(p), n, n)

    for i in range(n):
        M[i,:] = (A[0]**i).vector()

    try:
        Minv = M.inverse()
    except ZeroDivisionError:
        print 'booooooooooooo'
        return [M,A,B]

    facteur = 1
    while facteur < m:
        for i in range(n):
            N[i,:] = (((facteur*B)[0])**i).vector()

        C = Minv*N
        v = C[1,:]
        res = G(v[0])
        tab = []
        tab.append(res)
         
        if f(res) == 0:
            print 'CPU %s, Wall %s' % (cputime(c), walltime(w))
            return res
        facteur = facteur + 1
        
    print 'boohoo...'
    return [res, M, N, C, facteur, tab]

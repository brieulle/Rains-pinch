# -*- coding: utf-8 -*-

def isom_elliptic(k1, k2, k = None, Y_coordinates = False, bound = None):
    '''
    INPUT : 
    a base field k, 
    two extensions of the same degree k1, k2,
    a boolean Y_coordinates. 

    OUTPUT : a tuple with normal elements in k1 and k2 respectively.

    Use :

    isom_elliptic(k1, k2, k)

    returns two normal elements using abscissas of torsion point. If you want to
    use the ordinates, write 
    
    isom_elliptic(k, k1, k2, True)

    Algorithm:

    Given two extensions of the same base field and the same degree, we return 
    two normal elements the use of Gaussian elliptic period, via the function 
    find_unique_orbit_elliptic, on an curve E which is determined by the 
    function find_elliptic_curve.

    First we have to find an integer m and an elliptic curve E  over k such 
    that :

    - m divides the number of points of E over the extensions k1 and k2.

    - if m is power of p, we want that the trace t of the Frobenius on E over k 
      is of order n, the degree of the extension in (Z/m)*.

    - if m is a prime power, we want that the root of smallest order in (Z/m)*
      of the characteristic polynomial of the Frobenius over k is of order n,
      the degree of the extension, i.e.:

      For x**2 -tx + q = (x - a)(x - b) mod m, we want that ord_m(a) = n nd
      ord_m(a) < ord_m(b).

    - if m is composite, there is no algorithm yet.

    Once we have both of them, we have to compute two m torsion points P and Q 
    of order m in E/k1 and E/k2 respectively and then we compute their Gaussian
    elliptic periods with abscissas or ordinates depending on the boolean given 
    in arguments; that is done in find_unique_orbit_elliptic.
    
    Thus we found two unique elements alpha and beta such that there's exist an 
    isomorphism phi with :

    phi(alpha) = beta

    and you can find phi from there using the function convert as done in the 
    function isom_normal.
    '''
    c, w = cputime(), walltime()
    if k is None:
	    k = k1.base_ring()
    p = k.characteristic()
    n = k1.degree()  
    q = k1.cardinality()
    
    # We take the smallest m, we can change that later.
    m = find_m(n, bound)[0]
    
    # Finding the elliptic curve on which we can work. 
    E_s= find_elliptic_curve(k, k1, m) 
    G = find_group_gen(n, m, E_s[1])

    Ek1 = E_s[0].change_ring(k1)
    Ek2 = E_s[0].change_ring(k2)

    a, b = (find_unique_orbit_elliptic(Ek1, G, m, Y_coordinates), 
    find_unique_orbit_elliptic(Ek2, G, m, Y_coordinates))
    print 'CPU %s, Wall %s' % (cputime(c), walltime(w))

    return a, b

def find_unique_orbit_elliptic(E, G, m, Y_coordinates = False):
    '''
    INPUT : 
    
    an elliptic curve E, a tuple G, an intger m, a boolean Y_coordinates

    OUPUT : 
    
    a Gaussian elliptic periods of point P of order m in E/k, k a finite field

    Use : 
    find_unique_orbit_elliptic(E,G)

    returns a Gaussian elliptic periods using the X-coordinates of a point P of 
    order m. If you want to use the Y-coordinates, type :

    find_unique_orbit_elliptic(E,G,True)

    Algorithm:

    Function that given an elliptic curve and a group G returns the following 
    Gaussian elliptic periods :

    sum_{\sigma\in G}{([\sigma]P)_{T}}

    with T = X or Y depending on the boolean given in arguments.
    '''
    cofactor = E.cardinality()//m
    sum_P = []

    # Searching for a point of order exactly m.
    P = E(0)
    while any((m//i)*P == 0 for i in m.prime_divisors()):
        P = cofactor*E.random_point()

    for i in range(G[1]):
            sum_P.append(ZZ(G[0]**i)*P)

    if not Y_coordinates:
        return sum(P[0] for P in sum_P)
    else:
        raise NotImplementedError, 'No algorithm for Y-coordinates yet.'


def find_elliptic_curve(k, K, m):
    '''
    INPUT : 

    a base field k, an extension K, a integer m 

    OUTPUT : 
    
    an elliptic curve over k, a trace and eventually a tuple containing a
    element of order n in (Z/m)* and the class of the trace modulo m.

    Algorithm :

    Function that finds an elliptic curve with the required charateristics, 
    those given in the function isom_elliptic.

    First, we have to determine if m is composite, a prime power or a power of 
    p, the characteristic of the base field. The first case is not implemented 
    yet. 
    We also note that the m given should satisfies several conditions based
    on the characteristic and the degree of K. See the docstrings of 
    isom_elliptic for more information.

    - If m is a power of p, the charateristic of the base field k, then we shall
      proceed as follow :

        We pick a random curve E/k and we set down t = Tr_k(Fr_E), for the
        curve to be what we want, we need :

          - (Z/m)* = <t> x S or #<t> = n (or something)
          - m divides #E/K but not #E/L, for any intermediate extension L 
            of K/k; so we can construct points of order m such that their 
            abscissas or ordinates span exactly K. Or that we haven't any point
            of order m in any sub-extension.

        Then we test those conditions for both E and its quadratic twist, if one
        of them meet the requirements, we return it and its trace. If no 
        elliptic curves are found an error is returned.

    - If m is primer power, then we shall proceed as follow :

        We have m = l^r, for l a prime. For this method to work, we need l to 
        be an Elkies prime. A prime l is an Elkies prime for an elliptic curve 
        if the charateristic polynomial of the aforesaid elliptic curve splits 
        in GF(l).

        For now, we pick a random curve E/k and for it to work, if we set down 
        t = Tr_k(Fr_E), we need the following :

          - We have x**2 - tx + q = (x - a)(x - b) mod m, meaning the polynomial
            splits in Z/m,
          - (Z/m)* = <a> x S, with #<a> = n, 
          - ord_m(a) < ord_m(b),
          - m divides #E/K but not #E/L, for any intermediate extension L 
            of K/k; so we can construct points of order m such that their 
            abscissas or ordinates span exactly K.

        Once again, we test all that for both E and its quadratic twist; if one 
        them meet the requirements, we return it, its trace and a tuple
        containing the root a and t mod m. If none are found, there is 
        something wrong.


    - If m is composite, TODO.
    '''
    def test_curve(E, t, S, m_case):
        '''
        INPUT : An elliptic curve E, an integer t, a list S and an integer
        m_case

        OUTPUT : a boolean or a tuple with a bolean and an integer

        Function that is passed down to find_elliptic. It deterines if the curve
        meets the desired requirements depending on the nature of m : a power of
        p(1), a prime power(2) or a composite number(3).
        '''
        Zm = Zmod(m)

        # m = p^e
        if m_case == 0:
            # If the trace is none of the candidate for trace of order n in
            # (Z/m)*, then we don't want this curve.
            if all(Zm(t) != t_m for t_m in S):
                return False
            # If we can't find point of order m, we don't want this curve either
            elif E.change_ring(K).cardinality()%m != 0:
                return False
            # We want the point of order m to span exactly K/k and not any
            # sub-extension.
            elif any(E.change_ring(k.extension(n//d)).cardinality()%m == 0
                    for d in n.prime_divisors()):
                return False
            else:
                return True
        elif m_case == 1:
            # We're trying to find if t mod m is equal to one of the trace in
            # S (a list of tuple). Then we want to remember the index for
            # which it is right because we will need the root associated to
            # compute the Galois group.
            index = None
            
            for i in range(len(S)):
                if Zm(t) == S[i][1]:
                    index = i
                    break
            if index is None:
                return (False, None)
            elif E.change_ring(K).cardinality()%m != 0:
                return (False, None)
            elif any(E.change_ring(k.extension(n//d)).cardinality()%m == 0
                    for d in n.prime_divisors()):
                return (False, None)
            else:
            	return (True, index)
        elif m_case == 2:
            raise NotImplementedError, 'm composite is not implemented yet'

    p = k.characteristic()
    q = K.cardinality()
    n = K.degree()

    if not m.is_prime_power():
        raise NotImplementedError, 'Case m composite is not implemened yet.'
    else:
        if m%p == 0:
            # Picking the candidates class modulo m
            S_t = find_trace(n, m, k)
            E_rejected = []

            while True:
                # This method is far from optimal, but we assume that after q 
                # draws we have a good chance of trying enough curves.
                if len(E_rejected) > q:
                    raise RuntimeError, 'No suitable elliptic curves found.'

                E = EllipticCurve(j = k.random_element())
                while any(E == Ef for Ef in E_found):
                    E = EllipticCurve(j = k.random_element())

                t = E.trace_of_frobenius()

                # We want an ordinary curve. More precisely, we can't find an
                # order for a zero element.
	        if t%p == 0:
	            continue

                # We try to see if E or its quadratic twist meets the 
                # requirements
                for EE,tt in [(E,t), (E.quadratic_twist(), -t)]:
                    if test_curve(EE, tt, S_t, 0):
                        return EE, tt

                # We don't want to work on those curves anymore.
                E_rejected.append(E)
                E_rejected.append(E.quadratic_twist())

        else:
            S_at = find_trace(n,m,k)
            E_rejected = []

            while True:
                if len(E_rejected) > q:
                    raise RuntimeError, 'No suitable elliptic curves found.'

                E = EllipticCurve(j = k.random_element())
                while any(E == Ef for Ef in E_rejected):
                    E = EllipticCurve(j = k.random_element())

                t = E.trace_of_frobenius()

                for EE,tt in [(E,t), (E.quadratic_twist(), -t)]:
                    res = test_curve(EE, tt, S_at, 1)
                    if res[0]:
                        # Are we only interested in the class of a mod m ?
                        return EE, S_at[res[1]][0]

                E_rejected.append(E)
                E_rejected.append(E.quadratic_twist())
                    

def find_trace(n,m,k):
    '''
    INPUT : an integer n, an integer m, a base field k

    OUTPUT : a list of integer mod m or a list of a couple of integers mod m

    Function that gives a list of candidates for the trace.
    It returns a list of trace of order n in (Z/m)* if m%p = 0
    and a list of couple (a, t) when m is a prime power; where 
    a is the root of X² -tX + q of smallest order equal to n in
    (Z/m)*.
    We could possibly just return the trace in the second case.
    But I don't know yet how the group (Z/m)*/<a> = S will be 
    implemented.
    '''
    Zm = Zmod(m)
    p = k.characteristic()
    q = k.cardinality()

    # If m is a multiple of p, then we just need the trace to be of order 
    #exactly n in (Z/m)*
    if m%p == 0:
        sol = []
        for t in Zm:
            # We only want the trace to be of order exactly n in (Z/m)* and not 
            # to define supersingular curves.
            if not t.is_unit():
                continue
            elif (Zm(t).multiplicative_order() != n):
                continue
            else:
                sol.append(t)
        return sol
    # If m is prime (power), then we need to look at the roots
    # Probably a temporary condition
    elif m.is_prime_power():
        sol = []
        for a in Zm:
            # If a is not invertible in Z/m, we can't compute any order.
            # We'll probably have to look for an element of order exactlu n
            # We don't want a = 0 mod m. 
            if not a.is_unit():
                continue
            else:
                ord_a = Zm(a).multiplicative_order()
                ord_b = Zm(q/a).multiplicative_order()
                # We need an element of order n
                if (ord_a != n and ord_b !=n):
                    continue
                # We want ord_a != ord_b
                elif (ord_a == ord_b):
                    continue
                elif (ord_b != n): 
                    if (ord_a > ord_b):
                        continue
                    else:
                        sol.append((a, a + q/a)) #return (a, a + q/a, q, q/a)
                elif (ord_a != n):
                    if (ord_b > ord_a):
                        continue
                    else:
                        sol.append((q/a, a + q/a))
        return sol
        
    else:
        raise NotImplementedError, 'm composite is not implemented yet'

def find_m(degree, q, bound = None):
    '''
    INPUT : three integers degree & q & bound

    OUTPUT : an integer m

    Algorithm :

    Functions that given an integer n (degree of an extension) and a bound 
    returns all the candidates m, such that :

    - n|phi(m), the euler totient function,

    - (n, phi(m)/n) = 1,

    - Another one ? Maybe for q = p^d we'd want (n,d) = 1,

    - m <= bound, the bound is at most q^n + 2*sqrt(q^n) + 1 since we wish for 
    it to divide E(F_q^n). Also, this generally too big to handle.

    We can note that if m = r^e with (e-1,n) = 1 or e = 1, then r = a*n + 1 with
    (a,n) = 1 is a suitable form for m as then phi(m) = (a*n)(an + 1)^(e-1);
    n|phi(m) and n doesn't divide (an+1)^(e-1), it's even coprime with it (not 
    so sure about that, it'd be nice; or we could add that as a 
    condition if need be).
    It also works in the general case if all the prime factors of m are of the 
    form a*n + 1 with (a,n) = 1. You just have to apply that to them and 
    multiply the results.
    '''
    if bound is None:
        bound_a = 100  # Arbitrary value.  
    else:
        # if m = a*n + 1 < b, then a < (b- 1)/n.
        bound_a = (bound - 1) / degree 

    sol = []

    for a in range(bound_a):
        m = a*degree + 1
        if not m.is_prime_power():
            continue 
        else:
            if (euler_phi(m)//degree).gcd(degree) == 1:
                sol.append(m) 
        
    return sol


def find_group_gen(n, m, a_t):
    '''
    INPUT : an integer p, characteristic of the ambient base field,
        an integer q, the cardinality of the ambient base field,
        an integer n, the degree of the considered extension.

    OUTPUT : a generator of G and its order

    Algorithm :

    For now, we still focus solely on the case m a prime power. The 
    situation is as follow :

    We have the cyclic finite groupe (Z/m)* and a subgroup of order n 
    <alpha> or <t>, depending on the situation. We are trying to determine 
    the generator of the quotient subgroup (Z/m)*/<alpha>. It shall be of order 
    phi(m)/n.
    '''
    gen = Zmod(m).unit_gens()
    sol = []

    if not m.is_prime_power():
        raise NotImplementedError, 'm composite is not implemented yet.'
    else:
        order = euler_phi(m) / n
        sub_group = [a_t**i for i in range(n)]

        for i in range(euler_phi(m)):
            if any(gen[0]**i == s for s in sub_group):
                continue
            else:
                if any(((gen[0]**i)*s == sol[k] for s in sub_group)
                        for k in range(len(sol))):
                    continue
                else:
                    sol.append(gen[0]**i)# ord

        return sol[0], order #raise RuntimeError, 'No generator found.'


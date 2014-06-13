# -*- coding: utf-8 -*-
from rains import *

def isom_elliptic(k1, k2, k = None, Y_coordinates = False):
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

    First we have to find an integer m and  an elliptic curve E  over k such 
    that :

    - m divides the number of points of E over the extensions k1 and k2.

    - if m is power of p, we want that the trace t of the Frobenius on E over k 
      is of order n, the degree of the extension, in (Z/m)*.

    - if m is a prime power, we want that the root of smallest order in (Z/m)*
      of the characteristic polynomial of the Frobenius over k is of order n,
      the degree of the extension, i.e.:

      For x**2 -tx + q = (x - a)(x - b) mod m, we want that ord_m(a) = n.

    - if m is composite, there is no algorithm yet.

    This is the work of both functions find_root_order and find_elliptic_curve.

    Once we have both of these, we have to compute two m torsion points P and Q 
    of order m in E/k1 and E/k2 respectively and then we compute their Gaussian
    elliptic periods with abscissas or ordinates depending on the boolean given 
    in arguments; that is done in find_unique_orbit_elliptic.
    
    Thus we found two elements alpha and beta such that there's exist an 
    isomorphism phi with :

    phi(alpha) = beta

    and you can find phi from there using the function conver as done in the 
    function isom_normal.
    '''
    c, w = cputime(), walltime()
    if k is None:
	k = k1.base_ring()
    p = k.characteristic()
    n = k1.degree()  
    q = k1.cardinality()
    
    m = function_that_finds_m()
    
    # Finding the elliptic curve on which we can work. 
    E = find_elliptic_curve(k, k1, m)

    G = function_that_finds_G_from_m_and_something_else_maybe()

    Ek1 = E.change_ring(k1)
    Ek2 = E.change_ring(k2)

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

    # Searching for a point of order exactly m.
    P = E(0)
    while any((m//i)*P == 0 for i in m.prime_divisors()):
        P = cofactor*E.random_point()

    if not Y_coordinates:
        # Return the sum of ([a]P)_X for a in G.
        return sum((prod(ZZ(g)**e for (g, _), e in zip(G, exps))*P)[0]
            for exps in CProd(*map(lambda (_,x): xrange(x), G)))
    else:
        return sum((prod(ZZ(g)**e for (g, _), e in zip(G, exps))*P)[1]
            for exps in CProd(*map(lambda (_,x): xrange(x), G)))


def find_elliptic_curve(k, K, m):
    '''
    INPUT : 

    a base field k, an extension K, a integer m 

    OUTPUT : 
    
    an elliptic curve over k

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

        We still pick a random curve E/k and we set down t = Tr_k(Fr_E), for the
        curve to be what we want, we need :

          - (Z/m)* = <t> x S or #<t> = n (or something)
          - m divides #E/K but not #E/L, for any intermediate extension L 
            of K/k; so we can construct points of order m such that their 
            abscissas or ordinates span exactly K.

        Then we test those conditions for both E and its quadratic twist, if one
        of hem meet the requirements, we return it. If no elliptic curves are 
        found an error is returned.

    - If m is primer power, then we shall proceed as follow :

        We have m = l^r, for l a prime. For this method to work, we need l to 
        be an Elkies prime. A prime l is an Elkies prime for an elliptic curve 
        if the charateristic polynomial of the aforesaid elliptic curve splits 
        in GF(l).

        For now, we pick a random curve E/k and for it to work, if we set down 
        t = Tr_k(Fr_E), we need the following :

          - We have x**2 - tx + q = (x - a)(x - b) mod m, meaning the polynomial
            splits in Z/m,
          - (Z/m)* = <a> x S, with <a> the Galois group of the extension K/k or
            equivalently #<a> = n, 
            (!!! IMPORTANT !!! : I'm not yet sure about the second part of 
            this statement)
          - ord_m(a) < ord_m(b),
          - m divides #E/K but not #E/L, for any intermediate extension L 
            of K/k; so we can construct points of order m such that their 
            abscissas or ordinates span exactly K.

        Once again, we test all that for both E and its quadratic twist; if one 
        them meet the requirements, we return it. If none are found, there is 
        something wrong.


    - If m is composite, TODO.
    '''
    def test_curve(E, t, S, m_case):
        '''
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
            elif any(E.change_ring(k.extension(n//d)).cardinality()%m != 0
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
            elif any(E.change_ring(k.extension(n//d)).cardinality()%m != 0
                    for d in n.prime_divisors()):
                return (True, index)
        elif m_case == 2:
            raise NotImplementedError, 'm composite is not implemented yet'

    p = k.characteristic()
    q = K.cardinality()
    n = K.degree()

    if not m.is_prime_power():
        raise NotImplementedError, 'Case m composite is not implemened yet.'
    else:
        # This method is far from optimal, but we assume that after q draws we 
        # have a good chance of trying enough curves.

        if m%p == 0:
            # Picking the candidates class modulo m
            S_t = find_trace(n, m, k)
            E_rejected = []

            while True:
                if len(E_rejected) > q:
                    raise RuntimeError, 'No suitable elliptic curves found.'

                E = EllipticCurve(j = k.random_element())
                while any(E == Ef for Ef in E_found):
                    E = EllipticCurve(j = k.random_element())

                t = E.trace_of_frobenius()

                # We want an ordinary curve. More precisely, if t = 0, we can't 
                # compute its order in (Z/m)*
	        if t%p == 0:
	            continue

                # We try to see if E or its quadratic twist meets the 
                # requirements
                for EE,tt in [(E,t), (E.quadratic_twist(), -t)]:
                    if test_curve(EE, tt, S_t, 0):
                        return EE, tt

                # We don't want to work on those curve anymore.
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
                        # Are we only interested in the class of t mod m ?
                        return EE, tt, S_at[res[1]]

                E_rejected.append(E)
                E_rejected.append(E.quadratic_twist())
                    

def find_trace(n,m,k):
    '''
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

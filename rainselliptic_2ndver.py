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
    if k is None:
	    k = k1.base_ring()
    p = k.characteristic()
    n = k1.degree()  
    q = k1.cardinality()
    
    # We compute a list of candidates for m (i.e. such that n divides phi(m) 
    # and (phi(m)/n,n) = 1. It lacks the conditions on the trace.
    m_f = find_m(n, bound)
    
    # Finding the elliptic curve on which we can work. 
    E_s = find_elliptic_curve(k, k1, m_f) 

    Ek1 = E_s[0].change_ring(k1)
    Ek2 = E_s[0].change_ring(k2)

    a, b = (find_unique_orbit_elliptic(Ek1, E_s[1], Y_coordinates), 
    find_unique_orbit_elliptic(Ek2, E_s[1], Y_coordinates))

    return a, b

def find_unique_orbit_elliptic(E, m, Y_coordinates = False):
    '''
    INPUT : 
    
    an elliptic curve E, an integer m, a boolean Y_coordinates

    OUPUT : 
    
    a Gaussian elliptic periods of point P of order m in E/K, K a finite field

    Use : 
    find_unique_orbit_elliptic(E, m)

    returns a Gaussian elliptic periods using the X-coordinates of a point P of 
    order m. If you want to use the Y-coordinates, type :

    find_unique_orbit_elliptic(E, m, True)

    Algorithm:

    Function that given an elliptic curve and a group G returns the following 
    Gaussian elliptic periods :

    sum_{\sigma\in G}{([\sigma]P)_{T}}

    with T = X or Y depending on the boolean given in arguments.
    '''
    cofactor = E.cardinality()//m
    n = E.base_ring().degree()
    sum_P = []
    order = euler_phi(m)//n
    if not m.is_prime_power():
        raise NotImplementedError, 'case m composite not implemented yet.' 
    else:
        gen_G = Zmod(m).unit_gens()[0]

        # Searching for a point of order exactly m.
        P = E(0)
        while any((m//i)*P == 0 for i in m.prime_divisors()):
            P = cofactor*E.random_point()

        for i in range(order):
            sum_P.append(ZZ(gen_G**i)*P)

        if not Y_coordinates:
            return sum(P[0] for P in sum_P)
        else:
            raise NotImplementedError


def find_elliptic_curve(k, K, m_f):
    '''
    INPUT : 

    a base field k, an extension K, a list integer m_f

    OUTPUT : 
    
    an elliptic curve over k, an integer m 

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

          - t not zero,
          - (Z/m)* = <t> x S or #<t> = n 
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
    def _test_curve(E, t, S):
        '''
        INPUT : An elliptic curve E, an integer t, a list S of class modulo m 

        OUTPUT : a boolean

        Algorithm :

        Function that determines if a curve E has its trace t = t_m mod m for 
        t_m in S; plus that m divides only the number of point of E/K.
        '''
        Zm = Zmod(m)

        # If the trace is none of the candidate for trace of order n in
        # (Z/m)*, then we don't want this curve.
        if all(Zm(t) != t_m for t_m in S):
            return False
        else:
            return True

    p = k.characteristic()
    q = K.cardinality()
    n = K.degree()
    m = 0

    # We look for the proper m, i.e. the one such that we can find an eigenvalue
    # a such that ord_m(a) = n < ord_m(q/a) or the other way around, or a trace
    # such that ord(t) = n
    for M in m_f:
        S_t = find_trace(n, M, k)
        if len(S_t) != 0:
            m = M
            break

    if m == 0:
        raise RuntimeError, 'No proper m found.'


    if not m.is_prime_power():
        raise NotImplementedError, 'Case m composite is not implemened yet.'
    elif m%p == 0:
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
                E_rejected.append(E)
                E_rejected.append(E.quadratic_twist())
                continue

            # We try to see if E or its quadratic twist meets the 
            # requirements
            for EE,tt in [(E,t), (E.quadratic_twist(), -t)]:
                if _test_curve(EE, tt, S_t):
                    print 'Temps pour trouver E : CPU %s, Wall %s' %(
                                            cputime(c), walltime(w))
                    return EE, m

            # We don't want to work on those curves anymore.
            E_rejected.append(E)
            E_rejected.append(E.quadratic_twist())

    else:
        E_rejected = []

        while True:
            if len(E_rejected) > q:
                raise RuntimeError, 'No suitable elliptic curves found.'

            E = EllipticCurve(j = k.random_element())
            while any(E == Ef for Ef in E_rejected):
                E = EllipticCurve(j = k.random_element())

            t = E.trace_of_frobenius()

            for EE,tt in [(E,t), (E.quadratic_twist(), -t)]:
                if _test_curve(EE, tt, S_t):
                    return EE, m

            E_rejected.append(E)
            E_rejected.append(E.quadratic_twist())
                    

def find_trace(n,m,k):
    '''
    INPUT : an integer n, an integer m, a base field k

    OUTPUT : a list of integer mod m or a list of a couple of integers mod m

    Algorithm :

    If m is a power of p, then we look for class modulo m with order equal to n.
    Then, we return the list of all such class.

    If m is a power of prime different from p, we look for a in (Z/m)* such 
    that :

    - ord_m(a) < ord_m(q/a) and ord_m(a) = n,

    or

    - ord_m(q/a) < ord_a and ord_m(q/a) = n.

    And we return a + q/a.

    Here a plays the role of one of the two roots of the future characteristic 
    polynomial of the Frobenius of the elliptic curve we'll use; i.e.

    X^2 - (a + q/a)*X + a*(q/a) = X^2 - t*X + q

    if we write t = a + q/a. From that, we will pick elliptic curves which have 
    one of the t's as trace of its Frobenius.
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
            if not a.is_unit():
                continue
            else:
                ord_a = Zm(a).multiplicative_order()
                ord_b = Zm(q/a).multiplicative_order()
                # We need an element of order n
                if (ord_a == ord_b):
                    continue
                elif min(ord_a, ord_b) != n:
                    continue
                else:
                    sol.append(a + q/a)
        return sol
        
    else:
        raise NotImplementedError, 'm composite is not implemented yet'

def find_m(degree, bound = None):
    '''
    INPUT : an integers degree, an integer bound

    OUTPUT : a list of integer

    Algorithm :

    Functions that given an integer n (degree of an extension) and a bound 
    returns all the candidates m, such that :

    - n|phi(m), the euler totient function,

    - (n, phi(m)/n) = 1,

    - Another one ? Maybe for q = p^d we'd want (n,d) = 1,

    We can note that if m = r^e with (e-1,n) = 1 or e = 1, then r = a*n + 1 with
    (a,n) = 1 is a suitable form for m as then phi(m) = (a*n)(an + 1)^(e-1);

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
        elif (euler_phi(m)//degree).gcd(degree) != 1:
            continue
        else:
            sol.append(m) 
        
    return sol

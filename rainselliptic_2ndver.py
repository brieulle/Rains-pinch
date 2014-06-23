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
    #c_comp, w_comp = cputime(), walltime()
    #c, w = cputime(), walltime()
    m_t = find_m(n, k, bound)
    #print 'Computing the m\'s : CPU %s, Walltime %s' % (cputime(c), walltime(w))
    
    if m_t is None:
        raise RuntimeError, 'No suitable m found, increase your bound.'

    # Finding the elliptic curve on which we can work. 
    #c, w = cputime(), walltime()
    E, case, m = find_elliptic_curve(k, k1, m_t) 
    #print 'Finding E : CPU %s, Walltime %s' % (cputime(c), walltime(w))

    if E is None:
        raise RuntimeError, 'No suitable elliptic curve found, check your \
                                                                    parameters'

    #c, w = cputime(), walltime()
    Ek1 = E.change_ring(k1)
    Ek2 = E.change_ring(k2)
    #print 'Computing E in extension : CPU %s, Walltime %s' % (cputime(c),
    #                                                             walltime(w))

    #c, w = cputime(), walltime()
    a, b = (find_unique_orbit_elliptic(Ek1, m, Y_coordinates, case),
    find_unique_orbit_elliptic(Ek2, m, Y_coordinates, case))
    #print 'Computing periods : CPU %s, Walltime %s' % (cputime(c), walltime(w))
    #print 'Total time : CPU %s, Walltime %s' % (cputime(c_comp),
    #                                                         walltime(w_comp))
    #if not a or not b:
    #    return m_t, E

    return a, b

def find_unique_orbit_elliptic(E, m, Y_coordinates = False, case = None):
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

    case = None => default case : j != 0, 1728 or t = 0
    case = 1 => j = 1728 and t != 0
    case = 2 => j = 0 and t != 0

    Algorithm:

    Function that given an elliptic curve and a group G returns the following 
    Gaussian elliptic periods :

    sum_{\sigma\in G}{([\sigma]P)_{T}}

    with T = X or Y depending on the boolean given in arguments.
    '''
    cofactor = E.cardinality()//m
    n = E.base_ring().degree()

    # Searching for a point of order exactly m.
    P = E(0)
    while any((m//i)*P == 0 for i in m.prime_divisors()):
        P = cofactor*E.random_point()

    if not m.is_prime_power():
        raise NotImplementedError, 'case m composite not implemented yet.' 
    elif case is None:
        # Looking for a generator of order exactly phi(m)/n in 
        # phi(m)/something.
        gen_G = Zmod(m).unit_gens()[0]**n
        # Note : Ã  priori, l'ordre peut Ãªtre le mÃªme en utilisant les X ou 
        # les Y. Donc il faudrait quotienter par {+-1} mÃªme pour les Y. 
        # Sauf si la puissance carrÃ© dans la somme Ã©quivaut au final Ã  
        # quotienter par i.
        order = euler_phi(m)//(2*n)

        if not Y_coordinates:
            return sum((ZZ(gen_G**i)*P)[0] for i in range(order))
        else:
            return sum(((ZZ(gen_G**i)*P)[1])**2 for i in range(order))
    elif case == 1:
        gen_G = Zmod(m).unit_gens()[0]**n
        order = euler_phi(m)//(4*n)
        
        if not Y_coordinates:
            return sum(((ZZ(gen_G**i)*P)[1])**2 for i in range(order))
        else:
            raise NotImplementedError

    elif case == 2:
        if not Y_coordinates:
            raise NotImplementedError
        else:
            raise NotImplementedError



def find_elliptic_curve(k, K, m_t):
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
    p = k.characteristic()
    q = K.cardinality()
    n = K.degree()

    # a such that ord_m(a) = n < ord_m(q/a) or the other way around, or a trace
    # such that ord(t) = n

    j_rejected = []

    while True:
        # This method is far from optimal, but we assume that after q 
        # draws we have a good chance of trying enough curves.
        if len(j_rejected) > q:
            raise RuntimeError, 'No suitable elliptic curves found.'

        j = k.random_element()
        E = EllipticCurve(j = j)
        while any(j == jr for jr in j_rejected):
            j = k.random_element()
            E = EllipticCurve(j = j)

        t = E.trace_of_frobenius()

        # We try to see if E or its quadratic twist meets the 
        # requirements
        if j == 1728 :

            '''
            if t == 0:
                j_rejected.append(j)
                continue

            # We need a m with the additionnal property that 4 divides 
            # phi(m)
            for i in range(len(m_t)):
                if not m_t[i][0]%4:
                    m = None
                    continue
                else:
                    m = m_t[i][0]
                    S_t = m_t[i][1] 

            # TODO : A test on m to see if it's a prime power or not
            L = [(E,t), (E.quadratic_twist(), -t), (E.quartic_twist(?), ?t),
                    (E.quartic_twist(?), ?t)] 
            for EE, tt in L:
                if (Zmod(m)(tt) in S_t):
                    return EE, 1, m
            
            j_rejected.append(j)
            '''

            raise NotImplementedError
                
        elif j == 0 :
            raise NotImplementedError
        else:
            # There's no additional requirements for m in this case, so 
            # we just take the first one.
            m = m_t[0][0]
            S_t = m_t[0][1]

            # TODO : A test on m to see if it's a prime power or not
            L = [(E,t), (E.quadratic_twist(), -t)]
            for EE,tt in L:
                if (Zmod(m)(tt) in S_t):
                    print 'grrrrrrrrrrr'
                    return EE, None, m

        # We don't want to work on those curves anymore.
        j_rejected.append(j)

    # No elliptic curve found. At least, not in q random draws.
    return None

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
    if not m.is_prime_power():
        raise NotImplementedError
    elif m%p == 0:
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
    else:
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
        return set(sol)

def find_m(n, k, bound = None):
    '''
    INPUT : an integers n, a base field k, an integer bound

    OUTPUT : an integer

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
        bound_a = 100 # Arbitrary value.  
    else:
        # if m = a*n + 1 < b, then a < (b- 1)/n.
        bound_a = (bound - 1) / n 

    sol = []

    for a in range(bound_a):
        m = a*n + 1
        # m composite not implemented yet
        if not m.is_prime_power():
            continue 
        elif (euler_phi(m)//n).gcd(n) != 1:
            continue
        else:
            S_t = find_trace(n, m, k)
            if len(S_t) < 1:   # Some time in the future we'd like to have a 
                continue       # better bound than just 1.
            else:
                sol.append((m, S_t))

    if sol == []:
        return None
    else:
        return sol

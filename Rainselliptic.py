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
    
    isom_elliptic(k, k1, k2, true)

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

    # Computing the quotient group of the Galois group with which we will 
    # compute the Gaussian elliptic periods. It will also gives the m we 
    # need for the method to work.
    # See Luca De Feo's page for more informations :
    # https://www.github.com/defeo/ffisom/blob.master/rains.py
    G = find_root_order(p,n, accept_elliptic)[1]
    m = G[0][0].parent().order()
    
    # Finding the elliptic curve on which we can work. 
    E = find_elliptic_curve(k, k1, m)

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

    find_unique_orbit_elliptic(E,G,true)

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

    a base field k, an extension K, a integer m, a integer n, degree of 
    the extension

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
    counter = 0
    p = k.characteristic()
    q = K.cardinality()
    n = K.degree()
    Zm = Integers(m)
    #f = m.factor()

    if not m.is_prime_power():
        raise NotImplementedError, 'Case m composite not implemened yet.'
    else:
        # This method is far from optimal, but we assume that after q draws we 
        # have a good chance of trying enough curves.
        while counter < q:
            # We pick a random elliptic curve and we'll try to see if it 
            # corresponds to what we want (i.e. trace of E = a good t).
            E = EllipticCurve(j = k.random_element())
            t = E.trace_of_frobenius()
            #EK = E.change_ring(K)

            if m%p == 0:
                # We want an ordinary curve. More precisely, if t = 0, we can't 
                # compute its order in (Z/m)*
				if not (t%p):
					counter = counter + 1
					continue

				sucess = false
                # We try to see if E or its quadratic twist meets the 
                # requirements
				for EE,tt in [(E,t), (E.quadratic_twist(), -t)]:
					if Zm(tt).multiplicative_order() == n  and ( 
                     all(EE.change_ring(k.extension(n//d)).cardinality()%m != 0
                     for d in n.prime_divisors()))  and (
                        EE.change_ring(K).cardinality()%m == 0):
						success = true
						return EE
				if not success:
					counter = counter + 1 
            else:
                # For now, we focus on the case m exactly prime.
                PZm = PolynomialRing(Zm, 'X')
                roots = PZm(E.frobenius_polynomial()).roots()
        
                # If it splits in GF(l) then it splits in GF(l^r)
                if len(roots) == 0:
                    counter = counter + 1
                    continue

                alpha = roots[0][0]
                beta = roots[1][0]

                # TODO : pensez a gerer ce cas, peut-etre
                if(Zm(alpha).multiplicative_order() ==
                        Zm(beta).multiplicative_order()):
                    continue

                # Root of the characteristic polynomials for the quadratic 
                # twist of E.
                # If r & t are roots of X + aX + b, then -r & -t are roots of
                # X - aX + b
                alpha_t = -alpha
                beta_t = -beta


                # We want the root of smallest order in (Z/m)*.
                if (Zm(alpha).multiplicative_order() 
                                > Zm(beta).multiplicative_order()):
                    alpha = beta

                # We could probably deduce that from the above condition.
                if (Zm(alpha_t).multiplicative_order()
                                > Zm(beta).multiplicative_order()):
                    alpha_t = beta_t

                success = false 

                for EE, a, b in [(E, alpha, beta), (E.quadratic_twist(),
                    alpha_t, beta_t)]:
                    if Zm(a).multiplicative_order() == n and (
                     all(EE.change_ring(k.extension(n//d,
                         conway= true, prefix = 'z')).cardinality()%m != 0
                     for d in n.prime_divisors())) and (
                           EE.change_ring(K).cardinality()%m == 0) and (
                        Zm(a).multiplicative_order != 
                        Zm(b).multiplicative_order):
                           success = true
                           return EE
                if not success:
                    counter = counter + 1

        return RuntimeError, 'No appropriate elliptic curves found.'

def accept_elliptic(r, e, n):
    '''
    This function accepts only if :

    (1) the order of t in Z/r^e is a multiple of n;
    (2) phi(m)/n is prime to n.

    or

    (1) for a root alpha of X^2 - tX + q, alpha or alpha/t has for order 
    a multiple of n in Z/m, the polynome must be separable.
    (2) phi(m)/n is prime to n.

    '''
    if r!=2 and t != p:
        m = r**e
        ord = (r-1) * r**(e-1) # Euler/Carmichael function.
        Zm = Zmod(m)
        if m%p == 0:

            return ( (ord // n.expand()).gcd(n.expand()) == 1 and
                    all(Zm(t)**(ord // ell) != 1 for (ell, _) in n))

        else:
            PZm = PolynomialRing(Zm, 'X'); X = PZm.gen()
            f = X**2 - Zm(t)*X + Zm(q)
            roots = f.roots(multiplicities=false)

            # If we have a root of multiplicity two, we don't want this trace
            if len(roots) <= 1:
                return false

            # We are interested in the roots of smallest order in (Z/m)*
            if (Zm(roots[1]).multiplicative_order() >
                    Zm(roots[0]).multiplicative_order()):
                a = roots[0]
            else:
                a = roots[1]

            return( (ord // n.expand()).gcd(n.expand()) == 1 and
                    all(Zm(a)**(ord//ell) != 1 for (ell, _) in n))
    elif r == 2:
        raise NotImplementedError, 'm = 2 is not implemented yet'
    else:
        return false





















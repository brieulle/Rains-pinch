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


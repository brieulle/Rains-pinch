from sage.all import *
from ellipticrains import *
from sage.ffisom.ellipticrains import isom_elliptic as isom_elliptic_flint
from sage.rings.finite_rings.hom_finite_field import FiniteFieldHomomorphism_generic as FFH
from sage.misc.misc import cputime
from rains import *
from Rains_isom_normal import *
import sys, getopt

# Evoulution of certain parameters for a single characteristic when the degree 
# is varying.
#1, 2
def test_ell_param(fixed_var, bound, start_var = 7, leap = 1, fixed_degree = False):
    if not fixed_degree:
        p = fixed_var
        n = start_var
        print '#n m t_E ord_m period'

        for i in range(bound):
            R = PolynomialRing(GF(p), name = 'X')
            f = R.irreducible_element(n, algorithm='random')
            k1 = GF(p**n, name = 'x', modulus = f)
            e = k1.random_element()
            g = e.minpoly()
            while(g.degree() != n):
                e = k1.random_element()
                g = e.minpoly()
            k2 = GF(p**n, name = 'y', modulus = g)

        # Times to compute m, times to find an E, times to compute a point of 
        # order m, time to compute the periods, (value of m, won't be used in 
        # this routine).
            w_m, w_E, w_ordm, w_period, m, trashvar = isom_elliptic(k1, k2)
        
        # Skipping first line because of an anomaly that gives an execution time
        # higher than it should be for one of the parameters. The culprit is 
        # probably a module called the first time.
            if i:
                print '%s %s %s %s %s' % (n, w_m, w_E, w_ordm, w_period)
        
            for j in range(leap):
                n = next_prime(n)
    else:
        n = fixed_var
        p = start_var
        print '#p m t_E ord_m period'

        for i in range(bound):
            R = PolynomialRing(GF(p), name = 'X')
            f = R.irreducible_element(n, algorithm='random')
            k1 = GF(p**n, name = 'x', modulus = f)
            g = 1
            e = k1.random_element()
            g = e.minpoly()
            while(g.degree() != n):
                e = k1.random_element()
                g = e.minpoly()
            k2 = GF(p**n, name = 'y', modulus = g)

            w_m, w_E, w_ordm, w_period, m, trashvar = isom_elliptic(k1, k2)
        
            if i:
                print '%s %s %s %s %s' % (p, w_m, w_E, w_ordm, w_period)
        
            for j in range(leap):
                p = next_prime(p)

#3, 4
def test_E_nbtrace(fixed_var, bound, start_var = 3, leap = 1, fixed_degree = False):
# Evolution of other parameters with p fixed.
    if not fixed_degree:
        p = fixed_var
        n = start_var
        print '#n m len(S_t) t_E compteur'

        for i in range(borne_nbn):
            m, S_t = find_m(n, GF(p))
            k1 = GF(p**n, name = 'x', modulus =
                    PolynomialRing(GF(p), name='X').irreducible_element(n,
                    algorithm='random'))

            w = cputime()
            n_E = find_elliptic_curve(GF(p), k1, (m, S_t))[-1]
            w_t = cputime(w)
            if i :
                #The degree, the value of m, the number of trace candidates,
                # time to find an E, number of elliptic curves tested.
                print '%s %s %s %s %s' % (n, m, len(S_t), w_t, n_E)

            for j in range(leap):
                n = next_prime(n)
# Evolution of other parameters with n fixed.
    else:
        n = fixed_var
        p = start_var
        print '#p m len(S_t) t_E compteur'

        for i in range(borne_nbn):
            m, S_t = find_m(n, GF(p))
            k1 = GF(p**n, name = 'x', modulus =
                    PolynomialRing(GF(p), name='X').irreducible_element(n,
                    algorithm='random'))

            w = cputime()
            n_E = find_elliptic_curve(GF(p), k1, (m, S_t))[-1]
            w_t = cputime(w)
            if i :
                #The degree, the value of m, the number of trace candidates,
                # time to find an E, number of elliptic curves tested.
                print '%s %s %s %s %s' % (p, m, len(S_t), w_t, n_E)

            for j in range(leap):
                p = next_prime(p)

def find_ext(p, borne_nbn, start_var, with_exp, leap):
    '''
    Function designed to return a list of n depending on p for which the 
    cyclotomic method of Rains requires the use of another extension or not,
    depending on the value of with_exp.
    '''
    sol = []
    if not with_exp:
        n = start_var
        for j in range(borne_nbn):
            o = find_root_order(p,n)[0]
            if o == 1:
                sol.append((n,-1))
                for j in range(leap):
                    n = next_prime(n)
            else:
                for j in range(leap):
                    n = next_prime(n)

        return sol
    else:
        n = start_var
        for j in range(borne_nbn):
            o = find_root_order(p,n)[0]
            if o > 1:
                sol.append((n,o))
                for j in range(leap):
                    n = next_prime(n)
            else:
                for j in range(leap):
                    n = next_prime(n)

        return sol

# 5 : with_exp = flint = False, 6 : with_exp = True, flint = False, 7 : with_exp
# = False, flint = True, 8 : with_exp = flint = True
def cmp_EllCycl(p, borne_nbn, with_exp = False, start_var = 7, leap = 1, flint =
False):

    print '#n t_cylc t_ell o' 
    liste = find_ext(p, borne_nbn, ZZ(start_var), with_exp, leap)
    i = 0   # Failsafe, prevent the first loop to be reported 'cause it ain't
            # workin' correctly.

    for n in liste:
        R = PolynomialRing(GF(p), name = 'X')
        f = R.irreducible_element(n[0], algorithm='random')
        g = R.irreducible_element(n[0], algorithm='random')
        if g == f:
            g = R.irreducible_element(n[0], algorithm='random') 
        

        if flint is False:
            k1 = GF(p**n[0], name = 'x', modulus = f)
            k2 = GF(p**n[0], name = 'y', modulus = g)

            w = cputime()
            find_gens_cyclotomic(k1, k2)
            t_cycl = cputime(w)
            w = cputime()
            isom_elliptic(k1, k2)
            t_ell = cputime(w) 


            if i:
                print '%s %s %s %s' % (n[0], t_cycl, t_ell, n[1]) 
            else:
                i += 1
        else:
            k1 = GF(p**n[0], name = 'x', modulus = f, impl = 'flint_fq_nmod')
            k2 = GF(p**n[0], name = 'y', modulus = g, impl = 'flint_fq_nmod')

            w = cputime()
            find_gens_cyclotomic(k1, k2)
            t_cycl = cputime(w)
            w = cputime()
            isom_elliptic_flint(k1, k2)
            t_ell = cputime(w) 


            if i:
                print '%s %s %s %s' % (n[0], t_cycl, t_ell, n[1]) 
            else:
                i += 1

#TODO: Si jamais je veux comparer l'implem flin et l'implem normale
def cmp_flint(fixed_var, bound, start_var = 7, leap = 1, fixed_degree = False): 
    if not fixed_degree:
        print '#n t_ell t_ell_flint m nb_t'    
        p = fixed_var
        n = start_var

        for i in range(bound):
            R = PolynomialRing(GF(p), name = 'X')
            f = R.irreducible_element(n, algorithm='random')
            g = R.irreducible_element(n, algorithm='random')
            if g == f:
                g = R.irreducible_element(n, algorithm='random') 

            k1 = GF(p**n, name = 'x', modulus = f)
            k2 = GF(p**n, name = 'y', modulus = g)
            l1 = GF(p**n, name = 'a', modulus = f, impl = 'flint_fq_nmod')
            l2 = GF(p**n, name = 'b', modulus = g, impl = 'flint_fq_nmod')

            w = cputime()
            isom_elliptic_flint(k1, k2)
            t_ell = cputime(w)
            w = cputime()
            isom_elliptic_flint(l1, l2)
            t_ell_flint = cputime(w) 

            m_t = find_m(n, GF(p))


            if i:
                print '%s %s %s %s %s' % (n, t_ell, t_ell_flint, m_t[0], 
                len(m_t[1])) 

            for j in range(leap):
                n = next_prime(n)
    else:
        print '#p t_ell t_ell_flint m nb_t'    
        n = fixed_var
        p = start_var

        for i in range(bound):
            R = PolynomialRing(GF(p), name = 'X')
            f = R.irreducible_element(n, algorithm='random')
            g = R.irreducible_element(n, algorithm='random')
            if g == f:
                g = R.irreducible_element(n, algorithm='random') 

            k1 = GF(p**n, name = 'x', modulus = f)
            k2 = GF(p**n, name = 'y', modulus = g)
            l1 = GF(p**n, name = 'a', modulus = f, impl = 'flint_fq_nmod')
            l2 = GF(p**n, name = 'b', modulus = g, impl = 'flint_fq_nmod')

            w = cputime()
            isom_elliptic_flint(k1, k2)
            t_ell = cputime(w)
            w = cputime()
            isom_elliptic_flint(l1, l2)
            t_ell_flint = cputime(w) 

            m_t = find_m(n, GF(p))


            if i:
                print '%s %s %s %s %s' % (p, t_ell, t_ell_flint, m_t[0], 
                len(m_t[1])) 

            for j in range(leap):
                p = next_prime(p)


#11, 12, 13, 14
def cmp_ellFFH(fixed_var, bound, start_var = 7, leap = 1, fixed_degree = False,
    allombert = False, flint = False):
    if not fixed_degree:
        print '#n t_ell t_all/naive m nb_t'    
        p = fixed_var
        n = start_var

        for i in range(bound):
            R = PolynomialRing(GF(p), name = 'X')
            f = R.irreducible_element(n, algorithm='random')
            g = R.irreducible_element(n, algorithm='random')
            if g == f:
                g = R.irreducible_element(n, algorithm='random') 

            if not flint:
                k1 = GF(p**n, name = 'x', modulus = f)
                k2 = GF(p**n, name = 'y', modulus = g)

                if not allombert:
                    w = cputime()
                    isom_elliptic_flint(k1, k2)
                    t_ell = cputime(w)
                    w = cputime()
                    FFH(Hom(k1, k2))
                    t_all = cputime(w) 
                else:
                    w = cputime()
                    isom_elliptic_flint(k1, k2)
                    t_ell = cputime(w)
                    w = cputime()
                    FFH(Hom(k1, k2), algorithm = "allombert")
                    t_all = cputime(w) 


                if i:
                    print '%s %s %s' % (n, t_ell, t_all)

                for j in range(leap):
                    n = next_prime(n)
            else:
                k1 = GF(p**n, name = 'a', modulus = f, impl = 'flint_fq_nmod')
                k2 = GF(p**n, name = 'b', modulus = g, impl = 'flint_fq_nmod')

                if not allombert:
                    w = cputime()
                    isom_elliptic_flint(k1, k2)
                    t_ell = cputime(w)
                    w = cputime()
                    FFH(Hom(k1, k2))
                    t_all = cputime(w) 
                else:
                    w = cputime()
                    isom_elliptic_flint(k1, k2)
                    t_ell = cputime(w)
                    w = cputime()
                    FFH(Hom(k1, k2), algorithm = "allombert")
                    t_all = cputime(w) 
        #            w = cputime()
        #            find_gens_cyclotomic(k1,k2)
        #            t_cycl = cputime(w)


                if i:
                    print '%s %s %s' % (n, t_ell, t_all)

                for j in range(leap):
                    n = next_prime(n)
    else:
        print 'TODO'

# 15, 16
def test_cycl(fixed_var, bound, start_var = 7, leap = 1, fixed_degree = False):
    if not fixed_degree:
        p = fixed_var
        n = start_var
        print '#n t_cycl o'

        for i in range(bound):
            R = PolynomialRing(GF(p), name='X')
            f = R.irreducible_element(n, algorithm='random')
            g = R.irreducible_element(n, algorithm='random')
            if f == g:
                g = R.irreducible_element(n, algorithm='random')

            k1 = GF(p**n, name = 'x', modulus = f)
            k2 = GF(p**n, name = 'y', modulus = g)

            w = cputime()
            find_gens_cyclotomic(k1, k2)
            t_cycl = cputime(w)

            o = find_root_order(p, n)[0]

            if i:
                print '%s %s %s' % (n, t_cycl, o)

            for j in range(leap):
                n = next_prime(n)
    else:
        n = fixed_var
        p = start_var
        print '#p t_cycl o'

        for i in range(bound):
            R = PolynomialRing(GF(p), name='X')
            f = R.irreducible_element(n, algorithm='random')
            g = R.irreducible_element(n, algorithm='random')
            if f == g:
                g = R.irreducible_element(n, algorithm='random')

            k1 = GF(p**n, name = 'x', modulus = f)
            k2 = GF(p**n, name = 'y', modulus = g)

            w = cputime()
            find_gens_cyclotomic(k1, k2)
            t_cycl = cputime(w)

            o = find_root_order(p, n)[0]

            if i:
                print '%s %s %s' % (p, t_cycl, o)

            for j in range(leap):
                p = next_prime(p)
#17, 18
def test_ell(fixed_var, bound, start_var = 7, leap = 1, fixed_degree = False):
    if not fixed_degree:
        p = fixed_var
        n = start_var
        print '#n t_ell m'

        for i in range(bound):
            R = PolynomialRing(GF(p), name='X')
            f = R.irreducible_element(n, algorithm='random')
            g = R.irreducible_element(n, algorithm='random')
            if f == g:
                g = R.irreducible_element(n, algorithm='random')

            k1 = GF(p**n, name = 'x', modulus = f)
            k2 = GF(p**n, name = 'y', modulus = g)

            w = cputime()
            isom_elliptic(k1, k2)
            t_ell = cputime(w)

            m = find_m(n, GF(p))[0]

            if i:
                print '%s %s %s' % (n, t_ell, m)

            for j in range(leap):
                n = next_prime(n)
    else:
        n = fixed_var
        p = start_var
        print '#p t_ell m'

        for i in range(bound):
            R = PolynomialRing(GF(p), name='X')
            f = R.irreducible_element(n, algorithm='random')
            g = R.irreducible_element(n, algorithm='random')
            if f == g:
                g = R.irreducible_element(n, algorithm='random')

            k1 = GF(p**n, name = 'x', modulus = f)
            k2 = GF(p**n, name = 'y', modulus = g)

            w = cputime()
            isom_elliptic(k1, k2)
            t_ell = cputime(w)

            m = find_m(n, GF(p))[0]

            if i:
                print '%s %s %s' % (p, t_ell, m)

            for j in range(leap):
                p = next_prime(p)

def test_normal(fixed_var, bound, start_var = 7, leap = 1, fixed_degree = False):
    if not fixed_degree:
        p = fixed_var
        n = start_var
        print '%n t_coef t_isom'

        for i in range(bound):
            R = PolynomialRing(GF(p), name='X')
            f = R.irreducible_element(n, algorithm='random')
            g = R.irreducible_element(n, algorithm='random')
            if f == g:
                g = R.irreducible_element(n, algorithm='random')

            k1 = GF(p**n, name = 'x', modulus = f)
            k2 = GF(p**n, name = 'y', modulus = g)

            o = find_root_order(p, n)[0]

            if o == 1:
                tuple = find_gens_cyclotomic(k1, k2)
            else:
                tuple = isom_elliptic_flint(k1, k2)

            t_coef, t_isom =  isom_normal(tuple[0], tuple[1], k1, k2)

            if i:
                print '%s %s %s' % (n, t_coef, t_isom)

            for i in range(leap):
                n = next_prime(n)
    else:
        print 'TODO'
#20, 21, 22
def cmp_ellcyclFFH_nocase(fixed_var, bound, start_var = 7, leap = 1, 
        naive = False, allombert = False):
    if not FFH:
        p = fixed_var
        n = start_var
        print '#n t_cycl o'

        for i in range(bound):
            R = PolynomialRing(GF(p), name='X')
            f = R.irreducible_element(n, algorithm='random')
            g = R.irreducible_element(n, algorithm='random')
            if f == g:
                g = R.irreducible_element(n, algorithm='random')

            k1 = GF(p**n, name = 'x', modulus = f, impl="flint_fq_nmod")
            k2 = GF(p**n, name = 'y', modulus = g, impl="flint_fq_nmod")

            w = cputime()
            find_gens_cyclotomic(k1, k2)
            t_cycl = cputime(w)
            
            w = cputime()
            isom_elliptic_flint(k1, k2)
            t_ell = cputime(w)

            o = find_root_order(p, n)[0]

            if i:
                print '%s %s %s %s' % (n, t_cycl, t_ell, o)

            for j in range(leap):
                n = next_prime(n)
    else:
        p = fixed_var
        n = start_var
        print '#n t_cycl t_all/norm'
        

        for i in range(bound):
            R = PolynomialRing(GF(p), name='X')
            f = R.irreducible_element(n, algorithm='random')
            g = R.irreducible_element(n, algorithm='random')
            if f == g:
                g = R.irreducible_element(n, algorithm='random')

            k1 = GF(p**n, name = 'x', modulus = f)
            k2 = GF(p**n, name = 'y', modulus = g)

            if not allombert:
                w = cputime()
                find_gens_cyclotomic(k1, k1)
                t_cycl = cputime(w)

                w =  cputime()
                FFH(Hom(k1, k2))
                t_allnaive = cputime(w)
            else:
                w = cputime()
                find_gens_cyclotomic(k1, k1)
                t_cycl = cputime(w)

                w =  cputime()
                FFH(Hom(k1, k2), algorithm="allombert")
                t_allnaive = cputime(w)

            if i:
                print '%s %s %s' % (n, t_cycl, t_allnaive)

            for j in range(leap):
                n = next_prime(n)

if __name__ == '__main__':
    test = int(sys.argv[1])
    if test == 1:
        test_ell_param(ZZ(sys.argv[2]), int(sys.argv[3]), ZZ(sys.argv[4]))
    elif test == 2:
        test_ell_param(ZZ(sys.argv[2]), int(sys.argv[3]), ZZ(sys.argv[4]), 
        fixed_degree = True)
    elif test == 3:
        test_E_nbtrace(ZZ(sys.argv[2]), int(sys.argv[3]), ZZ(sys.argv[4]))
    elif test == 4:
        test_E_nbtrace(ZZ(sys.argv[2]), int(sys.argv[3]), ZZ(sys.argv[4]),
        fixed_degree = True)
    elif test == 5:
        cmp_EllCycl(ZZ(sys.argv[2]), int(sys.argv[3]), start_var =
        ZZ(sys.argv[4]))
    elif test == 6:
        cmp_EllCycl(ZZ(sys.argv[2]), int(sys.argv[3]), with_exp = True, 
        start_var = ZZ(sys.argv[4]))
    elif test == 7:
        cmp_EllCycl(ZZ(sys.argv[2]), int(sys.argv[3]), start_var =
        ZZ(sys.argv[4]), flint = True)
    elif test == 8:
        cmp_EllCycl(ZZ(sys.argv[2]), int(sys.argv[3]), with_exp = True, 
        start_var = ZZ(sys.argv[4]), flint = True)
    elif test == 9:
        cmp_flint(ZZ(sys.argv[2]), int(sys.argv[3]), 
        start_var = ZZ(sys.argv[4]))
    elif test == 10:
        cmp_flint(ZZ(sys.argv[2]), int(sys.argv[3]), 
        start_var = ZZ(sys.argv[4]), fixed_degree = True)
    elif test == 11:
        cmp_ellFFH(ZZ(sys.argv[2]), int(sys.argv[3]), start_var =
        ZZ(sys.argv[4]))
    elif test == 12:
        cmp_ellFFH(ZZ(sys.argv[2]), int(sys.argv[3]), start_var =
        ZZ(sys.argv[4]), flint = True)
    elif test == 13:
        cmp_ellFFH(ZZ(sys.argv[2]), int(sys.argv[3]), start_var =
        ZZ(sys.argv[4]), allombert = True)
    elif test == 14:
        cmp_ellFFH(ZZ(sys.argv[2]), int(sys.argv[3]), start_var =
        ZZ(sys.argv[4]), allombert = True, flint = True)
    elif test == 15:
        test_cycl(ZZ(sys.argv[2]), int(sys.argv[3]), start_var = 
                ZZ(sys.argv[4]))
    elif test == 16:
        test_cycl(ZZ(sys.argv[2]), int(sys.argv[3]), start_var = 
                ZZ(sys.argv[4]), 
        fixed_degree=True)
    elif test == 17:
        test_ell(ZZ(sys.argv[2]), int(sys.argv[3]), start_var = ZZ(sys.argv[4]))
    elif test == 18:
        test_ell(ZZ(sys.argv[2]), int(sys.argv[3]), start_var = ZZ(sys.argv[4]),
        fixed_degree=True)
    elif test == 19:
        test_normal(ZZ(sys.argv[2]), int(sys.argv[3]), start_var = 
                ZZ(sys.argv[4]))
    elif test == 20:
        cmp_ellcyclFFH_nocase(ZZ(sys.argv[2]), int(sys.argv[3]), start_var = 
                ZZ(sys.argv[4]))
    elif test == 21:
        cmp_ellcyclFFH_nocase(ZZ(sys.argv[2]), int(sys.argv[3]), start_var = 
                ZZ(sys.argv[4]), naive = True)
    elif test == 22:
        cmp_ellcyclnaive_nocase(ZZ(sys.argv[2]), int(sys.argv[3]), start_var = 
                ZZ(sys.argv[4]), naive = True, allombert = True)
    else:
        raise RuntimeError, 'There\'s only 16 tests for now.'

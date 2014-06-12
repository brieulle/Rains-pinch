# -*- coding: utf-8 -*-
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
    PZm = PolynomialRing(Zm, 'X')
    X = PZm.gen()
    p = k.characteristic()
    q = k.cardinality()


    # If m is a multiple of p, then we just need the trace to be of order 
    #exactly n in (Z/m)*
    if m%p == 0:
        sol = []
        for t in Zm
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
        #TODO : rÃ©flÃ©chir Ã  la forme que a pourrait Ã©ventuellement avoir..
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
                elif (ord_b > ord_a):
                    continue
                else:
                    sol.append((q/a, a + q/a)) #return (q/a, a + q/a)
        return sol
        
    else:
        raise NotImplementedError, 'm composite is not implemented yet'

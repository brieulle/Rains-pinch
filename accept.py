# -*- coding: utf-8 -*-
def find_trace(n,m,k):
    '''
    We look for t's in ZZ such that ord_m(t) = n if  m is a power of p or a's in
    CC such that ord_m(a) = n or ord_m(q/a) = n whether ord_m(a) < ord_m(q/a) or
    ord_m(q/a) < ord_m(a) and a + q/a is in ZZ if m is just a power of a prime.
    The element a will be one of the roots of the minimal polynomail of the
    frobenius : 

    X - (a + a/q)X + a.(q/a)

    So actually, we also need the above polynomial to be separable in Z/m.
    '''
    Zm = Zmod(m)
    PZm = PolynomialRing(Zm, 'X')
    X = PZm.gen()
    p = k.characteristic()
    q = k.cardinality()


    # If m is a multiple of p, then we just need the trace to be of order 
    #exactly n in (Z/m)*
    if m%p == 0:
        count = 0
        t = ZZ(1)
        while true:
            if count >= euler_phi(m):
                raise RuntimeError, 'No suitable trace found.'
            # We only want the trace to be of order exactly n in (Z/m)* and not 
            # to define supersingular curves.
            if t.gcd(m) != 1:
                count +=1
                continue
            elif (Zm(t).multiplicative_order() != n):
                t += 1
                count +=1
            else:
                return t 
    # If m is prime (power), then we need to look at the roots
    # Probably a temporary condition
    elif m.is_prime_power():
        #TODO : rÃ©flÃ©chir Ã  la forme que a pourrait Ã©ventuellement avoir..
        count = 0
        sol = []
        while true:
            for a in Zm:
                if (count >= 500):
                    return list(set(sol)) #raise RuntimeError, 'Did not found any workable trace'
                # If a is not invertible in Z/m, we can't compute any order.
                # We'll probably have to look for an element of order exactlu n
                # We don't want a = 0 mod m. 
                if not a.is_unit():
                    count += 1
                    continue
                else:
                    ord_a = Zm(a).multiplicative_order()
                    ord_b = Zm(q/a).multiplicative_order()
                    # We need an element of order n
                    if (ord_a != n and ord_b !=n):
                        count +=1
                        continue
                    # We want ord_a != ord_b
                    elif (ord_a == ord_b):
                        count +=1
                        continue
                    elif (ord_b != n): 
                        if (ord_a > ord_b):
                            count += 1
                            continue
                        else:
                            count += 1
                            sol.append(a) #return (a, a + q/a, q, q/a)
                    elif (ord_b > ord_a):
                        count += 1
                        continue
                    else:
                        count += 1
                        sol.append(q/a) #return (q/a, a + q/a)
                # If we pass through all invertible and none corresponds, then 
                # the m is probably incorrect.

    else:
        raise NotImplementedError, 'm composite is not implemented yet'

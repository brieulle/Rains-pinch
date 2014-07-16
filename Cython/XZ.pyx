cpdef getBit(object n, int k):
    return (n & (1 << k)) >> k

cpdef point(object P):
    if P[-1] == 0:
        return (0,0)
    else:
        return (P[0]/P[-1], 1)

# Using dbl-2002-bj
cpdef doubling(object P, object a, object b):
    Z1 = P[-1]
    X1 = P[0]

    T1 = X1*Z1**2
    X12 = X1**2

    X3 = (X1**2 - a*Z1**2)**2 - 8*b*T1*Z1
    Z3 = 4*Z1*(X12*X1 + a*T1 + b*Z1**3)

    return X3, Z3

# Using dadd-2002-it
cpdef dadd(object P, object Q, object diff, object a, object b):
    if Q[-1] == 0:
        return (P[0], P[-1])
    elif P[-1] == 0:
        return (Q[0], Q[-1])
    elif diff[-1] == 0:
        return doubling(P, a, b)
    else:
        Z1 = diff[-1]
        Z2 = P[-1]
        Z3 = Q[-1]
        X1 = diff[0]
        X2 = P[0]
        X3 = Q[0]

        T1 = X2*Z3
        S1 = X3*Z2

        X5 = Z1*((X2*X3 - a*Z2*Z3)**2 - 4*b*Z2*Z3*(T1 + S1))
        Z5 = X1*(T1 - S1)**2

        return X5, Z5

cpdef ladder(object P, object m, object a, object b, object E = None):
    if E is None:
        S = (0, 0)
    else:
        S = E(0)

    R = P

    for k from 0 <= k <= m.nbits():
        if getBit(m, m.nbits() - k) == 0:
            R = dadd(R, S, P, a, b)
            S = doubling(S, a, b)
        else:
            S = dadd(S, R, P, a, b)
            R = doubling(R, a, b)

    return point(S)

cpdef find_ordm(object E, object m):
    cofactor = E.cardinality()//m
    coprime = m.prime_divisors()
    size = len(coprime)

    P = E(0)

    while(1):
        count = 0
        for a in coprime:
            m_a = m//a
            if ladder(P, m_a, E.a4(), E.a6())[1] == 0:
                continue
            else:
                count += 1

        if count != size:
            P = ladder(E.random_point(), cofactor, E.a4(), E.a6())
        elif ladder(P, m, E.a4(), E.a6())[1] != 0:
            P = ladder(E.random_point(), cofactor, E.a4(), E.a6())
        else:
            return P

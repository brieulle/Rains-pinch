cdef point(object P):
    if P[1] == 0:
        return (0,0)
    else:
        return (P[0]/P[1], 1)

cdef dadd(object P, object Q, object diff, object A):
    if Q[1] == 0:
        return (P[0], P[1])
    elif P[1] == 0:
        return (Q[0], Q[1])
    elif diff[1] == 0:
        return doubling(P, A)
    else:
        Z1 = diff[1]
        Z2 = P[1]
        Z3 = Q[1]
        X1 = diff[0]
        X2 = P[0]
        X3 = Q[0]

        A = X2 + Z2
        B = X2 - Z2
        C = X3 + Z3
        D = X3 - Z3

        DA = D*A
        CB = C*B

        X5 = Z1*(DA + CB)**2
        Z5 = X1*(DA - CB)**2

        return (X5, Z5)


cdef doubling(object P, object A):
        A24 = (A + 2)/4
        X1 = P[0]
        Z1 = P[1]

        a = (X1 + Z1)**2
        b = (X1 - Z1)**2
        c = (a - b)
        
        X3 = a*b
        Z3 = c*(b + ((A+2)/4)*c)

        return (X3, Z3)


cdef ladder(object P, object m, object A, object E = None):
    if E is None:
        S = (0, 0)
    else:
        S = E(0)

    R = P
    bits = m.binary()

    for bit in bits:
        if bit == '0':
            R = dadd(R, S, P, A)
            S = doubling(S, A)
        else:
            S = dadd(S, R, P, A)
            R = doubling(R, A)

    return point(S)

cpdef find_ordm(object E, object m, object A):
    cofactor = E.cardinality()//m
    coprime = m.prime_divisors()
    size = len(coprime)

    while(1):
        count = 0
        for a in coprime:
            m_a = m//a
            if ladder((0,0), m_a, A)[1] == 0:
                continue
            else:
                count += 1

        if count != size:
            temp = E.random_point()
            P = ladder((temp[0], 1), cofactor, A)
        else:
            return P

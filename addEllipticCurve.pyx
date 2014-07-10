cdef class MontgomeryCurve(SageObject):
    def __init__(self, A, B):
        self.A = A
        self.B = B
        self.base_field = A._parent



cdef class MontgomeryPoint(SageObject):
    def __init__(self, MontgomeryCurve E, temptype x, temptype z, y = None):
        self.curve = E
        self.x = x
        self.z = z
        self.y = y 


    cpdef MontgomeryPoint dadd(self, MontgomeryPoint P, MontgomeryPoint diff):
        cdef temptype da, cb

        da = (self.x + self.z)*(P.x - P.z)
        cb = (self.x - self.z)*(P.x - P.z)

        return MontgomeryPoint(self.curve, diff.z*(da - cb).square(), diff.x*(da - cb).square())

    cpdef MontgomeryPoint doubling(self):
        cdef temptype a, b, c




cdef add(point P, int n, int m):
    cdef int* x = NULL
    cdef int* z = NULL

#cdef int* double(int x, int z, const int A):
#    return -1

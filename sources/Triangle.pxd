from Triangle cimport Triangle
from libcpp.vector cimport vector
cdef extern from "basic_types/triangle.h":
    cdef cppclass Triangle:
        Triangle() except +
        int TV(int pos)
        void TE(int pos, vector[int]& e)
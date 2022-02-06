from Mesh cimport Mesh
from libcpp.string cimport string
from libcpp cimport bool

# Declare the class Mesh with cdef
cdef extern from "io/reader.h":
    cdef cppclass Reader:
        @staticmethod
        bool read_mesh(Mesh& mesh, string path);        
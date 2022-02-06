# distutils: language = c++
from libcpp cimport bool
from libcpp.vector cimport vector
from cython.operator cimport dereference as deref


from Mesh cimport Mesh
cdef class PyMesh:
    cdef Mesh *c_mesh
    def __cinit__(self):
        self.c_mesh = new Mesh()

    def get_vertex(self, int pos):
        vertex = PyVertex() 
        vertex.c_vertex = self.c_mesh.get_vertex(pos)
        return vertex

    def get_triangle(self, int pos):
        triangle = PyTriangle() 
        triangle.c_triangle = self.c_mesh.get_triangle(pos)
        return triangle

    # def __dealloc__(self):  ## will cause segmentation fault. 
    #     del self.c_mesh

from Vertex cimport Vertex
cdef class PyVertex:
    cdef Vertex c_vertex

    def __cinit__(self):
        self.c_vertex = Vertex()
    
    def get_c(self, int pos):
        return self.c_vertex.get_c(pos)

from Triangle cimport Triangle
cdef class PyTriangle:
    cdef Triangle c_triangle

    def __cinit__(self):
        self.c_triangle = Triangle()
    
    def TV(self, int pos):
        return self.c_triangle.TV(pos)

    def TE(self, int pos):
        cdef vector[int] e
        self.c_triangle.TE(pos, e)
        return e

from Reader cimport Reader
cdef class PyReader:
    cdef Reader *c_reader  # declare as pointer to avoid using constructor, which is private 

    def Py_read_mesh(self, PyPRT_Tree tree, str path):
        path_new = bytes(path, encoding='utf8')
        self.c_reader.read_mesh(tree.c_pt_pr_tree.get_mesh(), path_new)

         
    
    # def Py_read_mesh_old(self, PyMesh mesh, str path):   # can work, but will cause Segmentation fault (core dumped). Probably because c_mesh is deleted after tree is deleted?
    #     path_new = bytes(path, encoding='utf8')
    #     self.c_reader.read_mesh(deref(mesh.c_mesh), path_new)


from Node_V cimport Node_V
cdef class PyNode_V:
    cdef Node_V *c_node_v 
    # cdef vector[set[int]] vvs

    def __cinit__(self):
        self.c_node_v = new Node_V()

    def get_v_start(self):
        return self.c_node_v.get_v_start()
    
    def get_v_end(self):
        return self.c_node_v.get_v_end()

    def get_VT(self, PyMesh mesh):
        cdef vector[vector[int]] vts
        self.c_node_v.get_VT(vts, deref(mesh.c_mesh))
        return vts

    def is_leaf(self):
        return self.c_node_v.is_leaf()

    def get_son(self,  int i):
        node_v = PyNode_V()
        if self.c_node_v.get_son(i) is NULL:
            print("Null")
            return None
        else:
            node_v.c_node_v = self.c_node_v.get_son(i)
            return node_v
            
    def indexes_vertices(self):
        return self.c_node_v.indexes_vertices()


    # def __dealloc__(self):
    #     del self.c_node_v

from PRT_Tree cimport PRT_Tree
cdef class PyPRT_Tree:
    cdef PRT_Tree *c_pt_pr_tree
    # cdef Mesh * c_mesh
    def __cinit__(self, int v_per_leaf, int division_type):
        self.c_pt_pr_tree = new PRT_Tree(v_per_leaf, division_type)

    def build_tree(self):
        self.c_pt_pr_tree.build_tree()

    def get_mesh(self, PyMesh mesh):
        mesh.c_mesh = &(self.c_pt_pr_tree.get_mesh()) 

    def get_root(self, PyNode_V node):
        node.c_node_v = &(self.c_pt_pr_tree.get_root()) 


    # def get_leaves_number(self):
    #     return self.c_pt_pr_tree.get_leaves_number()
    
    def __dealloc__(self):
        del self.c_pt_pr_tree    

from Reindexer cimport Reindexer
cdef class PyReindexer:
    cdef Reindexer c_reindexer
    cdef vector[int] original_vertex_indices
    cdef vector[int] original_triangle_indices

    def __cinit__(self):
        self.c_reindexer = Reindexer()

    def reindex_tree_and_mesh(self, PyPRT_Tree pr_tree, bool save_v_indices,bool save_t_indices):
        self.c_reindexer.reindex_tree_and_mesh(deref(pr_tree.c_pt_pr_tree),save_v_indices,self.original_vertex_indices, save_v_indices,self.original_triangle_indices)
    
    # def __dealloc__(self):
    #     del self.original_vertex_indices    
    #     del self.original_triangle_indices


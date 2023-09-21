#ifndef SEA_ICE_PROCESSOR_H
#define SEA_ICE_PROCESSOR_H
#include "morse/forman_gradient_aux_structure.h"
#include "basic_types/mesh.h"
#include <unordered_map>

class Sea_Ice_Processor{
    public:
        ~Sea_Ice_Processor(){};
        Sea_Ice_Processor(simplices_multimap& input){this->extracted_triangles = input;};
        simplices_multimap get_processed_triangles(coord_type mode, Mesh& mesh);
    private:    
        void add_triangle_to_simplex(const ivect& triangle, itype simplex);
        itype find_root(itype simplex);
        void connect(itype simplex1, itype simplex2);
        // Returns true if the triangle is lower/equal to the sea ice level (mode elevation). 
        // If at least one vertex of the triangle is above the sea ice level, returns false.
        bool is_lower_than_level(const ivect& triangle, coord_type mode, Mesh& mesh);
        void update_label();
        simplices_multimap extracted_triangles; 
        vector<itype> size_of_simplices;
        vector<itype> root_of_simplices;
        // unordered_map<itype, itype> size_of_simplices;
        // unordered_map<itype, itype> root_of_simplices;
        map<pair<itype, itype>, itype> root_of_edges;
};

#endif // SEA_ICE_PROCESSOR_H

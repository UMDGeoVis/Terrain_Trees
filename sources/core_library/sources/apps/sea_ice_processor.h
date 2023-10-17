#ifndef SEA_ICE_PROCESSOR_H
#define SEA_ICE_PROCESSOR_H
#include "morse/forman_gradient_aux_structure.h"
#include "basic_types/mesh.h"
#include <unordered_map>

class Sea_Ice_Processor{
    public:
        ~Sea_Ice_Processor(){};
        Sea_Ice_Processor(simplices_multimap& input, const double area_limit){this->extracted_triangles = input; this->area_limit = area_limit;};
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
        double area_limit;
        // unordered_map<itype, itype> size_of_simplices;
        // unordered_map<itype, itype> root_of_simplices;
        map<pair<itype, itype>, itype> root_of_edges;
        inline double compute_tri_area(const ivect& triangle, Mesh& mesh){
            coord_type area;
            
            Vertex &vi = mesh.get_vertex(triangle[0]);
            Vertex &vj = mesh.get_vertex(triangle[1]);
            Vertex &vk = mesh.get_vertex(triangle[2]);

            dvect ki={vi.get_x()-vk.get_x(),vi.get_y()-vk.get_y()};
            dvect ij={vj.get_x()-vi.get_x(),vj.get_y()-vi.get_y()};
                // compute the vectors of the edges after rotate by 90 degrees
            dvect ki_vert= {vk.get_y() - vi.get_y() , vi.get_x()-vk.get_x()};
            dvect ij_vert= {vi.get_y() - vj.get_y() , vj.get_x()-vi.get_x()};
            // compute the area of the triangle 
            area=abs(0.5*(ki[0]*ij[1]-ki[1]*ij[0]));
            return area;

        };
};

#endif // SEA_ICE_PROCESSOR_H

#ifndef SEA_ICE_PROCESSOR_H
#define SEA_ICE_PROCESSOR_H
#include "morse/forman_gradient_aux_structure.h"
#include "basic_types/mesh.h"
#include <unordered_map>

class Sea_Ice_Processor{
    public:
        ~Sea_Ice_Processor(){};
        Sea_Ice_Processor(simplices_multimap& input, const double length_limit, coord_type mode){
            for(auto tri:input){
                extracted_triangles.push_back(tri.first);
            }
            this->length_limit = length_limit;
            this->area_mode = mode;
            };
        simplices_multimap get_processed_triangles(Mesh& mesh);
        inline void enable_peak_elevation_filter(bool flag){this->peak_elevation_filter = flag;}
        vector<pair<Vertex, Vertex>> get_ridge_paths_edges(){return this->filtered_ridge_paths;}
        inline void enable_roughness_filter(bool flag, coord_type threshold) {this->roughness_filter = flag; this->roughness_limit = threshold;}
        vector<pair<Vertex, Vertex>> get_ridge_paths_edges_new(Mesh& mesh, map<itype, vector<ivect>>& valid_ridge_paths);
    
    private:    
        void add_triangle_to_simplex(const ivect& triangle, itype simplex);
        itype find_root(itype simplex);
        void connect(itype simplex1, itype simplex2, const ivect& tri1 ,const ivect& tri2, Mesh& mesh);
        // Returns true if the triangle is lower/equal to the sea ice level (mode elevation). 
        // If at least one vertex of the triangle is above the sea ice level, returns false.
        bool is_lower_than_level(const ivect& triangle, coord_type mode, Mesh& mesh);
        bool has_valid_roughness(const ivect& triangle, coord_type threshold, int roughness_fid, Mesh& mesh);
        simplices_multimap update_label();
        vector<ivect> extracted_triangles; 
        vector<coord_type> size_of_simplices;
        vector<itype> root_of_simplices;
        vector<coord_type> length_of_simplices;
        vector<coord_type> peak_elevation_of_simplices;
        vector<pair<Vertex, Vertex>> ridge_paths;
        vector<vector<itype>> edges_for_each_root;
        vector<pair<Vertex, Vertex>> filtered_ridge_paths;
        double length_limit;
        double roughness_limit;
        bool peak_elevation_filter = false; 
        bool roughness_filter = false;
        double area_mode; // level sea ice elevation;
        // unordered_map<itype, itype> size_of_simplices;
        // unordered_map<itype, itype> root_of_simplices;
        map<pair<itype, itype>, itype> root_of_edges;
        // Note that the triangles are saved as list of vids by 1ascending extraction function, 
        // so here both the area and the centroid cannot be directly calculated by geometry_wrapper functions
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
        inline Vertex get_centroid(const ivect& triangle, Mesh& mesh){
            Vertex &vi = mesh.get_vertex(triangle[0]);
            Vertex &vj = mesh.get_vertex(triangle[1]);
            Vertex &vk = mesh.get_vertex(triangle[2]);
            Vertex p;
            for(int i = 0; i < 3; i++)
                p.set_c(i, (vi.get_c(i) + vj.get_c(i) + vk.get_c(i)) / 3.0);
            return p;
        }
};

#endif // SEA_ICE_PROCESSOR_H

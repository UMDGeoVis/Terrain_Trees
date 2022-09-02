/*
    This file is part of the Terrain Trees library.

    Author(s): Riccardo Fellegara (riccardo.fellegara@gmail.com)

    This project has been supported by the Italian Ministry of Education and
    Research under the PRIN 2009 program, and by the National Science Foundation
    under grant number IIS-1116747.

    The Terrain Trees library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The Terrain Trees library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the Terrain Trees library.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _WRITER_H
#define	_WRITER_H

#include <string>
#include <set>

#include <fstream>
#include <queue>
#include <iostream>
#include <boost/function.hpp>

#include "terrain_trees/tree.h"
#include "statistics/index_statistics.h"
#include "statistics/full_query_statistics.h"
#include "basic_types/box.h"

#include "terrain_trees/node_v.h"
#include "terrain_trees/node_t.h"

using namespace std;

typedef struct{
    int v[4];
} box_4;

///A class that provides an interface for writing to file or standard output some data structures or statistics
class Writer {
public:
    ///A public method that writes to file the tree structure
    /*!
     * \param fileName a string argument, representing the output path where the tree structure will be written
     * \param root a N& argument, representing the root of the tree to visit
     * \param division a D& argument, representing the space subdivision type used for the spatial index
     */
    template<class N, class D> static void write_tree(string fileName, N& root, D& division);
    ///A public method that writes to standard output the spatial index statistics
    /*!
     * \param indexStats an IndexStatistics& argument, representing the statistics to save
     */
    static void write_tree_stats(IndexStatistics& indexStats);
    ///A public method that writes to standard output the query statistics
    /*!
     * \param size an integer argument, representing the number of query executed
     * \param fullQueryStats an FullQueryStatistics& argument, representing the statistics to save
     * \param hit_ratio an integer representing the number of successfully answered queries
     */
    static void write_queries_stats(int size, FullQueryStatistics& fullQueryStats, int hit_ratio);
    ///A public method that writes to file a series of points that will be used as point location input
    /*!
     * \param points a set<Point>& argument, representing the points list to save
     * \param fileName a string argument, representing the file name
     */
    static void write_point_queries(set<Point> &points, string fileName);
    ///A public method that writes to file a series of boxes that will be used as box query input
    /*!
     * \param boxes a set<Box>& argument, representing the boxes list to save
     * \param fileName a string argument, representing the file name
     */
    static void write_box_queries(set<Box>& boxes, string fileName);

    ///A public method that writes to file the tree structure (VTK format)
    /*!
     * \param fileName a string argument, representing the output path where the tree structure will be written
     * \param root a N& argument, representing the root of the tree to visit
     * \param division a D& argument, representing the space subdivision type used for the spatial index
     * \param mesh representing the triangle mesh
     */
    template<class N, class D> static void write_tree_VTK(string file_name, N &root, D &division, Mesh &mesh);

    /**
     * @brief A public method that writes to file a triangle mesh in OFF format
     *
     * @param mesh_name a string argument represents the output path of the mesh
     * @param operation_type a string containing the name of the application from which we have obtained the new mesh file
     * @param mesh representing the triangle mesh to save
     * @param extra_fields a boolean stating if there are extra fields than the elevation encoded by the triangle mesh
     */
    static void write_mesh(string mesh_name, string operation_type, Mesh &mesh, bool extra_fields); /// OFF format
    /**
     * @brief A public method that writes to file a triangle mesh in VTK format
     *
     * @param mesh_name a string argument represents the output path of the mesh
     * @param mesh representing the triangle mesh to save
     */
    static void write_mesh_VTK(string mesh_name, Mesh &mesh);
    /**
     * @brief A public method that writes to file a triangle mesh and the curvature values in VTK format
     *
     * @param mesh_name a string argument represents the output path of the mesh
     * @param mesh representing the triangle mesh to save
     * @param curvature_type a string containing the type of curvature used
     * @param c_pos an integer representing the curvature position in the vertex array
     */
    static void write_mesh_curvature_VTK(string mesh_name, Mesh &mesh, string curvature_type, int c_pos);

    static void write_filtered_points_cloud(string mesh_name, Mesh &mesh); /// SpatialHadoop format
    static void write_filtered_points_cloud_with_id(string mesh_name, Mesh &mesh); /// SpatialHadoop format with vertex index
    static void write_multifield_points_cloud(string mesh_name, vertex_multifield &multifield, Mesh &mesh);

    static void write_critical_points(string mesh_name, map<short, set<ivect> > &critical_simplices, Mesh &mesh);    

    static void write_edge_costs(string mesh_name, vector<pair<coord_type, int>>& edge_costs);


    template<class N, class D> static void write_edge_cost_distribution(N& root,D& division, string file_name, Mesh &mesh, string& csv_path);

    static void write_edge_costs_vtk(string mesh_name, Mesh &mesh, map<pair<int, int>, coord_type>& edge_costs);

    static void write_boundary_matrix_pair(string mesh_name, Mesh &mesh, map<int,int> &pairs, vector<vector<int>>& index_to_simplex);

protected:
    ///A constructor method
    Writer() {}
    ///A constructor method
    Writer(const Writer&) {}
    ///A destructor method
    virtual ~Writer() {}

    ///A private method that writes to an output stream a tree node information (Generic Version)
    /*!
     * \param output an ofstream& argument, representing the stream
     * \param n a Node_T* argument, representing the node to save
     */
    static void write_node(ofstream& output, Node_T* n);
    ///A private method that writes to an output stream a tree node information (P_Node Version)
    /*!
     * \param output an ofstream& argument, represents the stream
     * \param n a Node_V* argument, represents the node to save
     */
    static void write_node(ofstream& output, Node_V *n);

    template<class N, class D> static void extract_leaf_domains(N& n, D& division, Box& domain, int level, vector<Point>& all_points,
                                                                vector<box_4>& all_leaves, Mesh& mesh);
};

template<class N, class D> void Writer::write_tree(string fileName, N& root, D& division)
{
    ofstream output(fileName.c_str());
    output.unsetf( std::ios::floatfield ); // floatfield not set
    output.precision(15);
    queue<N*> coda;
    N* visited;

    coda.push(&root);
    bool is_root = true;

    while (!coda.empty()) {
        visited = coda.front();
        coda.pop();

        string begin;

        if (visited->is_leaf())
            begin = "L";
        else
            begin = "N";

        if(!is_root)
            output << endl;
        output << begin << " ";

        is_root = false;

        Writer::write_node(output, visited);

        if (visited->is_leaf() == false)
            for (int i = 0; i < division.son_number(); i++)
                coda.push(visited->get_son(i));
    }
    output.close();
}

template<class N, class D> void Writer::write_tree_VTK(string file_name, N &root, D &division, Mesh &mesh)
{
    ofstream output(file_name.c_str());
    output.unsetf( std::ios::floatfield ); // floatfield not set
    output.precision(15);

    vector<box_4> all_leaves;
    vector<Point> all_points; //with_duplicate

    output<<"# vtk DataFile Version 2.0" << endl << endl
         << "ASCII" << endl << "DATASET UNSTRUCTURED_GRID " <<  endl << endl;

    Writer::extract_leaf_domains(root,division,mesh.get_domain(), 0, all_points, all_leaves, mesh);

    output<< "POINTS " << all_points.size() << " float" << endl;

    for(unsigned int v=0; v<all_points.size(); v++)
        output<<all_points.at(v).get_x()<<" "<<all_points.at(v).get_y()<<" -20"<<endl; /// fixed constant

    output<<endl << "CELLS " << all_leaves.size() << " " << (all_leaves.size()*5) << endl;

    for(unsigned int b=0; b<all_leaves.size(); b++)
        output<<"4 "<<all_leaves.at(b).v[0]<<" "<<all_leaves.at(b).v[1]<<" "<<all_leaves.at(b).v[2]<<" "<<all_leaves.at(b).v[3]<<endl;

    output<< endl << "CELL_TYPES " << all_leaves.size() << endl;
    for (unsigned int i = 0; i < all_leaves.size(); ++i)
        output<< "9 ";
    output<< endl;

    output<< "POINT_DATA " << all_points.size() << endl << endl;
    output<< "FIELD FieldData 1" << endl << endl;
    output<< "fieldvalue 1 " << all_points.size() << " float" << endl;

    for (unsigned int i =0; i <all_points.size(); ++i)
        output<<i<< " ";
    output<<endl;

    output.close();
}

template<class N, class D> void Writer::extract_leaf_domains(N &n, D &division, Box &domain, int level, vector<Point> &all_points,
                                                             vector<box_4>& all_leaves, Mesh& mesh)
{
    if(n.is_leaf())
    {
        //implementazione non bella a vedersi ma che fa quello che deve fare..
        //da migliorare se si trova il tempo..
        Point min = domain.get_min();
        Point max = domain.get_max();

        vector<Point> tmp;

        tmp.push_back(min);
        tmp.push_back(Point(min.get_x(),max.get_y()));        
        tmp.push_back(max);
        tmp.push_back(Point(max.get_x(),min.get_y()));

        box_4 newOne;

        for(unsigned int i=0; i<tmp.size(); i++)
        {
            all_points.push_back(tmp.at(i));
            newOne.v[i] = all_points.size() - 1;
        }

        all_leaves.push_back(newOne);
    }
    else
    {
        int son_level = level + 1;

        for (int i = 0; i < division.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                Box son_dom = division.compute_domain(domain,level,i);
                Writer::extract_leaf_domains(*n.get_son(i),division,son_dom,son_level,all_points,all_leaves,mesh);
            }
        }
    }
}

template<class N, class D> void Writer::write_edge_cost_distribution(N &root, D &division, string file_name, Mesh &mesh, string& csv_path){


    ofstream output(file_name.c_str());
    output.unsetf( std::ios::floatfield ); // floatfield not set
    output.precision(15);
    vector<box_4> all_leaves;
    vector<Point> all_points; //with_duplicate

    ifstream csvFile(csv_path);
    string id_str;
    if(!csvFile.is_open()) throw std::runtime_error("Could not open file");
    getline(csvFile, id_str);
    map<int, int> edgeCount;
    while(getline(csvFile, id_str, ',')){
        //stringstream ss(line);
        int id, num;
        id = atoi(id_str.c_str());
        
        string count_str;
        getline(csvFile, count_str);
        num = atoi(count_str.c_str());
        
        //ss>>id>>tmp>>num;
        edgeCount[id] = num;
        // cout<<"id: "<<id<<" count: "<<num<<endl;
    }

    output<<"# vtk DataFile Version 2.0" << endl << endl
         << "ASCII" << endl << "DATASET UNSTRUCTURED_GRID " <<  endl << endl;

    Writer::extract_leaf_domains(root,division,mesh.get_domain(), 0, all_points, all_leaves, mesh);
    cout<< "number of leaves: "<<all_leaves.size()<<endl;
    output<< "POINTS " << all_points.size() << " float" << endl;

    for(unsigned int v=0; v<all_points.size(); v++)
        output<<all_points.at(v).get_x()<<" "<<all_points.at(v).get_y()<<" -20"<<endl; /// fixed constant

    output<<endl << "CELLS " << all_leaves.size() << " " << (all_leaves.size()*5) << endl;

    for(unsigned int b=0; b<all_leaves.size(); b++)
        output<<"4 "<<all_leaves.at(b).v[0]<<" "<<all_leaves.at(b).v[1]<<" "<<all_leaves.at(b).v[2]<<" "<<all_leaves.at(b).v[3]<<endl;

    output<< endl << "CELL_TYPES " << all_leaves.size() << endl;
    for (unsigned int i = 0; i < all_leaves.size(); ++i)
        output<< "9 ";
    output<< endl;

    output<< "CELL_DATA "<<all_leaves.size()<<endl;
    output<< "scalars count int"<<endl;
    output<< "LOOKUP_TABLE default"<<endl;
    //int leaf_id = 0;
    for(int b=0; b<all_leaves.size(); b++){
        float val=0;
        if( edgeCount.find(b) != edgeCount.end() )
        { 

            val = edgeCount[b];

        }
        output<<val<<" ";
    }

    output<< endl;
    output<< "scalars avg_cost float"<<endl;
    output<< "LOOKUP_TABLE default"<<endl;
    //int leaf_id = 0;
    for(int b=0; b<all_leaves.size(); b++){
        float val=0;
        if( edgeCount.find(b) != edgeCount.end() )
        { 
            Point corner_left = all_points.at(all_leaves.at(b).v[0]);
            Point corner_right = all_points.at(all_leaves.at(b).v[2]);
            
            double x_dif = abs(corner_left.get_x() - corner_right.get_x());
            double y_dif = abs(corner_left.get_y() - corner_right.get_y());
          //  cout<<"x dif: "<<x_dif<<", y dif: "<<y_dif<<endl;
            val = double(edgeCount[b])/x_dif / y_dif;

        }
        output<<val<<" ";
    }

    output.close();
}


#endif	/* _WRITER_H */


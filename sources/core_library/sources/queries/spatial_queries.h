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

#ifndef SPATIAL_QUERIES_H
#define SPATIAL_QUERIES_H

#include "statistics/statistics.h"
#include "utilities/timer.h"
#include "priority_queue_element.h"
#include "terrain_features/critical_points_extractor.h"

/**
 * @brief The Spatial_Queries class provides an interface for executing spatial queries on the Terrain Trees
 */
class Spatial_Queries
{
public:
    Spatial_Queries() {}

    ///A public method that excutes point locations, reading the points from file
    /*!
     * This method prints the results on standard output
     *
     * \param tree a T& argument, represents the tree where the statistics are executed
     * \param query_path a string argument, representing the file path of the query input
     * \param stats a Statistics& argument, representing the object for computing the associated statistics
     */
    template<class T> void exec_point_locations(T& tree, string query_path, Statistics &stats);
    ///A public method that excutes point locations
    /*!
     * This method prints the results on standard output
     *
     * \param tree a T& argument, represents the tree where the statistics are executed
     * \param points a set containing the points for the queries
     * \param stats a Statistics& argument, representing the object for computing the associated statistics
     */
    template<class T> void exec_point_locations(T& tree, set<Point> &points, Statistics &stats);
    ///A public method that excutes box queries, reading the boxes from file
    /*!
     * This method prints the results on standard output
     *
     * \param tree a T& argument, represents the tree where the statistics are executed
     * \param query_path a string argument, representing the file path of the query input
     * \param stats a Statistics& argument, representing the object for computing the associated statistics
     */
    template<class T> void exec_box_queries(T& tree, string query_path, Statistics &stats);
    ///A public method that excutes box queries
    /*!
     * This method prints the results on standard output
     *
     * \param tree a T& argument, represents the tree where the statistics are executed
     * \param boxes a set containing the boxes for the queries
     * \param stats a Statistics& argument, representing the object for computing the associated statistics
     */
    template<class T> void exec_box_queries(T& tree, set<Box> &boxes, Statistics &stats);
    ///A public method that excutes incremental nearest neighbor queries
    /*!
     * This method prints the results on standard output
     * The algorithm implemented follows the procedure described in:
     * "Ranking in spatial databases", G. R. Hjaltason, H. Samet, 1995, Advances in Spatial Databases
     *
     * \param tree a T& argument, represents the tree where the statistics are executed
     * \param query_path a string argument, representing the file path of the query input
     * \param key a Point_type variable, representing the the nearest key that is searched by the queries
     * \param critical_points, a vector variable, containing the point types
     */
    template<class T> void exec_incremental_nearest_neighbor_queries(T& tree, string query_path, Point_Type key, vector<Point_Type> &critical_points);
    ///A public method that excutes incremental nearest neighbor queries
    /*!
     * This method prints the results on standard output
     * The algorithm implemented follows the procedure described in:
     * "Ranking in spatial databases", G. R. Hjaltason, H. Samet, 1995, Advances in Spatial Databases
     *
     * \param tree a T& argument, represents the tree where the statistics are executed
     * \param points a set containing the points for the queries
     * \param key a Point_type variable, representing the the nearest key that is searched by the queries
     * \param critical_points, a vector variable, containing the point types
     */
    template<class T> void exec_incremental_nearest_neighbor_queries(T& tree, set<Point> &points, Point_Type key,
                                                                     vector<Point_Type> &critical_points);

    template<class T> bool exec_point_interpolation(T& tree, Point& p, Mesh& mesh, Spatial_Subdivision& division, coord_type& result);
private:
    ///A private method that executes a single point location on a Terrain tree
    /*!
     * \param n a N& argument, representing the current node to visit
     * \param dom a Box& argument, representing the node domain
     * @param level an integer argument representing the level of n in the hierarchy
     * \param p a Point& argument, representing the point query
     * \param qS a QueryStatistics& argument, representing the object in which the statistics are saved
     * \param mesh a Mesh& argument, representing the current mesh
     * \param division a Spatial_Subdivision& argument, representing the tree subdivision
     */
    template<class N> void exec_point_query(N& n, Box& dom, int level, Point& p, QueryStatistics& qS, Mesh& mesh, Spatial_Subdivision& division);
    ///A private method that executes a single box query on a Terrain tree
    /*!
     * \param n a N& argument, representing the actual node to visit
     * \param dom a Box& argument, representing the node domain
     * @param level an integer argument representing the level of n in the hierarchy
     * \param b a Box& argument, representing the box query
     * \param qS a QueryStatistics& argument, representing the object in which the statistics are saved
     * \param mesh a Mesh& argument, representing the current mesh
     * \param division a Spatial_Subdivision& argument, representing the tree subdivision type
     */
    template<class N> void exec_box_query(N& n, Box& dom, int level, Box& b, QueryStatistics& qS, Mesh& mesh, Spatial_Subdivision& division,
                                          bool get_stats);

    ///A private method that executes a point location in a leaf block
    /*!
     * \param n a N& argument, representing the current leaf
     * \param p a Point& argument, representing the point query
     * \param qS a QueryStatistics& argument, representing the object in which the statistics are saved
     * \param mesh a Mesh& argument, representing the current mesh
     */
    template<class N> void exec_point_query_leaf(N& n, Point& p, QueryStatistics& qS, Mesh& mesh);
    /**
     * @brief A private method executing a point-in-triangle test on a triangle
     * @param tet_id an integer representing the current triangle to test
     * @param p a Point& argument, representing the point query
     * @param qS a QueryStatistics& argument, representing the object in which the statistics are saved
     * @param mesh a Mesh& argument, representing the current mesh
     * @return true if an intersection exists, false otherwise
     */
    bool atomic_point_in_triangle_test(itype t_id, Point& p, QueryStatistics& qS, Mesh& mesh);
    ///A private method that executes a box query in a leaf
    /*!
     * \param n a N& argument, representing the actual leaf
     * \param b a Box& argument, representing the box query
     * \param qS a QueryStatistics& argument, representing the object in which the statistics are saved
     * \param mesh a Mesh& argument, representing the current mesh
     * \param get_stats a boolean, true if statistics must be computed, false otherwise
     */
    template<class N> void exec_box_query_leaf_test(N& n, Box& b, QueryStatistics& qS, Mesh& mesh, bool get_stats);
    /**
     * @brief A private method executing a triangle-in-box test on a triangle
     *
     * @param t_id an integer representing the current triangle to test
     * @param b a Box& argument, representing the box query
     * @param qS a QueryStatistics& argument, representing the object in which the statistics are saved
     * @param mesh a Mesh& argument, representing the current mesh
     * @param get_stats a boolean, true if statistics must be computed, false otherwise
     */
    void atomic_triangle_in_box_test(itype t_id, Box& b, QueryStatistics& qS, Mesh& mesh, bool get_stats);
    /**
     * @brief A private method that adds the triangles in a leaf to the result set
     * NOTA: this procedures simply add all the tetrahedra as the domain of the leaf is completely contained by the box QueryStatistics
     *
     * @param n a N& argument, representing the actual leaf
     * @param qS a QueryStatistics& argument, representing the object in which the statistics are saved
     * @param get_stats a boolean, true if statistics must be computed, false otherwise
     */
    template<class N> void add_tetrahedra_to_box_query_result(N& n, QueryStatistics& qS, bool get_stats);

    // incremental nearest neighbor auxiliary procedures
    ///A public method that excutes single incremental nearest neighbor queries
    /*!
     * The algorithm implemented follows the procedure described in:
     * "Ranking in spatial databases", G. R. Hjaltason, H. Samet, 1995, Advances in Spatial Databases
     *
     * @param root a N& argument, representing the actual leaf
     * \param dom a Box& argument, representing the node domain
     * @param p a Point& argument, representing the point query
     * @param mesh a Mesh& argument, representing the current mesh
     * \param division a Spatial_Subdivision& argument, representing the tree subdivision type
     * \param key a Point_type variable, representing the the nearest key that is searched by the queries
     * \param critical_points, a vector variable, containing the point types
     */
    template<class N> int exec_incremental_nearest_neighbor_query(N &root, Box& dom, Point& p, Mesh& mesh, Spatial_Subdivision& division,
                                                                  Point_Type &key, vector<Point_Type> &critical_points);
    /**
     * @brief A private procedure that process a leaf block (Node_V version)
     * @param queue a priority queue variable, containing the blocks and vertices in increasing order
     * @param element a PQueue_Element representing the element to process
     * @param p a Point& argument, representing the point query
     * @param mesh a Mesh& argument, representing the current mesh
     */
    void insert_nearest_vertices(priority_iNN_queue<Node_V> &queue, PQueue_Element<Node_V> &element, Point &p, Mesh &mesh);
    /**
     * @brief A private procedure that process a leaf block (Node_T version)
     * @param queue a priority queue variable, containing the blocks and vertices in increasing order
     * @param element a PQueue_Element representing the element to process
     * @param p a Point& argument, representing the point query
     * @param mesh a Mesh& argument, representing the current mesh
     */
    void insert_nearest_vertices(priority_iNN_queue<Node_T> &queue, PQueue_Element<Node_T> &element, Point &p, Mesh &mesh);
};

template<class T> void Spatial_Queries::exec_point_locations(T& tree, string query_path, Statistics &stats)
{
    QueryStatistics qS = QueryStatistics();
    vector<Point> points;
    Reader::read_queries(points,query_path);

    Timer time;
    coord_type tot_time = 0;
    int hit_ratio = 0;

    for(unsigned i=0;i<points.size();i++)
    {
        time.start();
        this->exec_point_query(tree.get_root(),tree.get_mesh().get_domain(),0,points[i],qS, tree.get_mesh(),tree.get_subdivision());
        time.stop();
        tot_time += time.get_elapsed_time();

        //debug print
        if(qS.triangles.size()>0)
            cout<<"found triangle: "<<qS.triangles[0]<<" "<<tree.get_mesh().get_triangle(qS.triangles[0])<<" for point "<<points[i]<<endl;
        else
            cout<<"nothing found for point "<<points[i]<<endl;

        hit_ratio += stats.compute_queries_statistics(qS);
        qS.reset();
    }
    cerr<<"[TIME] exec point locations "<<tot_time<<endl;

    Writer::write_queries_stats(points.size(),stats.get_query_statistics(),hit_ratio);
    points.clear();
}

template<class T> void Spatial_Queries::exec_point_locations(T& tree, set<Point> &points, Statistics &stats)
{
    QueryStatistics qS = QueryStatistics();

    Timer time;
    coord_type tot_time = 0;
    int hit_ratio = 0;

    for(set<Point>::iterator it=points.begin(); it!=points.end(); ++it)
    {
        Point p = *it;
        time.start();
        this->exec_point_query(tree.get_root(),tree.get_mesh().get_domain(),0,p,qS, tree.get_mesh(),tree.get_subdivision());
        time.stop();
        tot_time += time.get_elapsed_time();

        //debug print
        if(qS.triangles.size()>0)
            cout<<"found triangle: "<<qS.triangles[0]<<" "<<tree.get_mesh().get_triangle(qS.triangles[0])<<" for point "<<p<<endl;
        else
            cout<<"nothing found for point "<<p<<endl;

        hit_ratio += stats.compute_queries_statistics(qS);
        qS.reset();
    }
    cerr<<"[TIME] exec point locations "<<tot_time<<endl;

    Writer::write_queries_stats(points.size(),stats.get_query_statistics(),hit_ratio);
}

template<class T> void Spatial_Queries::exec_box_queries(T& tree, string query_path, Statistics &stats)
{
    QueryStatistics qS = QueryStatistics(tree.get_mesh().get_triangles_num(),4);

    vector<Box> boxes;
    Reader::read_queries(boxes,query_path);

    Timer time;
    coord_type tot_time = 0;
    int hit_ratio = 0;

    for(unsigned j=0;j<boxes.size();j++)
    {
        // exec for timings
        time.start();
        this->exec_box_query(tree.get_root(),tree.get_mesh().get_domain(),0,boxes[j],qS, tree.get_mesh(),tree.get_subdivision(),false);
        time.stop();
        tot_time += time.get_elapsed_time();

        // exec again for stats
        qS.reset(false);
        this->exec_box_query(tree.get_root(),tree.get_mesh().get_domain(),0,boxes[j],qS, tree.get_mesh(),tree.get_subdivision(),true);

        //debug print
        cout<<qS.triangles.size()<<" triangles intersect box "<<j<<endl;

        hit_ratio += stats.compute_queries_statistics(qS);
        qS.reset(true);
    }
    cerr<<"[TIME] exec box queries "<<tot_time<<endl;

    Writer::write_queries_stats(boxes.size(),stats.get_query_statistics(),hit_ratio);
    boxes.clear();
}

template<class T> void Spatial_Queries::exec_box_queries(T& tree, set<Box> &boxes, Statistics &stats)
{
    QueryStatistics qS = QueryStatistics(tree.get_mesh().get_triangles_num(),4);

    Timer time;
    coord_type tot_time = 0;
    int hit_ratio = 0;
    int counter = 0;

    for(set<Box>::iterator it=boxes.begin(); it!=boxes.end();++it)
    {
        Box b = *it;
        // exec for timings
        time.start();
        this->exec_box_query(tree.get_root(),tree.get_mesh().get_domain(),0,b,qS, tree.get_mesh(),tree.get_subdivision(),false);
        time.stop();
        tot_time += time.get_elapsed_time();

        // exec again for stats
        qS.reset(false);
        this->exec_box_query(tree.get_root(),tree.get_mesh().get_domain(),0,b,qS, tree.get_mesh(),tree.get_subdivision(),true);

        //debug print
        cout<<qS.triangles.size()<<" triangles intersect box "<<counter<<endl;

        hit_ratio += stats.compute_queries_statistics(qS);
        qS.reset(true);
        counter++;
    }
    cerr<<"[TIME] exec box queries "<<tot_time<<endl;

    Writer::write_queries_stats(boxes.size(),stats.get_query_statistics(),hit_ratio);
}

template<class N> void Spatial_Queries::exec_point_query(N &n, Box &dom, int level, Point &p, QueryStatistics &qS, Mesh &mesh, Spatial_Subdivision &division)
{
    qS.numNode++;

    if (n.is_leaf())
    {
        qS.numLeaf++;
        this->exec_point_query_leaf(n,p,qS,mesh);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(dom,level,i);
            int son_level = level +1;
            if(son_dom.contains(p,mesh.get_domain().get_max()))
            {
                this->exec_point_query(*n.get_son(i),son_dom,son_level,p,qS,mesh,division);
                break;
            }
        }
    }
}

template<class N> void Spatial_Queries::exec_point_query_leaf(N &n, Point &p, QueryStatistics &qS, Mesh &mesh)
{
    Box bb;
    pair<itype,itype> run;

    for(ivect_iter it=n.get_t_array_begin(); it!=n.get_t_array_end(); ++it)
    {
        if(n.get_run_bounding_box(it,bb,mesh,run))
        {
            if(bb.contains(p,mesh.get_domain().get_max()))
            {
                for(itype t_id=run.first; t_id<=run.second; t_id++)
                {
                    if(atomic_point_in_triangle_test(t_id,p,qS,mesh))
                        return;
                }
            }
        }
        else
        {
            if(atomic_point_in_triangle_test(*it,p,qS,mesh))
                return;
        }
    }
}

template<class T> bool Spatial_Queries::exec_point_interpolation(T& tree,  Point& p, Mesh& mesh, Spatial_Subdivision& division, coord_type& result){
    QueryStatistics qS;
    this->exec_point_query(tree.get_root(),tree.get_mesh().get_domain(),0,p,qS, tree.get_mesh(), tree.get_subdivision());
    if(qS.triangles.size() == 0) return false;
    result = Geometry_Wrapper::elevation_interpolate_tri(qS.triangles[0], p, tree.get_mesh());
    return true;
}   

template<class N> void Spatial_Queries::exec_box_query(N &n, Box &dom, int level, Box &b, QueryStatistics &qS, Mesh &mesh, Spatial_Subdivision &division, bool get_stats)
{
    if(get_stats)
        qS.numNode++;

    if (!dom.intersects(b))
        return;

    if (n.is_leaf())
    {
        if(get_stats)
            qS.numLeaf++;
        if(b.completely_contains(dom))
        {
            this->add_tetrahedra_to_box_query_result(n,qS,get_stats);
        }
        else
            this->exec_box_query_leaf_test(n,b,qS,mesh,get_stats);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(dom,level,i);
            int son_level = level +1;
            this->exec_box_query(*n.get_son(i), son_dom, son_level, b, qS, mesh,division, get_stats);
        }
    }
}

template<class N> void Spatial_Queries::exec_box_query_leaf_test(N &n, Box &b, QueryStatistics &qS, Mesh &mesh, bool get_stats)
{
    Box bb;
    pair<itype,itype> run;

    for(ivect_iter it=n.get_t_array_begin(); it!=n.get_t_array_end(); ++it)
    {
        if(n.get_run_bounding_box(it,bb,mesh,run))
        {
            if(b.completely_contains(bb))
            {
                for(itype t_id=run.first; t_id<=run.second; t_id++)
                {
                    if(get_stats)
                        qS.increase_tri_counter(t_id);

                    if(!qS.is_checked(t_id))
                    {
                        qS.set_checked(t_id);
                        qS.triangles.push_back(t_id);
                    }
                }

            }
            else if(b.intersects(bb))
            {
                for(int t_id=run.first; t_id<=run.second; t_id++)
                    atomic_triangle_in_box_test(t_id,b,qS,mesh,get_stats);
            }
        }
        else
        {
            atomic_triangle_in_box_test(*it,b,qS,mesh,get_stats);
        }
    }

}

template<class N> void Spatial_Queries::add_tetrahedra_to_box_query_result(N& n, QueryStatistics& qS, bool get_stats)
{
    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        if(get_stats)
            qS.increase_tri_counter(*t_id);

        if(!qS.is_checked(*t_id))
        {
            qS.set_checked(*t_id);
            qS.triangles.push_back(*t_id);
        }
    }
}

template<class T> void Spatial_Queries::exec_incremental_nearest_neighbor_queries(T& tree, string query_path, Point_Type key, vector<Point_Type> &critical_points)
{
    Mesh &mesh = tree.get_mesh();

    vector<Point> points;
    Reader::read_queries(points,query_path);

    Timer time;
    coord_type tot_time = 0;

    for(utype i=0;i<points.size();i++)
    {
        time.start();
        int nearest_vertex = exec_incremental_nearest_neighbor_query(tree.get_root(),mesh.get_domain(), 0,points[i],mesh,tree.get_division(),key,critical_points);
        time.stop();
        tot_time += time.get_elapsed_time();
        cerr<<"nearest critical vertex of point "<<points[i]<<" is vertex "<<nearest_vertex<<" -> "<<mesh.get_vertex(nearest_vertex)<<endl;
    }
    cerr<<"[TIME] exec incremental nearest neighbor searches "<<tot_time<<endl;
    points.clear();
}

template<class T> void Spatial_Queries::exec_incremental_nearest_neighbor_queries(T& tree, set<Point> &points, Point_Type key, vector<Point_Type> &critical_points)
{
    Mesh &mesh = tree.get_mesh();

    Timer time;
    coord_type tot_time = 0;

    for(set<Point>::iterator it=points.begin(); it!=points.end(); ++it)
    {
        Point p = *it;
        time.start();
        int nearest_vertex = exec_incremental_nearest_neighbor_query(tree.get_root(),mesh.get_domain(),p,mesh,tree.get_subdivision(),key,critical_points);
        time.stop();
        tot_time += time.get_elapsed_time();
        cerr<<"nearest critical vertex of point "<<p<<" is vertex "<<nearest_vertex<<" -> "<<mesh.get_vertex(nearest_vertex)
           <<" distance: "<<p.distance(mesh.get_vertex(nearest_vertex))<<endl;
    }
    cerr<<"[TIME] exec incremental nearest neighbor searches "<<tot_time<<endl;
}

template<class N> int Spatial_Queries::exec_incremental_nearest_neighbor_query(N &root, Box& root_dom, Point& p, Mesh& mesh, Spatial_Subdivision &division,
                                                                               Point_Type &key, vector<Point_Type> &critical_points)
{
    priority_iNN_queue<N> queue;
    PQueue_Element<N> elem = PQueue_Element<N>(0.0,&root,0,root_dom);
    queue.push(elem);

    while(!queue.empty())
    {
        PQueue_Element<N> current = queue.top();
        queue.pop();

        if(current.is_block())
        {
            N* n = current.get_node();
            if(n->is_leaf())
            {
                // here we have to distinguish between Node_V and Node_T
                this->insert_nearest_vertices(queue,current,p,mesh);
            }
            else
            {
                // we push in queue the sons
                for (int i = 0; i < division.son_number(); i++)
                {
                    Box son_dom = division.compute_domain(current.get_domain(),current.get_value(),i);
                    elem.set_domain(son_dom);
                    elem.set_value(current.get_value() +1);
                    if(son_dom.contains(p,root_dom.get_max()))
                        elem.set_distance(0.0);
                    else
                        elem.set_distance(son_dom.min_distance(p));
                    elem.set_node(n->get_son(i));
                    queue.push(elem);
                }
            }
        }
        else
        {
            /// here we have two steps.. [1] check if the point is a critical point (minimum/saddle/maximum --> this depend on the value of key valiable)
            ///    if so return it as the solution
            /// otherwise if not remove all the duplicate (there are any of them?)
            if(critical_points[current.get_value()-1] == key)
                return current.get_value();
            else
            {
                while(1)
                {
                    elem = queue.top();
                    if(current.get_value() == elem.get_value())
                        queue.pop();
                    else
                        break;
                }
            }
        }
    }
}

#endif // SPATIAL_QUERIES_H

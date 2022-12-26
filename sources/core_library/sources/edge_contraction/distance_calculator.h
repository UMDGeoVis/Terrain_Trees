#ifndef DISTANCE_CALCULATOR_H
#define DISTANCE_CALCULATOR_H

#pragma once

#include "terrain_trees/node_v.h"
#include "terrain_trees/node_t.h"
#include "terrain_trees/prt_tree.h"
#include "statistics/statistics.h"
#include "utilities/usage.h"
#include "utilities/cli_parameters.h"
#include "queries/topological_queries.h"
#include "simplification_aux_structures.h"
#include "utilities/timer.h"
#include <limits>

class Distance_Calculator
{
public:
    Distance_Calculator();
    ~Distance_Calculator();
    template<class T> void vertical_distance(T& tree, vector<Vertex> &points);
    inline void print_stats(){    
        cout<<"avg: "<<stat_avg<<endl;
        cout<<"max: "<<stat_max<<endl;
        };
private:
    template<class N> coord_type point_interpolate(N& n, Box& dom, int level, Point& p, Mesh& mesh, Spatial_Subdivision& division, bool& external);
    template<class N> coord_type point_interpolate_leaf(N &n, Point &p, Mesh &mesh, bool& external);

    inline coord_type interpolate_triangle(int t_id, Point& p, Mesh& mesh){
        Triangle& t = mesh.get_triangle(t_id);
        Vertex v1 = mesh.get_vertex(t.TV(0));
        Vertex v2 = mesh.get_vertex(t.TV(1));
        Vertex v3 = mesh.get_vertex(t.TV(2));
        coord_type divider = (v2.get_y() - v3.get_y()) * (v1.get_x() - v3.get_x()) + (v3.get_x() - v2.get_x()) * (v1.get_y() - v3.get_y());
        coord_type w1 = (v2.get_y() - v3.get_y()) * (p.get_x() - v3.get_x()) + (v3.get_x() - v2.get_x()) * (p.get_y() - v3.get_y());
        w1 /= divider;
        coord_type w2 = (v3.get_y() - v1.get_y()) * (p.get_x() - v3.get_x()) + (v1.get_x() - v3.get_x()) * (p.get_y() - v3.get_y());
        w2 /= divider;
        coord_type w3 = 1 - w1 - w2;
        coord_type elevation = w1 * v1.get_z() + w2 * v2.get_z() + w3 * v3.get_z();
        return elevation;
    };

    coord_type stat_max;
    coord_type stat_avg;
    coord_type stat_min;
    int count_external;
    coord_type elev_max;
    coord_type elev_min;
};

template<class T> void Distance_Calculator::vertical_distance(T& tree, vector<Vertex> &points){
    coord_type stat_sum = 0;
    elev_max = points[0].get_z();
    elev_min = points[0].get_z();
    #pragma omp parallel for reduction(max:stat_max) reduction(+:stat_sum) reduction(+:count_external) reduction(min:stat_min) reduction(min:elev_min) reduction(max:elev_max)
    for(int i = 0; i < points.size(); i++){
        bool external = false;
        elev_max = points[i].get_z() > elev_max ? points[i].get_z():elev_max;
        elev_min = points[i].get_z() < elev_min ? points[i].get_z():elev_min;
        coord_type z = this->point_interpolate(tree.get_root(),tree.get_mesh().get_domain(),0,points[i], tree.get_mesh(),tree.get_subdivision(), external);
        if(external) {
            count_external++;
            // cout<<points[i].get_x()<<", "<<points[i].get_y()<<endl;
            continue;
        }
        coord_type dist = abs(z - points[i].get_z());
        dist = dist > Zero ? dist:0;
        stat_max = dist > stat_max ? dist:stat_max;
        stat_min = dist < stat_min ? dist:stat_min;
        stat_sum += dist;
        
    }
    
    stat_avg = stat_sum/double(points.size() - count_external);
    cout<<"Original mesh statistics (max, min): " <<endl;
    cout<<elev_max<<"  "<<elev_min<<endl;
    cout<<"=============Vertical distance statistics==========="<<endl;
    cout<<"avg: "<<stat_avg<<endl;
    cout<<"max: "<<stat_max<<endl;
    cout<<"min: "<<stat_min<<endl;
    cout<<"count external: "<<count_external<<endl;
    cout<<"diagonal: "<<tree.get_mesh().get_domain().get_diagonal()<<endl;
    cout<<"max wrt. diagonal "<< stat_max/tree.get_mesh().get_domain().get_diagonal()<<endl;
}

template<class N> coord_type Distance_Calculator::point_interpolate(N &n, Box &dom, int level, Point& p, Mesh &mesh, Spatial_Subdivision &division, bool& external)
{
    if (n.is_leaf())
    {
        return this->point_interpolate_leaf(n,p,mesh,external);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(dom,level,i);
            int son_level = level +1;
            if(son_dom.contains(p,mesh.get_domain().get_max()))
            {
                return this->point_interpolate(*n.get_son(i),son_dom,son_level,p,mesh,division,external);
            }
        }
        external = true;
        return 0;
    }
}

template<class N> coord_type Distance_Calculator::point_interpolate_leaf(N &n, Point &p, Mesh &mesh, bool& external)
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
                    if(Geometry_Wrapper::point_in_triangle(t_id,p,mesh))
                    {
                        coord_type z = interpolate_triangle(t_id, p, mesh);
                        return z;
                    }

                }
            }
        }
        else
        {
            if(Geometry_Wrapper::point_in_triangle(*it,p,mesh))
            {
                coord_type z = interpolate_triangle(*it, p, mesh);
                return z;
            }
        }
    }

    external = true;
    return 0;
}

#endif
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

#include "forman_gradient_features_extractor.h"

#include <fstream>

void Forman_Gradient_Features_Extractor::extract_incidence_graph(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, Spatial_Subdivision &division,
                                                                 OpType operation, int cache_size)
{
    /// reset the statistical counters
    reset_extraction_critica_counters();

    mig_cache cache = mig_cache(cache_size);
    ig_paths paths;

    this->incidence_graph_extraction(n,mesh,gradient,division,n,operation,cache,paths);
    if(!paths.visited_all())
    {
        cout<<"[ERROR] the procedure to not visit fully the dangling paths"<<endl;
    }

    //debug print
    if(operation == OUTPUT)
    {
        forman_ig.print_stats(false);
    }
}

void Forman_Gradient_Features_Extractor::extract_incidence_graph(Node_T &n, Box &n_dom, Mesh &mesh, Forman_Gradient &gradient, Spatial_Subdivision &division, OpType operation, int cache_size)
{
    /// reset the statistical counters
    reset_extraction_critica_counters();

    mig_cache cache = mig_cache(cache_size);
    ig_paths paths;

    this->incidence_graph_extraction(n,n_dom,0,mesh,gradient,division,n,operation,cache,paths);
    if(!paths.visited_all())
    {
        cout<<"[ERROR] the procedure to not visit fully the dangling paths"<<endl;
    }

    //debug print
    if(operation == OUTPUT)
    {
        forman_ig.print_stats(false);
    }
}

void Forman_Gradient_Features_Extractor::incidence_graph_extraction(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, Spatial_Subdivision &division,  Node_V &root,
                                                                    OpType operation, mig_cache &cache, ig_paths &paths)
{
    if (n.is_leaf())
    {
        this->incidence_graph_extraction_leaf(n,mesh,gradient,root,division,operation,cache,paths);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                this->incidence_graph_extraction(*n.get_son(i),mesh,gradient,division,root,operation,cache,paths);
            }
        }
    }
}

void Forman_Gradient_Features_Extractor::incidence_graph_extraction(Node_T &n, Box &n_dom, int n_level, Mesh &mesh, Forman_Gradient &gradient, Spatial_Subdivision &division,  Node_T &root,
                                                                    OpType operation, mig_cache &cache, ig_paths &paths)
{
    if (n.is_leaf())
    {
        this->incidence_graph_extraction_leaf(n,n_dom,mesh,gradient,root,division,operation,cache,paths);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(n_dom,n_level,i);
            int son_level = n_level +1;
            if(n.get_son(i)!=NULL)
            {
                this->incidence_graph_extraction(*n.get_son(i),son_dom,son_level,mesh,gradient,division,root,operation,cache,paths);
            }
        }
    }
}

void Forman_Gradient_Features_Extractor::incidence_graph_extraction_leaf(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, Node_V &root, Spatial_Subdivision &division,
                                                                         OpType operation, mig_cache &cache, ig_paths &paths)
{
    /// if there are no vertices in the leaf we have nothing to do..
    if(!n.indexes_vertices())
        return;
    Timer time;
    local_VTstar_ET local_rels;
    if(operation == TIME_VERBOSE)
        time.start();
    Forman_Gradient_Topological_Relations::get_VTstar_ET(local_rels,n,mesh,gradient);
    if(operation == TIME_VERBOSE)
    {
        time.stop();
        rels_gathering += time.get_elapsed_time();
        time.start();
    }
    get_new_IG_paths(n,mesh,gradient,forman_ig,local_rels,root,division,operation,cache,paths);
    if(operation == TIME_VERBOSE)
    {
        time.stop();
        new_paths += time.get_elapsed_time();
        time.start();
    }
    get_dangling_IG_paths(n,mesh,gradient,forman_ig,local_rels,root,division,operation,cache,paths);

    if(operation == TIME_VERBOSE)
    {
        time.stop();
        u_paths += time.get_elapsed_time();
        time.start();
    }
    //add in cache the local topological data structres
    if(mesh.get_vertices_num() > n.get_v_end())
    {
        cache.addItem(n.get_v_end()+n.get_v_start(),local_rels);
    }
    if(operation == TIME_VERBOSE)
    {
        time.stop();
        cache_t += time.get_elapsed_time();
    }

    //get stats
    if(operation == OUTPUT)
    {
        stats.set_cache_stats(cache);
        stats.set_all_dangling_paths(paths);
    }
}

void Forman_Gradient_Features_Extractor::incidence_graph_extraction_leaf(Node_T &n, Box &n_dom, Mesh &mesh, Forman_Gradient &gradient, Node_T &root, Spatial_Subdivision &division,
                                                                         OpType operation, mig_cache &cache, ig_paths &paths)
{
    itype v_start;
    itype v_end;

    n.get_v_range(v_start,v_end,n_dom,mesh); // we need to gather the vertices range..

    if(v_start == v_end) //no internal vertices..
        return;

    Timer time;
    local_VTstar_ET local_rels;
    if(operation == TIME_VERBOSE)
        time.start();
    Forman_Gradient_Topological_Relations::get_VTstar_ET(local_rels,n,v_start,v_end,mesh,gradient);
    if(operation == TIME_VERBOSE)
    {
        time.stop();
        rels_gathering += time.get_elapsed_time();
        time.start();
    }
    get_new_IG_paths(n,v_start,v_end,mesh,gradient,forman_ig,local_rels,root,division,operation,cache,paths);
    if(operation == TIME_VERBOSE)
    {
        time.stop();
        new_paths += time.get_elapsed_time();
        time.start();
    }
    get_dangling_IG_paths(n,v_start,v_end,mesh,gradient,forman_ig,local_rels,root,division,operation,cache,paths);

    if(operation == TIME_VERBOSE)
    {
        time.stop();
        u_paths += time.get_elapsed_time();
        time.start();
    }
    //add in cache the local topological data structres
    if(mesh.get_vertices_num() > v_end)
    {
        cache.addItem(v_end+v_start,local_rels);
    }
    if(operation == TIME_VERBOSE)
    {
        time.stop();
        cache_t += time.get_elapsed_time();
    }

    //get stats
    if(operation == OUTPUT)
    {
        stats.set_cache_stats(cache);
        stats.set_all_dangling_paths(paths);
    }
}

void Forman_Gradient_Features_Extractor::get_new_IG_paths(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, IG &ig, local_VTstar_ET &local_rels, Node_V &root,
                                                          Spatial_Subdivision &division, OpType operation, mig_cache &cache, ig_paths &paths)
{
    for(leaf_ET::iterator it_e=local_rels.begin_ETs(); it_e!=local_rels.end_ETs(); ++it_e)
    {        
        if(gradient.is_edge_critical(it_e->first,it_e->second,mesh))
        {
            /// intro stuff.. ///
            itype max_field_v = get_max_elevation_vertex(it_e->first);

            Triangle& t_first = mesh.get_triangle(it_e->second.first);

            /// edge and edge_vertices position in the two triangles
            short e_pos_first = t_first.edge_index(it_e->first);
            short e_pos_second = -1;

            short v1_pos = t_first.vertex_index(it_e->first[0]);
            short v2_pos = t_first.vertex_index(it_e->first[1]);

            if(it_e->second.second != -1)
            {
                Triangle& t_second = mesh.get_triangle(it_e->second.second);
                e_pos_second = t_second.edge_index(it_e->first);
            }

            iNode* node = new iNode(max_field_v);

            /// for output and simulation purpose ///
            if(operation == OUTPUT)
            {
                /// for the saddle we are sure that we insert these directly and one time
                /// because each saddle is encoded in one leaf node
                if(it_e->second.second != -1)
                {
                    ig.add_saddle(it_e->second,node);
                    node->add_edge_id(it_e->second.first,it_e->second.second);
                }
                else
                {
                    pair<itype,itype> p = make_pair(it_e->second.first, -e_pos_first - 1);
                    ig.add_saddle(p,node);
                    node->add_edge_id(p.first,p.second);
                }
            }

            /// visit the descending 1 cells ///
            /// NEW --> we visit the two extremes with two distinct visits
            Edge *next_e1 = NULL;
            if(get_next_edge(n,it_e->first[0],next_e1,local_rels.get_VTstars(),cache.get_vtstar_cache(),root,division,mesh,gradient))
                explore_desc1cell_mig(n, next_e1, it_e->first[0], v1_pos, mesh, gradient, ig, local_rels, node, root, division,
                        operation, cache, paths);
            else if(operation == OUTPUT) /// I have already found a complete path (only save it if we will output a vtk file)
                pair_minimum_saddle(node,v1_pos,it_e->first[0],it_e->first[0],ig);

            Edge *next_e2 = NULL;
            if(get_next_edge(n,it_e->first[1],next_e2,local_rels.get_VTstars(),cache.get_vtstar_cache(),root,division,mesh,gradient))
                explore_desc1cell_mig(n, next_e2, it_e->first[1], v2_pos, mesh, gradient, ig, local_rels, node, root, division,
                        operation, cache, paths);
            else if(operation == OUTPUT) /// I have already found a complete path (only save it if we will output a vtk file)
                pair_minimum_saddle(node,v2_pos,it_e->first[1],it_e->first[1],ig);

            /// visit the ascending 1 cells ///
            pair<itype,short> coppia = make_pair(it_e->second.first,e_pos_first);
            explore_asc1cell_mig(n, coppia, node, it_e->second.first, local_rels.get_ETs(),
                                 mesh, gradient, ig, root, division,operation, paths.get_asc1cells_path(), cache.get_et_cache());
            if(it_e->second.second != -1)
            {
                pair<itype,short> coppia2 = make_pair(it_e->second.second,e_pos_second);
                explore_asc1cell_mig(n, coppia2, node, it_e->second.second, local_rels.get_ETs(),
                                     mesh, gradient, ig, root, division,operation, paths.get_asc1cells_path(), cache.get_et_cache());
            }
        }
    }
}

void Forman_Gradient_Features_Extractor::get_new_IG_paths(Node_T &n, itype v_start, itype v_end, Mesh &mesh, Forman_Gradient &gradient, IG &ig, local_VTstar_ET &local_rels, Node_T &root,
                                                          Spatial_Subdivision &division, OpType operation, mig_cache &cache, ig_paths &paths)
{
    for(leaf_ET::iterator it_e=local_rels.begin_ETs(); it_e!=local_rels.end_ETs(); ++it_e)
    {
        if(gradient.is_edge_critical(it_e->first,it_e->second,mesh))
        {
            /// intro stuff.. ///
            itype max_field_v = get_max_elevation_vertex(it_e->first);

            Triangle& t_first = mesh.get_triangle(it_e->second.first);

            /// edge and edge_vertices position in the two triangles
            short e_pos_first = t_first.edge_index(it_e->first);
            short e_pos_second = -1;

            short v1_pos = t_first.vertex_index(it_e->first[0]);
            short v2_pos = t_first.vertex_index(it_e->first[1]);

            if(it_e->second.second != -1)
            {
                Triangle& t_second = mesh.get_triangle(it_e->second.second);
                e_pos_second = t_second.edge_index(it_e->first);
            }

            iNode* node = new iNode(max_field_v);

            /// for output and simulation purpose ///
            if(operation == OUTPUT)
            {
                /// for the saddle we are sure that we insert these directly and one time
                /// because each saddle is encoded in one leaf node
                if(it_e->second.second != -1)
                {
                    ig.add_saddle(it_e->second,node);
                    node->add_edge_id(it_e->second.first,it_e->second.second);
                }
                else
                {
                    pair<itype,itype> p = make_pair(it_e->second.first, -e_pos_first - 1);
                    ig.add_saddle(p,node);
                    node->add_edge_id(p.first,p.second);
                }
            }

            /// visit the descending 1 cells ///
            /// NEW --> we visit the two extremes with two distinct visits
            Edge *next_e1 = NULL;
            if(get_next_edge(n,v_start,v_end,it_e->first[0],next_e1,local_rels.get_VTstars(),cache.get_vtstar_cache(),root,division,mesh,gradient))
                explore_desc1cell_mig(n, v_start, v_end, next_e1, it_e->first[0], v1_pos, mesh, gradient, ig, local_rels, node, root, division,
                        operation, cache, paths);
            else if(operation == OUTPUT) /// I have already found a complete path (only save it if we will output a vtk file)
                pair_minimum_saddle(node,v1_pos,it_e->first[0],it_e->first[0],ig);

            Edge *next_e2 = NULL;
            if(get_next_edge(n,v_start,v_end,it_e->first[1],next_e2,local_rels.get_VTstars(),cache.get_vtstar_cache(),root,division,mesh,gradient))
                explore_desc1cell_mig(n, v_start, v_end, next_e2, it_e->first[1], v2_pos, mesh, gradient, ig, local_rels, node, root, division,
                        operation, cache, paths);
            else if(operation == OUTPUT) /// I have already found a complete path (only save it if we will output a vtk file)
                pair_minimum_saddle(node,v2_pos,it_e->first[1],it_e->first[1],ig);

            /// visit the ascending 1 cells ///
            pair<itype,short> coppia = make_pair(it_e->second.first,e_pos_first);
            explore_asc1cell_mig(n, v_start, v_end, coppia, node, it_e->second.first, local_rels.get_ETs(),
                                 mesh, gradient, ig, root, division,operation, paths.get_asc1cells_path(), cache.get_et_cache());
            if(it_e->second.second != -1)
            {
                pair<itype,short> coppia2 = make_pair(it_e->second.second,e_pos_second);
                explore_asc1cell_mig(n, v_start, v_end, coppia2, node, it_e->second.second, local_rels.get_ETs(),
                                     mesh, gradient, ig, root, division,operation, paths.get_asc1cells_path(), cache.get_et_cache());
            }
        }
    }
}

void Forman_Gradient_Features_Extractor::get_dangling_IG_paths(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, IG &ig, local_VTstar_ET &local_rels,
                                                               Node_V &root, Spatial_Subdivision &division, OpType operation, mig_cache &cache, ig_paths &paths)
{
    itype key = n.get_v_start()+n.get_v_end();

    leaves_1_desc_mig_map::iterator it1 = paths.find_desc1cell(key); // Using key to find saved paths
    if(it1 != paths.end_desc1cell())
    {
        set_1paths &unfinished = it1->second;

        for(set_1paths::iterator it_e = unfinished.begin(); it_e != unfinished.end(); ++it_e)
        {
            explore_desc1cell_mig(n, it_e->get_edge(), it_e->get_last_v(), it_e->get_start_v(), mesh, gradient, ig,
                                  local_rels, it_e->get_node(), root, division, operation, cache, paths);
        }
        unfinished.clear();
        paths.erase(it1);
    }

    leaves_1_asc_mig_map::iterator it = paths.find_asc1cells(key);
    if(it != paths.end_asc1cells())
    {
        set_asc_quadruple &pairs = it->second;
        for(set_asc_quadruple::iterator it_p = pairs.begin(); it_p != pairs.end();++it_p)
        {
//            if(n.has_all_vertices_visited(mesh.get_triangle(it_p->get_first())))
                explore_asc1cell_mig(n, it_p->get_pair(), it_p->get_node(), it_p->get_label(), local_rels.get_ETs(),
                                     mesh, gradient, ig, root, division,operation, paths.get_asc1cells_path(), cache.get_et_cache());
//            else
//                cout<<"[ERROR][get_dangling_IG_paths] inserted in wrong place"<<endl;
        }
        pairs.clear();
        paths.erase(it);
    }
}

void Forman_Gradient_Features_Extractor::get_dangling_IG_paths(Node_T &n, itype v_start, itype v_end, Mesh &mesh, Forman_Gradient &gradient, IG &ig, local_VTstar_ET &local_rels,
                                                               Node_T &root, Spatial_Subdivision &division, OpType operation, mig_cache &cache, ig_paths &paths)
{
    itype key = v_start+v_end;

    leaves_1_desc_mig_map::iterator it1 = paths.find_desc1cell(key);
    if(it1 != paths.end_desc1cell())
    {
        set_1paths &unfinished = it1->second;

        for(set_1paths::iterator it_e = unfinished.begin(); it_e != unfinished.end(); ++it_e)
        {
            explore_desc1cell_mig(n, v_start, v_end, it_e->get_edge(), it_e->get_last_v(), it_e->get_start_v(), mesh, gradient, ig,
                                  local_rels, it_e->get_node(), root, division, operation, cache, paths);
        }
        unfinished.clear();
        paths.erase(it1);
    }

    leaves_1_asc_mig_map::iterator it = paths.find_asc1cells(key);
    if(it != paths.end_asc1cells())
    {
        set_asc_quadruple &pairs = it->second;
        for(set_asc_quadruple::iterator it_p = pairs.begin(); it_p != pairs.end();++it_p)
        {
//            if(mesh.get_triangle(it_p->get_first()).maxindex() < v_end)
                explore_asc1cell_mig(n, v_start, v_end, it_p->get_pair(), it_p->get_node(), it_p->get_label(), local_rels.get_ETs(),
                                     mesh, gradient, ig, root, division,operation, paths.get_asc1cells_path(), cache.get_et_cache());
//            else
//                cout<<"[ERROR][get_dangling_IG_paths] inserted in wrong place"<<endl;
        }
        pairs.clear();
        paths.erase(it);
    }
}

void Forman_Gradient_Features_Extractor::explore_desc1cell_mig(Node_V &n, Edge *edge, itype first_v, short label, Mesh &mesh, Forman_Gradient &gradient, IG &ig,
                                                               local_VTstar_ET &local_rels, iNode *saddle_node, Node_V &root,
                                                               Spatial_Subdivision &division, OpType operation, mig_cache &cache, ig_paths &paths)
{
    queue<Edge*> coda;
    coda.push(edge);

    itype last_v = first_v; /// last vertex visited

    while(!coda.empty())
    {
        Edge* e = coda.front();
        coda.pop();

        if(n.has_all_vertices_visited(*e))
        {
            itype v = e->EV(0) == last_v ? e->EV(1) : e->EV(0);
            Edge *edg = NULL;

            if(get_next_edge(n,v,edg,local_rels.get_VTstars(),cache.get_vtstar_cache(),root,division,mesh,gradient))
            {
                if(*edg != *e)
                {
                    coda.push(edg);
                    last_v = v;
                }
                else
                    delete edg;
            }
            else if(operation == OUTPUT)
            {
                pair_minimum_saddle(saddle_node,label,v,last_v,ig);
            }

            delete e;
        }
        else
        {
            Dangling_Paths_Handler::add_dangling_path(root,division,paths.get_desc1cells_path(),e,last_v,label,saddle_node);
        }
    }
}

void Forman_Gradient_Features_Extractor::explore_desc1cell_mig(Node_T &n, itype v_start, itype v_end, Edge *edge, itype first_v, short label, Mesh &mesh, Forman_Gradient &gradient, IG &ig,
                                                               local_VTstar_ET &local_rels, iNode *saddle_node, Node_T &root,
                                                               Spatial_Subdivision &division, OpType operation, mig_cache &cache, ig_paths &paths)
{
    queue<Edge*> coda;
    coda.push(edge);

    itype last_v = first_v; /// last vertex visited

    while(!coda.empty())
    {
        Edge* e = coda.front();
        coda.pop();

        if(e->maxindex() < v_end)
        {
            itype v = e->EV(0) == last_v ? e->EV(1) : e->EV(0);
            Edge *edg = NULL;

            if(get_next_edge(n,v_start,v_end,v,edg,local_rels.get_VTstars(),cache.get_vtstar_cache(),root,division,mesh,gradient))
            {
                if(*edg != *e)
                {
                    coda.push(edg);
                    last_v = v;
                }
                else
                    delete edg;
            }
            else if(operation == OUTPUT)
            {
                pair_minimum_saddle(saddle_node,label,v,last_v,ig);
            }

            delete e;
        }
        else
        {
            Dangling_Paths_Handler::add_dangling_path(root,mesh,division,paths.get_desc1cells_path(),e,last_v,label,saddle_node);
        }
    }
}

/// return true if the edge has been found and thus the vertex is not a minimum, false otherwise
bool Forman_Gradient_Features_Extractor::get_next_edge(Node_V &n, itype v, Edge *&e, leaf_VTstar &all_vtstar, vtstar_cache &cache, Node_V &root,
                                                       Spatial_Subdivision &division, Mesh &mesh, Forman_Gradient &gradient)
{
    itype vtstar = Forman_Gradient_Topological_Relations::get_VTstar(n,v,all_vtstar,cache,root,division,mesh,gradient);

    if(!gradient.is_vertex_critical(v,vtstar,mesh)) //*it is paired
    {
        Triangle &t = mesh.get_triangle(vtstar);
        short v1i = gradient.convert_compressed_to_expand(vtstar).get_vertex_pair(t.vertex_index(v));
        itype v1 = t.TV(v1i);
        e = new Edge(v,v1);
        return true;
    }
    return false;
}

/// return true if the edge has been found and thus the vertex is not a minimum, false otherwise
bool Forman_Gradient_Features_Extractor::get_next_edge(Node_T &n, itype v_start, itype v_end, itype v, Edge *&e, leaf_VTstar &all_vtstar, vtstar_cache &cache, Node_T &root,
                                                       Spatial_Subdivision &division, Mesh &mesh, Forman_Gradient &gradient)
{
    itype vtstar = Forman_Gradient_Topological_Relations::get_VTstar(n,v_start,v_end,v,all_vtstar,cache,root,division,mesh,gradient);

    if(!gradient.is_vertex_critical(v,vtstar,mesh)) //*it is paired
    {
        Triangle &t = mesh.get_triangle(vtstar);
        short v1i = gradient.convert_compressed_to_expand(vtstar).get_vertex_pair(t.vertex_index(v));
        itype v1 = t.TV(v1i);
        e = new Edge(v,v1);
        return true;
    }
    return false;
}

void Forman_Gradient_Features_Extractor::explore_asc1cell_mig(Node_V &n, const pair<itype, short> &tf, iNode *saddle_node, itype first_t_id,
                                                              leaf_ET &local_ef, Mesh &mesh, Forman_Gradient &gradient, IG &ig,
                                                              Node_V &root, Spatial_Subdivision &division, OpType operation,
                                                              leaves_1_asc_mig_map &paths, et_cache &cache)
{
    queue<pair<itype, short> > coda;
    coda.push(tf);
    pair<itype, short> current_pair;

    ivect edge;
    itype last_t = tf.first; /// last vertex visited

    while(!coda.empty())
    {
        current_pair = coda.front();
        coda.pop();

        Triangle& t = mesh.get_triangle(current_pair.first);

        if(n.has_all_vertices_visited(t))
        {
            short f1 = gradient.convert_compressed_to_expand(current_pair.first).get_face_pair();

            if(f1 != current_pair.second)
            {
                if(f1 > -1)
                {
                    t.TE(f1,edge);
                    pair<itype,itype> ft = Forman_Gradient_Topological_Relations::get_ET(n,edge,local_ef,cache,root,division,mesh);

                    itype next = (ft.first == current_pair.first) ? ft.second : ft.first;
                    if(next != -1)
                    {
                        f1 = mesh.get_triangle(next).edge_index(edge);
                        coda.push(make_pair(next,f1));
                        last_t = current_pair.first;
                    }
                }
                else if (operation == OUTPUT && gradient.is_triangle_critical(current_pair.first))
                {
                    //cout<<"Triangle "<<current_pair.first<<" is critical. "<<endl;
                    pair_saddle_maximum(saddle_node,first_t_id,current_pair.first,last_t,ig);
                }
            }
        }
        else
        {
            itype max_id = t.maxindex();
            Dangling_Paths_Handler::add_dangling_path(root,division,paths,max_id,current_pair,first_t_id, saddle_node);

        }
    }

}

void Forman_Gradient_Features_Extractor::explore_asc1cell_mig(Node_T &n, itype v_start, itype v_end, const pair<itype, short> &tf, iNode *saddle_node, itype first_t_id,
                                                              leaf_ET &local_ef, Mesh &mesh, Forman_Gradient &gradient, IG &ig, Node_T &root, Spatial_Subdivision &division,
                                                              OpType operation, leaves_1_asc_mig_map &paths, et_cache &cache)
{
    queue<pair<itype, short> > coda;
    coda.push(tf);
    pair<itype, short> current_pair;

    ivect edge;
    itype last_t = tf.first; /// last vertex visited

    while(!coda.empty())
    {
        current_pair = coda.front();
        coda.pop();

        Triangle& t = mesh.get_triangle(current_pair.first);

        if(t.maxindex() < v_end)
        {
            short f1 = gradient.convert_compressed_to_expand(current_pair.first).get_face_pair();

            if(f1 != current_pair.second)
            {
                if(f1 > -1)
                {
                    t.TE(f1,edge);
                    pair<itype,itype> ft = Forman_Gradient_Topological_Relations::get_ET(n,v_start,v_end,edge,local_ef,cache,root,division,mesh);

                    itype next = (ft.first == current_pair.first) ? ft.second : ft.first;
                    if(next != -1)
                    {
                        f1 = mesh.get_triangle(next).edge_index(edge);
                        coda.push(make_pair(next,f1));
                        last_t = current_pair.first;
                    }
                }
                else if (operation == OUTPUT && gradient.is_triangle_critical(current_pair.first))
                {
                    pair_saddle_maximum(saddle_node,first_t_id,current_pair.first,last_t,ig);
                }
            }
        }
        else
        {
            itype max_id = t.maxindex();
            Dangling_Paths_Handler::add_dangling_path(root,mesh,division,paths,max_id,current_pair,first_t_id, saddle_node);

        }
    }

}

void Forman_Gradient_Features_Extractor::getBoundarySimplices(ivect key_id, int dim, vector<vector<int>>& boundary){

    boundary.clear();
    set<vector<int>> bd_set;
    if(dim == 1){
        pair<int, int> et = make_pair(key_id[0], key_id[1]);
        iNode* saddle = forman_ig.find_saddle(et); 
        set<Arc*> &arcs_up = saddle->getArcs(true);
        for(set<Arc*>::const_iterator it=arcs_up.begin(); it!=arcs_up.end(); ++it)
        {
            int minimum = ((*it)->getNode_i())->get_critical_index();
            bd_set.insert(ivect(1,minimum));
        }
    }
    else if(dim == 2){
        nNode* maximum = forman_ig.find_maximum(key_id[0]); 
        set<Arc*> &arcs = maximum->getArcs();
        for(set<Arc*>::const_iterator it=arcs.begin(); it!=arcs.end(); ++it)
        {

            


            pair<int,int> et = ((iNode*)(*it)->getNode_i())->get_edge_id();
            vector<int> saddle = {et.first, et.second};
            bd_set.insert(saddle);
        }
    }
    else if(dim == 0){
        // should not happen if we have checked it before calling this function
        // we return nothing
        return;
    }
    
    boundary = vector<ivect>(bd_set.begin(), bd_set.end());
}
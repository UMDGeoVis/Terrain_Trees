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

#include "forman_gradient_simplifier.h"

void Forman_Gradient_Simplifier::exec_topological_simplification(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, Spatial_Subdivision &division, int cache_size, coord_type persistence)
{
    /// reset the statistical counters
    reset_extraction_critica_counters();
    mig_cache cache = mig_cache(cache_size);
    ig_paths paths;
    Timer time;
    time.start();
    if(!second_stage){
        forman_ig.init();                                                                       /// init again the MIG structures
        this->incidence_graph_extraction(n, mesh, gradient, division, n, OUTPUT, cache, paths); /// we need to encode explicitly the MIG
        if (!paths.visited_all())
        {
            cout << "[ERROR] the labels are not correctly assigned." << endl;
        }
    }
    cache.reset();
    
    print_arc_nums();
    // priority_arcs_queue queue;
    regular_queue = priority_arcs_queue();
    second_stage_queue = priority_arcs_queue_second();
    this->build_persistence_queue(n, division, mesh, gradient, cache, persistence);
    if (!second_stage)
    {
        cout << "Queue of simplices is built, its size is " << regular_queue.size() << endl;

        while (!regular_queue.empty())
        {
            Topo_Sempl sempl = regular_queue.top();
            regular_queue.pop();
            if (sempl.arc->getLabel() != 1)
            {
                if (sempl.arc->getLabel() == -1)
                    delete sempl.arc;
                continue;
            }
            iNode *saddle = NULL;
            if (sempl.lvl == 0)
            {
                saddle = (iNode *)sempl.arc->getNode_j();
            }
            else
            {
                saddle = (iNode *)sempl.arc->getNode_i();
            }
            itype t1 = saddle->get_edge_id().first;
            itype t2 = saddle->get_edge_id().second;

            ivect critical_edge;
            if (t2 < 0)
            {
                mesh.get_triangle(t1).TE(-t2 - 1, critical_edge);
            }
            else
            {
                Triangle &second = mesh.get_triangle(t2);
                for (int i = 0; i < 3; i++)
                {
                    if (!(second.has_vertex(mesh.get_triangle(t1).TV(i))))
                    {
                        mesh.get_triangle(t1).TE(i, critical_edge);
                        break;
                    }
                }
            }
            this->topological_simplification(critical_edge, sempl, n, mesh, gradient, division, n, cache, persistence);
        }
    }
    else
    {
        cout << "Queue of simplices is built, its size is " << second_stage_queue.size() << endl;
        while (!second_stage_queue.empty())
        {
            Topo_Sempl sempl = second_stage_queue.top();
            second_stage_queue.pop();
            iNode *saddle = NULL;
            if (sempl.lvl == 0)
            {
                saddle = (iNode *)sempl.arc->getNode_j();
            }
            else
            {
                saddle = (iNode *)sempl.arc->getNode_i();
            }
            itype t1 = saddle->get_edge_id().first;
            itype t2 = saddle->get_edge_id().second;

            if (sempl.arc->getLabel() != 1)
            {
                if (sempl.arc->getLabel() == -1)
                    delete sempl.arc;
                continue;
            }
            ivect critical_edge;
            if (t2 < 0)
            {
                mesh.get_triangle(t1).TE(-t2 - 1, critical_edge);
            }
            else
            {
                Triangle &second = mesh.get_triangle(t2);
                for (int i = 0; i < 3; i++)
                {
                    if (!(second.has_vertex(mesh.get_triangle(t1).TV(i))))
                    {
                        mesh.get_triangle(t1).TE(i, critical_edge);
                        break;
                    }
                }
            }
            this->topological_simplification(critical_edge, sempl, n, mesh, gradient, division, n, cache, persistence);
        }
    }

    time.stop();
    time.print_elapsed_time("TOT MIG + topological simplification time: ");
    bool sea_ice_mode = true;
    if (!sea_ice_mode || !this->second_stage)
    {
        cout <<"simplified topo: "<<refined_topo<<endl;
        return;
    }
    /// FOR SEA ICE ONLY
    // vector<Arc *> valid_ridge_paths;
    set<Arc *> &arcs = forman_ig.getLevelArcs(1);

    for (set<Arc *>::iterator it = arcs.begin(); it != arcs.end(); ++it)
    {
        if ((*it)->getLabel() != -1)
        {
            Vertex &v1 = mesh.get_vertex((*it)->getNode_i()->get_critical_index());
            Vertex &v2 = mesh.get_vertex(get_max_elevation_vertex(mesh.get_triangle((*it)->getNode_j()->get_critical_index())));
            if (v2.get_z() >= 0.6)
            {
                // cout << "add a path to the output"<<endl;
                // valid_ridge_paths.push_back(*it);
                extract_ridge_paths(n, division, *it, mesh, gradient, cache, n, (*it)->getNode_j()->get_critical_index(), false);
            }
        }
    }
    cout <<"simplified topo: "<<refined_topo<<endl;
    cout << "simplified in second stage: "<<second_stage_counter<<endl;

    /// after the topological simplification we check if we have some MIG arc that have a persistence below the target
    /// and could be simplified (i.e. saddle arcs == 2)
    /// for debug

    /////// ==== DISABLED FOR SEA ICE APPLICATION ==== /////
    // {
    // int not_simpl = 0;
    // for (int i = 0; i < 2; i++)
    // {
    //     for (set<Arc *>::iterator it = forman_ig.getLevelArcs(i).begin(); it != forman_ig.getLevelArcs(i).end(); ++it)
    //     {
    //         if ((*it)->getLabel() == 1)
    //         {
    //             Vertex &v1 = mesh.get_vertex((*it)->getNode_i()->get_critical_index());
    //             if (i == 1)
    //             {
    //                 Vertex &v2 = mesh.get_vertex(get_max_elevation_vertex(mesh.get_triangle((*it)->getNode_j()->get_critical_index())));
    //                 if (abs(v1.get_z() - v2.get_z()) <= persistence)
    //                 {
    //                     iNode *saddle = (iNode *)(*it)->getNode_i();
    //                     if (saddle->getArcs(false).size() == 2)
    //                     {
    //                         cout << saddle->get_edge_id().first << ", " << saddle->get_edge_id().second << endl;
    //                         not_simpl++;
    //                         cout << " Arc not simplified " << (*it)->getLabel() << endl;
    //                         cout << *(*it) << endl;
    //                         cout << *it << endl;
    //                     }
    //                 }
    //             }
    //             else
    //             {
    //                 Vertex &v2 = mesh.get_vertex((*it)->getNode_j()->get_critical_index());
    //                 if (abs(v1.get_z() - v2.get_z()) <= persistence)
    //                 {
    //                     iNode *saddle = (iNode *)(*it)->getNode_j();
    //                     if (saddle->getArcs(true).size() == 2)
    //                     {
    //                         cout << saddle->get_edge_id().first << ", " << saddle->get_edge_id().second << endl;
    //                         not_simpl++;
    //                         cout << " Arc not simplified " << (*it)->getLabel() << endl;
    //                         cout << *(*it) << endl;
    //                         cout << *it << endl;
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }
    // if (not_simpl > 0)
    //     cerr << "SIMPLIFICATION WARNING: " << not_simpl << " arcs can be simplified, but are not.." << endl;
    // }
}

void Forman_Gradient_Simplifier::topological_simplification(const ivect &critical_edge, Topo_Sempl &sempl, Node_V &n, Mesh &mesh, Forman_Gradient &gradient, Spatial_Subdivision &division, Node_V &root, mig_cache &cache, coord_type persistence)
{
    if (n.is_leaf())
    {
        /// if there are no vertices in the leaf we have nothing to do..
        if (!n.indexes_vertices() || !n.indexes_vertex(critical_edge[1]))
            return;
        this->simplify(critical_edge, sempl, n, mesh, gradient, cache, root, division, persistence);
    }
    else
    {
        if (!n.indexes_vertex(critical_edge[1]))
            return;
        for (int i = 0; i < division.son_number(); i++)
        {
            if (n.get_son(i) != NULL)
            {
                topological_simplification(critical_edge, sempl, *n.get_son(i), mesh, gradient, division, root, cache, persistence);
            }
        }
    }
}

void Forman_Gradient_Simplifier::simplify(const ivect &critical_edge, Topo_Sempl &sempl, Node_V &n, Mesh &mesh, Forman_Gradient &gradient, mig_cache &cache, Node_V &root, Spatial_Subdivision &division, coord_type persistence)
{
    iNode *saddle = NULL;
    nNode *extrema = NULL;
    int lvl0 = 0, lvl1 = 0, lvl2 = 0;
    local_VTstar_ET local_rels;
    // Forman_Gradient_Topological_Relations::get_VTstar_ET(local_rels, n, mesh, gradient);
    // cout <<"start simplify"<<endl;
    Forman_Gradient_Topological_Relations::get_VTstar_ET(local_rels, n, mesh, gradient, cache);
    if (sempl.lvl == 0)
    {
        lvl0++;
        saddle = (iNode *)sempl.arc->getNode_j();
        extrema = (nNode *)sempl.arc->getNode_i();
        if (saddle->getArcs(true).size() != 2)
            return;
        // cout<<"[Contraction] persistence value:"<<sempl.val<<" edge filtration: "<<sempl.filt_s0<<", "<<sempl.filt_s1<<"; extreme filtration: "<<sempl.filt_ex[0]<<endl;
        // cout<<sempl.arc->getSimplexi()"Persistence value:"<<sempl.val<<endl;
        contraction(critical_edge, extrema, saddle, mesh, gradient, local_rels, cache, n, root, division, persistence);
        refined_topo++;
    }
    else
    {
        lvl1++;
        saddle = (iNode *)sempl.arc->getNode_i();
        extrema = (nNode *)sempl.arc->getNode_j();

        // cout << "number of arcs: "<< extrema->size()<<endl;
        // for(auto maximum_arc:extrema->getArcs()){
        //     cout << *maximum_arc << "; ";
        // }
        // cout << endl;
        if (saddle->getArcs(false).size() != 2)
            return;
        //    cout<<"[Removal] persistence value:"<<sempl.val<<" edge filtration: "<<sempl.filt_s0<<", "<<sempl.filt_s1<<"; ";
        //    cout<<"extreme filtration: "<<sempl.filt_ex[0]<<", "<<sempl.filt_ex[1]<<", "<<sempl.filt_ex[2]<<endl;

        /// FOR SEA ICE ONLY
        if (this->second_stage)
        {
            vector<Arc *> arcs = saddle->get_vector_Arcs(false);
            nNode *other_extrema;
            itype triangle;
            itype ending_path_simplex;
            Arc *arc_to_add = NULL;
            if (arcs.size() != 2)
            {
                return;
            }
            
            if (arcs[0]->getNode_j() == extrema)
            {
                other_extrema = (nNode *)arcs[1]->getNode_j();
                arc_to_add = arcs[0];
            }
            else
            {
                other_extrema = (nNode *)arcs[0]->getNode_j();
                arc_to_add = arcs[1];
            }


            Vertex &maximum = mesh.get_vertex(get_max_elevation_vertex(mesh.get_triangle(extrema->get_critical_index())));
            Vertex &other_maximum = mesh.get_vertex(get_max_elevation_vertex(mesh.get_triangle(other_extrema->get_critical_index())));
            Vertex &saddle_vertex = mesh.get_vertex(saddle->get_critical_index());
            if(other_maximum.get_z() < 0.6 || other_maximum.get_z() <  maximum.get_z()){
                return;
            }
            if(other_maximum.get_z() < 2 * abs(other_maximum.get_z() - saddle_vertex.get_z())){
                return;
            }
            if(extrema->size() == 1){
                // cout << "Maximum cannot be removed" << endl;
                extract_ridge_paths(n, division, arc_to_add, mesh, gradient, cache, root, other_extrema->get_critical_index(), true);
                return;
            }
               
        }

        removal(critical_edge, extrema, saddle, mesh, gradient, local_rels, cache, n, root, division, persistence);
        refined_topo++;
        second_stage_counter++;
    }
    // cout<<"Level 0:"<<lvl0<<", Level 1:"<<lvl1<<", Level 2:"<<lvl2<<endl;
}

void Forman_Gradient_Simplifier::contraction(const ivect &critical_edge, nNode *extrema, iNode *saddle, Mesh &mesh, Forman_Gradient &gradient, local_VTstar_ET &local_rels, mig_cache &cache, Node_V &n, Node_V &root, Spatial_Subdivision &division, coord_type persistence)
{
    vector<Arc *> arcs = saddle->get_vector_Arcs(true);
    pair<itype, itype> critical_edge_tetra = saddle->get_edge_id();
    // ivect critical_edge;
    Triangle &first = mesh.get_triangle(critical_edge_tetra.first);
    nNode *other_extrema;
    itype vertex, next_vertex;
    itype ending_path_simplex;
    // SETUp
    if (arcs[0]->getNode_i() == extrema)
    {
        other_extrema = (nNode *)arcs[1]->getNode_i();
        next_vertex = first.TV(arcs[0]->getSimplexj());
        vertex = first.TV(arcs[1]->getSimplexj());
        ending_path_simplex = arcs[1]->getSimplexi();
    }
    else
    {
        other_extrema = (nNode *)arcs[0]->getNode_i();
        next_vertex = first.TV(arcs[1]->getSimplexj());
        vertex = first.TV(arcs[0]->getSimplexj());
        ending_path_simplex = arcs[0]->getSimplexi();
    }

    /// ----------------- GRADIENT UPDATES ----------------- ///
    contraction_update_gradient(vertex, extrema->get_critical_index(), next_vertex, critical_edge, mesh, gradient, local_rels, cache, n, root, division);

    /// ----------------- MIG UPDATES ----------------- ///
    arcs[0]->setLabel(-1);
    arcs[1]->setLabel(-1);
    /// remove maxima-saddle arcs
    remove_saddle_arcs(saddle, false, forman_ig);
    /// remove extrema arcs
    remove_extreme_arcs(extrema, saddle, other_extrema, ending_path_simplex, true, n, mesh, persistence);
    /// remove minima-saddle arcs
    remove_saddle_arcs(saddle, true, forman_ig);
    /// remove the saddle and the minimum
    forman_ig.remove_saddle(saddle->get_edge_id(), saddle);
    forman_ig.remove_minimum(extrema->get_critical_index(), extrema);
}

void Forman_Gradient_Simplifier::removal(const ivect &critical_edge, nNode *extrema, iNode *saddle, Mesh &mesh, Forman_Gradient &gradient, local_VTstar_ET &local_rels,
                                         mig_cache &cache, Node_V &n, Node_V &root, Spatial_Subdivision &division, coord_type persistence)
{
    vector<Arc *> arcs = saddle->get_vector_Arcs(false);
    nNode *other_extrema;
    itype triangle;
    itype ending_path_simplex;

    if (arcs[0]->getNode_j() == extrema)
    {
        other_extrema = (nNode *)arcs[1]->getNode_j();
        triangle = arcs[0]->getSimplexi();
        ending_path_simplex = arcs[1]->getSimplexj();
    }
    else
    {
        other_extrema = (nNode *)arcs[0]->getNode_j();
        triangle = arcs[1]->getSimplexi();
        ending_path_simplex = arcs[0]->getSimplexj();
    }

    /// ----------------- GRADIENTE UPDATES ----------------- ///
    removal_update_gradient(triangle, saddle, mesh, gradient, local_rels.get_ETs(), cache, n, root, division);

    /// ----------------- MIG UPDATES ----------------- ///
    arcs[0]->setLabel(-1);
    arcs[1]->setLabel(-1);

    /// remove minima-saddle arcs
    remove_saddle_arcs(saddle, true, forman_ig);
    /// remove extrema arcs
    remove_extreme_arcs(extrema, saddle, other_extrema, ending_path_simplex, false, n, mesh, persistence);
    /// remove maxima-saddle arcs
    remove_saddle_arcs(saddle, false, forman_ig);
    /// remove the saddle and the maximum
    forman_ig.remove_saddle(saddle->get_edge_id(), saddle);
    forman_ig.remove_maximum(extrema->get_critical_index(), extrema);
}

void Forman_Gradient_Simplifier::remove_extreme_arcs(nNode *extrema, iNode *saddle, nNode *other_extrema, itype ending_path_simplex,
                                                     bool is_minimum, Node_V &n, Mesh &mesh, coord_type persistence)
{
    iNode *node_saddle1 = NULL;
    itype starting_path_simplex;
    coord_type val;
    set<Arc *>::iterator it = extrema->begin();
    while (it != extrema->end())
    {
        if (is_minimum)
        {
            node_saddle1 = ((iNode *)(*it)->getNode_j());
            starting_path_simplex = (*it)->getSimplexj();
        }
        else
        {
            node_saddle1 = ((iNode *)(*it)->getNode_i());
            starting_path_simplex = (*it)->getSimplexi();
        }
        (*it)->setLabel(-1);
        forman_ig.removeArc(!is_minimum, *it);
        node_saddle1->removeArc(is_minimum, *it);
        extrema->removeArc(it);
        if (node_saddle1 != saddle)
        {
            Arc *existing_arc = forman_ig.already_connected(other_extrema, node_saddle1);
            if (existing_arc == NULL)
            {
                Arc *arco = NULL;

                if (is_minimum)
                    arco = forman_ig.addArc(other_extrema, ending_path_simplex, node_saddle1, starting_path_simplex, 0);
                else
                    arco = forman_ig.addArc(node_saddle1, starting_path_simplex, other_extrema, ending_path_simplex, 1);

                /// the arc can be simplified AND the saddle is indexed into a leaf node that has already been visited
                /// if the saddle will be processed in a successive leaf node
                if (arco->getLabel() == 1)
                {
                    itype index_i, index_j;
                    itype filt_i, filt_j;
                    ivect filt_ex;
                    if (is_minimum)
                    {
                        index_i = arco->getNode_i()->get_critical_index();
                        index_j = arco->getNode_j()->get_critical_index();
                        filt_ex.push_back(filtration[index_i - 1]);
                        pair<itype, itype> critical_edge_tetra = ((iNode *)arco->getNode_j())->get_edge_id();
                        ivect critical_edge;
                        Triangle &first = mesh.get_triangle(critical_edge_tetra.first);
                        if (critical_edge_tetra.second < 0)
                        {
                            first.TE(-critical_edge_tetra.second - 1, critical_edge);
                        }
                        else
                        {
                            Triangle &second = mesh.get_triangle(critical_edge_tetra.second);
                            for (int i = 0; i < 3; i++)
                            {
                                if (!(second.has_vertex(first.TV(i))))
                                {
                                    first.TE(i, critical_edge);
                                    break;
                                }
                            }
                        }
                        filt_i = (filtration[critical_edge[0] - 1] > filtration[critical_edge[1] - 1]) ? filtration[critical_edge[0] - 1] : filtration[critical_edge[1] - 1];
                        filt_j = (filtration[critical_edge[0] - 1] > filtration[critical_edge[1] - 1]) ? filtration[critical_edge[1] - 1] : filtration[critical_edge[0] - 1];
                    }
                    else
                    {
                        index_i = arco->getNode_i()->get_critical_index();
                        pair<itype, itype> critical_edge_tetra = ((iNode *)arco->getNode_i())->get_edge_id();
                        ivect critical_edge;
                        Triangle &first = mesh.get_triangle(critical_edge_tetra.first);
                        if (critical_edge_tetra.second < 0)
                        {
                            first.TE(-critical_edge_tetra.second - 1, critical_edge);
                        }
                        else
                        {
                            Triangle &second = mesh.get_triangle(critical_edge_tetra.second);
                            for (int i = 0; i < 3; i++)
                            {
                                if (!(second.has_vertex(first.TV(i))))
                                {
                                    first.TE(i, critical_edge);
                                    break;
                                }
                            }
                        }
                        // filt_i=filtration[critical_edge[0]-1];
                        // filt_j=filtration[critical_edge[1]-1];
                        filt_i = (filtration[critical_edge[0] - 1] > filtration[critical_edge[1] - 1]) ? filtration[critical_edge[0] - 1] : filtration[critical_edge[1] - 1];
                        filt_j = (filtration[critical_edge[0] - 1] > filtration[critical_edge[1] - 1]) ? filtration[critical_edge[1] - 1] : filtration[critical_edge[0] - 1];
                        Triangle t = mesh.get_triangle(arco->getNode_j()->get_critical_index());
                        for (int i = 0; i < 3; i++)
                            filt_ex.push_back(filtration[t.TV(i) - 1]);
                        index_j = get_max_elevation_vertex(t);
                    }
                    val = abs(mesh.get_vertex(index_i).get_z() - mesh.get_vertex(index_j).get_z());

                    // === SEA ICE ONLY VERSION ===
                    if (this->second_stage)
                    {
                        if (!is_minimum)
                        {
                            // cout<<"Check if updated arc should be added to the queue"<<endl;
                            vector<Arc *> arcs_of_saddle1 = node_saddle1->get_vector_Arcs(false);
                            Vertex maximum_to_connect; // The other maximum connected to the new saddle
                            Arc * arc_other_side = NULL; 
                            if (arcs_of_saddle1.size() == 2)
                            {
                                if (arcs_of_saddle1[0] == arco)
                                {
                                    arc_other_side = arcs_of_saddle1[1];
                                    nNode *new_extrema_to_connect = (nNode *)arcs_of_saddle1[1]->getNode_j();
                                    // cout<<"1:";
                                    // cout <<get_max_elevation_vertex(mesh.get_triangle(new_extrema_to_connect->get_critical_index()))<<endl;
                                    maximum_to_connect = mesh.get_vertex(get_max_elevation_vertex(mesh.get_triangle(new_extrema_to_connect->get_critical_index())));
                                }
                                else
                                {
                                    arc_other_side = arcs_of_saddle1[0];
                                    nNode *new_extrema_to_connect = (nNode *)arcs_of_saddle1[0]->getNode_j();
                                    // cout<<"2: ";
                                    // cout <<get_max_elevation_vertex(mesh.get_triangle(new_extrema_to_connect->get_critical_index()))<<endl;
                                    maximum_to_connect = mesh.get_vertex(get_max_elevation_vertex(mesh.get_triangle(new_extrema_to_connect->get_critical_index())));
                                }

                                sort(filt_ex.begin(), filt_ex.end(), greater<int>());
                                Topo_Sempl ts = Topo_Sempl(arco, maximum_to_connect.get_z(), !is_minimum, filt_i, filt_j, filt_ex);
                                second_stage_queue.push(ts);

                                // coord_type dif_other_side = abs(maximum_to_connect.get_z() - mesh.get_vertex(index_i).get_z());
                                // // In an update case, it is possible that either maximum is higher. So we need to handle both cases. 
                                // if ((maximum_to_connect.get_z() > mesh.get_vertex(index_j).get_z() && maximum_to_connect.get_z() >= 0.6
                                //  && maximum_to_connect.get_z() > 2* dif_other_side)){

                                //     sort(filt_ex.begin(), filt_ex.end(), greater<int>());
                                //     Topo_Sempl ts = Topo_Sempl(arco, maximum_to_connect.get_z(), !is_minimum, filt_i, filt_j, filt_ex);
                                //     second_stage_queue.push(ts);

                                // }else if(maximum_to_connect.get_z() <= mesh.get_vertex(index_j).get_z() && mesh.get_vertex(index_j).get_z() >= 0.6
                                // &&  mesh.get_vertex(index_j).get_z() > 2 * val ){
                                //     Triangle t_other_side = mesh.get_triangle(arc_other_side->getNode_j()->get_critical_index());
                                //     filt_ex.clear();
                                //     for (int i = 0; i < 3; i++)
                                //         filt_ex.push_back(filtration[t_other_side.TV(i) - 1]);
                                //     sort(filt_ex.begin(), filt_ex.end(), greater<int>());

                                //     Topo_Sempl ts = Topo_Sempl(arc_other_side, mesh.get_vertex(index_j).get_z(), 
                                //     !is_minimum, filt_i, filt_j, filt_ex);
                                //     second_stage_queue.push(ts);
                                // }
                            }
                        }
                    }
                    else if (val <= persistence)
                    /// NEW <= instead of < (uniform execution pattern)
                    {
                        /// this must be enable if we simplify only topologically!! (the same holds in the removal function)
                        // third parameter is lvl, minimum: 0 and maximum:1
                        sort(filt_ex.begin(), filt_ex.end(), greater<int>());
                        Topo_Sempl ts = Topo_Sempl(arco, val, !is_minimum, filt_i, filt_j, filt_ex);
                        regular_queue.push(ts);
                    }
                }
            }
            else
            {
                existing_arc->setLabel(2);
                // cout << "Double connected"<<endl;
                // cout << *existing_arc<<endl;
            }
        }

        /// set again to begin
        it = extrema->begin();
    }
}

void Forman_Gradient_Simplifier::build_persistence_queue(Node_V &n, Spatial_Subdivision &division, Mesh &mesh, Forman_Gradient &gradient, mig_cache &cache, coord_type persistence)
{
    if (n.is_leaf())
    {
        /// if there are no vertices in the leaf we have nothing to do..
        if (!n.indexes_vertices())
            return;
        local_VTstar_ET local_rels;
        Forman_Gradient_Topological_Relations::get_VTstar_ET(local_rels, n, mesh, gradient);
        this->build_persistence_queue_leaf(n, local_rels.get_ETs(), mesh, gradient, persistence);
        if (mesh.get_vertices_num() > n.get_v_end())
        {
            cache.addItem(n.get_v_end() + n.get_v_start(), local_rels);
        }
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            if (n.get_son(i) != NULL)
            {
                this->build_persistence_queue(*n.get_son(i), division, mesh, gradient, cache, persistence);
            }
        }
    }
}

/// push into the priority queue only those arc below the avg
void Forman_Gradient_Simplifier::build_persistence_queue_leaf(Node_V &n, leaf_ET &local_et, Mesh &mesh, Forman_Gradient &gradient, coord_type persistence)
{
    for (leaf_ET::iterator it_e = local_et.begin(); it_e != local_et.end(); ++it_e)
    {
        if (!gradient.is_edge_critical(it_e->first, it_e->second, mesh))
        {
            continue;
        }
        /// we get the internal saddle
        iNode *saddle;
        if (it_e->second.second != -1)
            saddle = forman_ig.find_saddle(it_e->second);
        else
        {
            Triangle &top_first = mesh.get_triangle(it_e->second.first);
            short e_pos_first = top_first.edge_index(it_e->first);
            pair<itype, itype> p = make_pair(it_e->second.first, -e_pos_first - 1);
            saddle = forman_ig.find_saddle(p);
        }

        /// and the arcs that link it with minima and maxima
        if (saddle != NULL)
        {
            coord_type val;
            set<Arc *> &arc_down = saddle->getArcs(false); /// we get the arcs from the saddle to the connected maxima
            for (set<Arc *>::iterator it = arc_down.begin(); it != arc_down.end(); ++it)
            {
                if ((*it)->getLabel() == 1)
                {
                    Vertex &v1 = mesh.get_vertex((*it)->getNode_i()->get_critical_index());
                    Vertex &v2 = mesh.get_vertex(get_max_elevation_vertex(mesh.get_triangle((*it)->getNode_j()->get_critical_index())));

                    // elevation difference
                    val = abs(v1.get_z() - v2.get_z());

                    // ====== FOR SEA ICE APPLICATION =====//
                    if (this->second_stage)
                    {
                        vector<Arc *> arcs = saddle->get_vector_Arcs(false);
                        if (arcs.size() != 2)
                        {
                            continue;
                        }
                        Vertex maximum_to_connect;
                        if (arcs[0] == (*it))
                        {
                            nNode *new_extrema_to_connect = (nNode *)arcs[1]->getNode_j();
                            maximum_to_connect = mesh.get_vertex(get_max_elevation_vertex(mesh.get_triangle(new_extrema_to_connect->get_critical_index())));
                        }
                        else
                        {
                            nNode *new_extrema_to_connect = (nNode *)arcs[0]->getNode_j();
                            maximum_to_connect = mesh.get_vertex(get_max_elevation_vertex(mesh.get_triangle(new_extrema_to_connect->get_critical_index())));
                        }
                        // cout<<"Find the other side"<<endl;
                        // if (maximum_to_connect.get_z() < v2.get_z())
                        // {
                        //     continue;
                        // }
                        coord_type dif_other_side = abs(maximum_to_connect.get_z() - v1.get_z());
                        // if (maximum_to_connect.get_z() > 2 * dif_other_side && maximum_to_connect.get_z() >= 0.6)
                        // {
                        Triangle t = mesh.get_triangle((*it)->getNode_j()->get_critical_index());
                        // filt0 and filt1 are the filtrations of two vertices of critical edge.
                        int filt0 = (filtration[it_e->first[0] - 1] > filtration[it_e->first[1] - 1]) ? filtration[it_e->first[0] - 1] : filtration[it_e->first[1] - 1];
                        int filt1 = (filtration[it_e->first[0] - 1] > filtration[it_e->first[1] - 1]) ? filtration[it_e->first[1] - 1] : filtration[it_e->first[0] - 1];
                        ivect filt_ex; // filtration of extreme, i.e., minima or maxima
                        for (int i = 0; i < 3; i++)
                        {
                            filt_ex.push_back(filtration[t.TV(i) - 1]);
                        }
                        sort(filt_ex.begin(), filt_ex.end(), greater<int>());
                        // In second stage ts stores the elevation of the maximum on the other side as val to rank.
                        Topo_Sempl ts = Topo_Sempl(*it, maximum_to_connect.get_z(), 1, filt0, filt1, filt_ex);
                        second_stage_queue.push(ts);
                        // }
                    }
                    else if (val <= persistence)
                    {
                        Triangle t = mesh.get_triangle((*it)->getNode_j()->get_critical_index());
                        // filt0 and filt1 are the filtrations of two vertices of critical edge.
                        int filt0 = (filtration[it_e->first[0] - 1] > filtration[it_e->first[1] - 1]) ? filtration[it_e->first[0] - 1] : filtration[it_e->first[1] - 1];
                        int filt1 = (filtration[it_e->first[0] - 1] > filtration[it_e->first[1] - 1]) ? filtration[it_e->first[1] - 1] : filtration[it_e->first[0] - 1];
                        ivect filt_ex; // filtration of extreme, i.e., minima or maxima
                        for (int i = 0; i < 3; i++)
                        {
                            filt_ex.push_back(filtration[t.TV(i) - 1]);
                        }
                        sort(filt_ex.begin(), filt_ex.end(), greater<int>());
                        Topo_Sempl ts = Topo_Sempl(*it, val, 1, filt0, filt1, filt_ex);
                        regular_queue.push(ts);
                    }
                }
            }
            if (!second_stage)
            {
                set<Arc *> &arc_up = saddle->getArcs(true); /// we get the arcs from the connected minima to the saddle
                for (set<Arc *>::iterator it = arc_up.begin(); it != arc_up.end(); ++it)
                {
                    if ((*it)->getLabel() != 1)
                        continue;
                    Vertex &v1 = mesh.get_vertex((*it)->getNode_i()->get_critical_index());
                    Vertex &v2 = mesh.get_vertex((*it)->getNode_j()->get_critical_index());
                    val = abs(v1.get_z() - v2.get_z());
                    if (val <= persistence)
                    {
                        int filt0 = (filtration[it_e->first[0] - 1] > filtration[it_e->first[1] - 1]) ? filtration[it_e->first[0] - 1] : filtration[it_e->first[1] - 1];
                        int filt1 = (filtration[it_e->first[0] - 1] > filtration[it_e->first[1] - 1]) ? filtration[it_e->first[1] - 1] : filtration[it_e->first[0] - 1];
                        ivect filt_ex;
                        filt_ex.push_back(filtration[(*it)->getNode_i()->get_critical_index() - 1]);
                        Topo_Sempl ts = Topo_Sempl(*it, val, 0, filt0, filt1, filt_ex);
                        regular_queue.push(ts);
                    }
                }
            }
        }
        else
        {
            cout << "[build_persistence_queue] something wrong here...." << endl;
            cout << "missing saddle --> edge: " << it_e->first[0] << " " << it_e->first[1];
            cout << " -- et: " << it_e->second.first << " " << it_e->second.second << endl;
            cout << "et0: " << mesh.get_triangle(it_e->second.first) << endl;
            cout << "et1: " << mesh.get_triangle(it_e->second.second) << endl;
            pair<itype, itype> opposite = make_pair(it_e->second.second, it_e->second.first);
            if (forman_ig.find_saddle(opposite) == NULL)
                cout << "     the opposite is not contained..." << endl;
            else
                cout << "     the opposite is contained into the MIG" << endl;
            int a;
            cin >> a;
        }
    }
}

void Forman_Gradient_Simplifier::extract_ridge_paths(Node_V &n, Spatial_Subdivision &division,
                                                     Arc *arc, Mesh &mesh, Forman_Gradient &gradient, mig_cache &cache,
                                                     Node_V &root, itype assigned_maximum, bool add_other_saddle_tri)
{

    itype saddle = arc->getNode_i()->get_critical_index();
    // cout << saddle <<endl;
    if (n.is_leaf())
    {
        /// if there are no vertices in the leaf we have nothing to do..
        if (!n.indexes_vertices() || !n.indexes_vertex(saddle))
            return;
        this->extract_ridge_paths_leaf(n, division, arc, mesh, gradient, cache, root, assigned_maximum, add_other_saddle_tri);
    }
    else
    {
        if (!n.indexes_vertex(saddle))
            return;
        for (int i = 0; i < division.son_number(); i++)
        {
            if (n.get_son(i) != NULL)
            {
                extract_ridge_paths(*n.get_son(i), division, arc, mesh, gradient, cache, root, assigned_maximum, add_other_saddle_tri);
            }
        }
    }
}

void Forman_Gradient_Simplifier::extract_ridge_paths_leaf(Node_V &n, Spatial_Subdivision &division, Arc *arc,
                                                          Mesh &mesh, Forman_Gradient &gradient, mig_cache &cache, Node_V &root, itype assigned_maximum, bool add_other_saddle_tri)
{
    // cout << "extract ridge"<<endl;
    local_VTstar_ET local_rels;
    Forman_Gradient_Topological_Relations::get_VTstar_ET(local_rels, n, mesh, gradient, cache);
    // itype maximum = arc->getNode_j()->get_critical_index();
    itype triangle = arc->getSimplexi();
    iNode *saddle_node = ((iNode *)arc->getNode_i());
    pair<itype, itype> twotriangles = saddle_node->get_edge_id();
    Triangle &tri = mesh.get_triangle(twotriangles.first);
    // ivect pre_edge; // starts from saddle edge
    ivect cur_edge;
    itype the_other_tri;
    vector<itype> path;
    if (add_other_saddle_tri)
    {
        if (twotriangles.second >= 0)
        {
            if (twotriangles.second == triangle)
            {
                path.push_back(twotriangles.first);
            }
            else
            {
                path.push_back(twotriangles.second);
            }
        }
        else
        {
            return;
        }
    }

    //     ivect e;
    //     Triangle &tri2 = mesh.get_triangle(twotriangles.second);
    //     for (int j = 0; j < tri2.vertices_num(); j++)
    //     {
    //         tri.TE(j, e);
    //         if (tri2.has_edge(e))
    //         {
    //             pre_edge = e;
    //             break;
    //         }
    //     }
    // }
    // else
    // {
    //     tri.TE(-twotriangles.second - 1, pre_edge);
    // }

    itype next_triangle = triangle;

    while (!gradient.is_triangle_critical(next_triangle))
    {
        path.push_back(next_triangle);
        Triangle &tri = mesh.get_triangle(next_triangle);
        int i = get_paired_edge(gradient, tri, next_triangle, cur_edge);
        triangle = next_triangle;
        pair<itype, itype> et = Forman_Gradient_Topological_Relations::get_ET(n, cur_edge, local_rels.get_ETs(), cache.get_et_cache(), root, division, mesh);
        next_triangle = (et.first == triangle) ? et.second : et.first;
    }
    path.push_back(next_triangle);
    this->valid_ridge_cells[assigned_maximum].push_back(path);
}


coord_type Forman_Gradient_Simplifier::get_average_persistence_value(IG &ig, Mesh &mesh)
{
    coord_type avg = 0.0;

    for (int i = 0; i < 2; i++)
    {
        set<Arc *> &arcs = ig.getLevelArcs(i);
        for (set<Arc *>::iterator it = arcs.begin(); it != arcs.end(); it++)
        {
            if ((*it)->getLabel() == 1)
            {
                Vertex &v1 = mesh.get_vertex((*it)->getNode_i()->get_critical_index());

                if (i == 1)
                {
                    Vertex &v2 = mesh.get_vertex(get_max_elevation_vertex(mesh.get_triangle((*it)->getNode_j()->get_critical_index())));
                    avg += fabs(v1.get_z() - v2.get_z());
                }
                else
                {
                    Vertex &v2 = mesh.get_vertex((*it)->getNode_j()->get_critical_index());
                    avg += fabs(v1.get_z() - v2.get_z());
                }
            }
        }
    }

    avg /= ig.get_arcs_number();

    /// for debug: check how many are below the average
    cerr << "total arc number: " << ig.get_arcs_number() << endl;
    int arc_below_counter = 0;
    for (int i = 0; i < 2; i++)
    {
        set<Arc *> &arcs = ig.getLevelArcs(i);
        for (set<Arc *>::iterator it = arcs.begin(); it != arcs.end(); it++)
        {
            if ((*it)->getLabel() == 1)
            {
                Vertex &v1 = mesh.get_vertex((*it)->getNode_i()->get_critical_index());

                if (i == 1)
                {
                    Vertex &v2 = mesh.get_vertex(get_max_elevation_vertex(mesh.get_triangle((*it)->getNode_j()->get_critical_index())));
                    if (avg > fabs(v1.get_z() - v2.get_z()))
                        arc_below_counter++;
                }
                else
                {
                    Vertex &v2 = mesh.get_vertex((*it)->getNode_j()->get_critical_index());
                    if (avg > fabs(v1.get_z() - v2.get_z()))
                        arc_below_counter++;
                }
            }
        }
    }
    cerr << "arcs below the average: " << arc_below_counter << endl;
    cerr << "average_persistence: " << avg << endl;

    return avg;
}

void Forman_Gradient_Simplifier::contraction_update_gradient(itype vertex, itype ex_minimum, itype next_vertex, const ivect &critical_edge, Mesh &mesh,
                                                             Forman_Gradient &gradient, local_VTstar_ET &local_rels, mig_cache &cache,
                                                             Node_V &n, Node_V &root, Spatial_Subdivision &division)
{
    ivect old_edge = critical_edge;
    ivect edge;
    pair<itype, itype> old_ef = Forman_Gradient_Topological_Relations::get_ET(n, old_edge, local_rels.get_ETs(), cache.get_et_cache(), root, division, mesh);
    while (next_vertex != ex_minimum)
    {
        // trovo il nuovo edge;
        itype vtstar = Forman_Gradient_Topological_Relations::get_VTstar(n, next_vertex, local_rels.get_VTstars(), cache.get_vtstar_cache(), root, division, mesh, gradient);
        Triangle &t = mesh.get_triangle(vtstar);
        short v1i = gradient.convert_compressed_to_expand(vtstar).get_vertex_pair(t.vertex_index(next_vertex));

        /// no checking vli != -1 because it must be paired
        itype v1 = t.TV(v1i);

        edge = {min(next_vertex, v1), max(next_vertex, v1)};
        vertex = next_vertex;
        next_vertex = v1;

        // azzero la sua adiacenza;
        pair<itype, itype> ef = Forman_Gradient_Topological_Relations::get_ET(n, edge, local_rels.get_ETs(), cache.get_et_cache(), root, division, mesh);
        gradient.free_VE(vertex, next_vertex, ef, mesh);

        // accoppio il vecchio edge;
        if (old_edge[0] == vertex)
        {
            gradient.set_VE(old_edge[0], old_edge[1], old_ef, mesh);
            Forman_Gradient_Topological_Relations::set_VTstar(n, old_edge[0], old_ef.first, local_rels.get_VTstars(), cache.get_vtstar_cache(), root, division);
        }
        else
        {
            gradient.set_VE(old_edge[1], old_edge[0], old_ef, mesh);
            Forman_Gradient_Topological_Relations::set_VTstar(n, old_edge[1], old_ef.first, local_rels.get_VTstars(), cache.get_vtstar_cache(), root, division);
        }
        /// save for the next step
        old_edge = edge;
        old_ef = ef;
    }
    // assert(next_vertex == ex_minimum);
    gradient.set_VE(next_vertex, vertex, old_ef, mesh);
    Forman_Gradient_Topological_Relations::set_VTstar(n, next_vertex, old_ef.first, local_rels.get_VTstars(), cache.get_vtstar_cache(), root, division);
}

void Forman_Gradient_Simplifier::removal_update_gradient(itype triangle, iNode *saddle, Mesh &mesh, Forman_Gradient &gradient, leaf_ET &local_ef, mig_cache &cache,
                                                         Node_V &n, Node_V &root, Spatial_Subdivision &division)
{
    ivect old_edge;
    ivect edge;
    pair<itype, itype> twotriangles = saddle->get_edge_id();
    Triangle &tri = mesh.get_triangle(twotriangles.first);

    if (twotriangles.second >= 0)
    {
        ivect e;
        Triangle &tri2 = mesh.get_triangle(twotriangles.second);
        for (int j = 0; j < tri2.vertices_num(); j++)
        {
            tri.TE(j, e);
            if (tri2.has_edge(e))
            {
                old_edge = e;
                break;
            }
        }
    }
    else
    {
        tri.TE(-twotriangles.second - 1, old_edge);
    }

    /// ----------------- GRADIENTE UPDATES ----------------- ///
    itype next_triangle = triangle;

    while (!gradient.is_triangle_critical(next_triangle))
    {
        Triangle &tri = mesh.get_triangle(next_triangle);
        int i = get_paired_edge(gradient, tri, next_triangle, edge);

        triangle = next_triangle;

        pair<itype, itype> et = Forman_Gradient_Topological_Relations::get_ET(n, edge, local_ef, cache.get_et_cache(), root, division, mesh);

        next_triangle = (et.first == triangle) ? et.second : et.first;
        int vertex = mesh.get_triangle(triangle).TV(i);

        gradient.free_ET(i, triangle);

        i = tri.edge_index(old_edge);
        vertex = mesh.get_triangle(triangle).TV(i);
        gradient.set_ET(i, triangle);

        old_edge = edge;
    }

    Triangle &tri_next = mesh.get_triangle(next_triangle);
    int i = tri_next.edge_index(old_edge);

    gradient.set_ET(i, next_triangle);
}

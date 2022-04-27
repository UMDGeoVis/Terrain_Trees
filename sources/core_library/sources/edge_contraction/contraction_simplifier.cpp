/*
    This file is part of the Stellar library.

    Author(s): Riccardo Fellegara (riccardo.fellegara@gmail.com)

    This project has been supported by the Italian Ministry of Education and
    Research under the PRIN 2009 program, and by the National Science Foundation
    under grant number IIS-1116747.

    The Stellar library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The Stellar library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the Stellar library.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "contraction_simplifier.h"


void Contraction_Simplifier::find_candidate_edges(Node_V &n, Mesh &mesh, leaf_VT &vts, edge_queue &edges, contraction_parameters &params)
{
    map<ivect, coord_type> lengths;
    ivect e;

    for (RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const &t_id = itPair.first;

        if (mesh.is_triangle_removed(*t_id))
        {
            //  cout<<"triangle removed"<<endl;
            continue;
        }
        Triangle &t = mesh.get_triangle(*t_id);

        for (int i = 0; i < t.vertices_num(); i++)
        {
            t.TE(i, e);

            if (n.indexes_vertex(e[1])) // e (v1,v2) is a candidate edge if at least v2 is in n
            {
                map<ivect, coord_type>::iterator it = lengths.find(e);
                // cout<<e[0]<<" and "<<e[1]<<endl;
                if (it == lengths.end())
                {

                    coord_type length;
                    Vertex &v1 = mesh.get_vertex(e[0]);
                    Vertex &v2 = mesh.get_vertex(e[1]);
                    dvect dif = {v1.get_x() - v2.get_x(), v1.get_y() - v2.get_y(), v1.get_z() - v2.get_z()};
                    //  cout<<dif[0]<<", "<<dif[1]<<", "<<dif[2]<<endl;
                    length = sqrt(dif[0] * dif[0] + dif[1] * dif[1] + dif[2] * dif[2]);

                    //  Edge e((*it)[0],(*it)[1]);
                    lengths[e] = length;
                    //Edge edge_obj(e[0],e[1]);
                    if (params.no_limit())
                    {
                        edges.push(new Geom_Edge(e, length));
                    }
                    else
                    {
                        if (length - params.get_maximum_limit() < Zero)
                        {
                            // cout<<"["<<e[0]<<","<<e[1]<<"]  Edge length: "<<length<<endl;
                            edges.push(new Geom_Edge(e, length));
                            //     cout<<"ENQUEUE"<<endl;
                        }
                    }
                }
            }
        }
    }
    //  cout<<"======NEXT NODE======"<<endl;
}

void Contraction_Simplifier::find_candidate_edges_parallel(Node_V &n, Mesh &mesh, leaf_VT &vts, edge_queue &edges, contraction_parameters &params, bool is_cross)
{
    map<ivect, coord_type> lengths;
    ivect e;

    for (RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const &t_id = itPair.first;

        if (mesh.is_triangle_removed(*t_id))
        {
            //  cout<<"triangle removed"<<endl;
            continue;
        }
        Triangle &t = mesh.get_triangle(*t_id);

        for (int i = 0; i < t.vertices_num(); i++)
        {
            t.TE(i, e);

            if (n.indexes_vertex(e[1])) // e (v1,v2) is a candidate edge if at least v2 is in n
            {
                map<ivect, coord_type>::iterator it = lengths.find(e);
                // cout<<e[0]<<" and "<<e[1]<<endl;
                if (it == lengths.end())
                {

                    coord_type length;
                    Vertex &v1 = mesh.get_vertex(e[0]);
                    Vertex &v2 = mesh.get_vertex(e[1]);
                    dvect dif = {v1.get_x() - v2.get_x(), v1.get_y() - v2.get_y(), v1.get_z() - v2.get_z()};
                    //  cout<<dif[0]<<", "<<dif[1]<<", "<<dif[2]<<endl;
                    length = sqrt(dif[0] * dif[0] + dif[1] * dif[1] + dif[2] * dif[2]);

                    //  Edge e((*it)[0],(*it)[1]);
                    lengths[e] = length;
                    //Edge edge_obj(e[0],e[1]);
                    if (params.no_limit())
                    {
                        edges.push(new Geom_Edge(e, length));
                    }
                    else
                    {
                        if (length - params.get_maximum_limit() < Zero)
                        {
                            // cout<<"["<<e[0]<<","<<e[1]<<"]  Edge length: "<<length<<endl;
                            edges.push(new Geom_Edge(e, length));
                            //     cout<<"ENQUEUE"<<endl;
                        }
                    }
                }
            }
        }
    }
    //}

    //  cout<<"= =====NEXT NODE======"<<endl;
}


void Contraction_Simplifier::contract_edge(ivect &e, ET &et, VT &vt0, VT &vt1, Node_V &outer_v_block, edge_queue &edges,
                                           Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params, map<vector<int>, double>& updated_edges)
{
    //cout<<"[EDGE CONTRACTION] v1 and v2:"<<e[0]<<", "<<e[1]<<endl;
    // cout<<"[NOTICE] Contract Edge"<<endl;

    ivect et_vec;
    et_vec.push_back(et.first);
    if (et.second != -1)
        et_vec.push_back(et.second);

    difference_of_vectors(vt0, et_vec); // vt0 now contains the difference VT0 - ET
    difference_of_vectors(vt1, et_vec); // vt1 now contains the difference VT1 - ET

    // cout<<"VT1 size: "<<vt0.size()<<" VT2 size: "<<vt1.size()<<endl;
    // contract v1 to v0.

    /// prior checking the d-1 faces we update
    /// (1) the corresponding vt0 relation (by adding the triangles in vt1-et to vt0)
    /// (2) then the outer_v_block (if e is a cross edge)
    /// (3) and, finally, we add the new edges crossing b to the edge queue

    Contraction_Simplifier::update(e, vt0, vt1, n, outer_v_block, edges, mesh, params,updated_edges);

    // we remove v2 and the triangles in et
    Contraction_Simplifier::remove_from_mesh(e[1], et, mesh, params);
    // finally we clear the VT(v2)
    vt1.clear();
    //et.clear();
}

void Contraction_Simplifier::contract_edge(ivect &e, ET &et, VT &vt0, VT &vt1, Node_V &outer_v_block, edge_queue &edges,
                                           Node_V &n, Mesh &mesh, contraction_parameters &params, map<vector<int>, double> &updated_edges)
{
    //cout << "[EDGE CONTRACTION] v1 and v2:" << e[0] - 1 << ", " << e[1] - 1 << endl;
    // cout<<"[NOTICE] Contract Edge"<<endl;
    // omp_set_lock(&(v_locks[e[0] - 1]));
    // omp_set_lock(&(v_locks[e[1] - 1]));
    ivect et_vec;
    et_vec.push_back(et.first);
    if (et.second != -1)
        et_vec.push_back(et.second);

    difference_of_vectors(vt0, et_vec); // vt0 now contains the difference VT0 - ET
    difference_of_vectors(vt1, et_vec); // vt1 now contains the difference VT1 - ET

    // cout<<"VT1 size: "<<vt0.size()<<" VT2 size: "<<vt1.size()<<endl;
    // contract v1 to v0.

    /// prior checking the d-1 faces we update
    /// (1) the corresponding vt0 relation (by adding the triangles in vt1-et to vt0)
    /// (2) then the outer_v_block (if e is a cross edge)
    /// (3) and, finally, we add the new edges crossing b to the edge queue


    Contraction_Simplifier::update_parallel(e, vt0, vt1, n, outer_v_block, edges, mesh, params, updated_edges);


    // we remove v2 and the triangles in et
    Contraction_Simplifier::remove_from_mesh(e[1], et, mesh, params);
    // finally we clear the VT(v2)
    vt1.clear();
    //et.clear();
    // omp_unset_lock(&(v_locks[e[0] - 1]));
    // omp_unset_lock(&(v_locks[e[1] - 1]));
}



void Contraction_Simplifier::get_ET(ivect &e, ET &et, Node_V &n, Mesh &mesh, leaf_VT &vts)
{
    int other_v, local_v_id;
    if (n.indexes_vertex(e[0]))
    {
        other_v = e[1];
        local_v_id = e[0] - n.get_v_start();
    }
    else
    {
        other_v = e[0];
        local_v_id = e[1] - n.get_v_start();
    }

    VT &vt = vts[local_v_id];

    ivect et_tmp;
    // cout<<"vt size:"<<vt.size()<<endl;
    //  #pragma omp critical
    {
        for (unsigned i = 0; i < vt.size(); i++)
        {

            Triangle &t = mesh.get_triangle(vt[i]);
            if (t.has_vertex(other_v) && (!mesh.is_triangle_removed(vt[i])))
            {
                et_tmp.push_back(vt[i]);

 
            }
        }
    }

    if (et_tmp.size() == 2)
        et = make_pair(et_tmp[0], et_tmp[1]);
    else if (et_tmp.size() == 0) // Can be deleted, since the default is (-1,-1)
    {
        et = make_pair(-1, -1);
    }
    else
    {
        et = make_pair(et_tmp[0], -1);
    }
}

void Contraction_Simplifier::clean_coboundary(VT &cob, Mesh &mesh)
{

    for (ivect_iter it = cob.begin(); it != cob.end();)
    {
        int tid = *it;
        if (mesh.is_triangle_removed(*it))
            cob.erase(it);
        else
            ++it;
    }
}

void Contraction_Simplifier::update(const ivect &current_edge, VT &vt, VT &difference, Node_V &n, Node_V &v_block,
                                    edge_queue &edges, Mesh &mesh, contraction_parameters &params, map<vector<int>, double>& updated_edges)
{
    set<ivect> e_set; /// we insert the new edges first in this set to avoid duplicate insertions in the queue

    if (params.is_QEM())
    {
        initialQuadric[current_edge[0]] += initialQuadric[current_edge[1]]; // Update the QEM of e[0]
        set<ivect> v1_related_set;
        for (int i = 0; i < vt.size(); i++)
        { //Store VE(v1)
            ivect new_e(2, 0);
            Triangle &t = mesh.get_triangle(vt[i]);
            int v1_pos = t.vertex_index(current_edge[0]);
            for (int j = 0; j < t.vertices_num(); j++)
            {
                if (j != v1_pos)
                {
                    t.TE(j, new_e);
                    //     if(n.indexes_vertex(new_e[1])) /// we process an edge only if it has all the extrema already processed
                    v1_related_set.insert(new_e);
                }
            }
        }
        for (auto it = v1_related_set.begin(); it != v1_related_set.end(); ++it)
        {
            //Calculate updated edge costs of VE(v1)
            ivect new_edge(2, -1);
            Vertex &v1 = mesh.get_vertex((*it)[0]);
            Vertex &v2 = mesh.get_vertex((*it)[1]);
            int new_vertex_pos = -1;
            double error = compute_error((*it)[0], (*it)[1], mesh, new_vertex_pos);
            //    cout<<"[DEBUG] calculated error: "<<error<<endl;
            assert(new_vertex_pos != -1);

            if (updated_edges.find(*it) != updated_edges.end())
                updated_edges[*it] = error;
            else
            {
                pair<ivect, double> updated_edge(*it, error);
                updated_edges.insert(updated_edge);
            }
            if (new_vertex_pos == 1)
                {
                    new_edge = {(*it)[1], (*it)[0]};
                }
            else
                new_edge = {(*it)[0], (*it)[1]};
            bool error_condition;
            if(params.no_limit())
            error_condition = true;
            else
            error_condition =(error + Zero <params.get_maximum_limit());
            if ( error_condition&& n.indexes_vertex(new_edge[1]))
            {
                  //cout<<"["<<e[0]-1<<","<<e[1]-1<<"]  Error will be introduced: "<<error<<endl;
                edges.push(new Geom_Edge(new_edge, error));
            }
        }
    }

        for (ivect_iter it = difference.begin(); it != difference.end(); ++it)
        {
            
            int pos = -1;
            //   bool not_valid=false;
            Triangle &t = mesh.get_triangle(*it);

            /// then we update the triangle changing e[1] with e[0]
            pos = t.vertex_index(current_edge[1]);
            t.setTV_keep_border(pos, current_edge[0]);

            /// before updating the triangle, we check
            /// if the leaf block indexing e[0] does not contain the current triangle we have to add it
            /// NOTA: there is one possible case.. as leaf block n already indexes the triangle in e[1]
            /// NOTA2: we have just to check that n and v_block are different (as if they are equal the edge is internal in n)
            // Note that in the triangle *it, the e[1] has been replaced by e[0] already

            if (n.get_v_start() != v_block.get_v_start())
            {
                 if (/*current_edge[0] < current_edge[1] &&*/ !v_block.index_triangle(*it, mesh, current_edge[0]))
                 {
                    v_block.add_triangle(*it);
                 }
                // else if (current_edge[0] > current_edge[1] && !n.index_triangle(*it, mesh, current_edge[0]))
                // {
                //     n.add_triangle(*it);
                // }
            }

           
            // if (params.is_QEM())
            // {
                ivect new_e(2, 0); // new_e.assign(2,0);
                for (int i = 0; i < 3; i++)
                {
                    if (i != pos)
                    {
                        t.TE(i, new_e); //t.TE(new_e,pos,i); need to check
                        e_set.insert(new_e);
                    }
                }
            // }
            // else
            // {
            //     /// we have to add the new edges in the queue
            //     ivect new_e;
            //     new_e.assign(2, 0);

            //         for (int i = 0; i < 3; i++)
            //         {
            //             if (i != pos)
            //             {
            //                 mesh.get_triangle(*it).TE(i, new_e); //t.TE(new_e,pos,i); need to check
            //                 e_set.insert(new_e);
            //             }
            //         }
                
            // }
        }
    

    /// we push the new "unique" edges in the queue
    for (auto it = e_set.begin(); it != e_set.end(); ++it)
    {
        //Calculate length
        
        ivect new_edge(2, -1);

        if (params.is_QEM())
        {
            int new_vertex_pos = -1;
            double error = compute_error((*it)[0], (*it)[1], mesh, new_vertex_pos);
            //    cout<<"[DEBUG] calculated error: "<<error<<endl;
            assert(new_vertex_pos != -1);

            if (updated_edges.find(*it) != updated_edges.end()) // updated_edges keep all the edges that have new cost values
                updated_edges[*it] = error;
            else
            {
                pair<ivect, double> updated_edge(*it, error);
                updated_edges.insert(updated_edge);
            }
            if (new_vertex_pos == 1)
                {
                    new_edge = {(*it)[1], (*it)[0]};
                }
            else
                 new_edge = {(*it)[0], (*it)[1]};
            bool error_condition=false;
            if(params.no_limit())
            error_condition = true;
            else
            error_condition = (error + Zero <params.get_maximum_limit() );

            if (error_condition && n.indexes_vertex(new_edge[1]))
            {

          //       cout<<"["<<e[0]-1<<","<<e[1]-1<<"]  Error will be introduced: "<<error<<endl;

                edges.push(new Geom_Edge(new_edge, error));
            }
        }
        else
        {
                double value;
                Vertex &v1 = mesh.get_vertex((*it)[0]);
                Vertex &v2 = mesh.get_vertex((*it)[1]);
                dvect dif = {v1.get_x() - v2.get_x(), v1.get_y() - v2.get_y(), v1.get_z() - v2.get_z()};
                value = sqrt(dif[0] * dif[0] + dif[1] * dif[1] + dif[2] * dif[2]);
                new_edge = {(*it)[0], (*it)[1]};

                 bool error_condition;
            if(params.no_limit())
            error_condition = true;
            else
            error_condition = (value + Zero < params.get_maximum_limit()) ;

                if (error_condition&& n.indexes_vertex(new_edge[1]))
                {
                    // Geom_Edge new_edge(e,length);
                   // cout << "[" << e[0] << "," << e[1] << "]  Edge length: " << value << endl;
                    edges.push(new Geom_Edge(new_edge, value));
                }
            
        }
    }

    /// finally we update the VT relation of e[0]
    unify_vectors(vt, difference);
}

void Contraction_Simplifier::remove_from_mesh(int to_delete_v, ET &et, Mesh &mesh, contraction_parameters &params)
{
    if (et.first != -1)
    {
        //omp_set_lock(&(t_locks[et.first - 1]));
        mesh.remove_triangle(et.first);
        //  omp_unset_lock(&(t_locks[et.first - 1]));
        //params.increment_counter();
    }
    if (et.second != -1)
    {
        // omp_set_lock(&(t_locks[et.second - 1]));
        mesh.remove_triangle(et.second);
        //  omp_unset_lock(&(t_locks[et.second - 1]));
        //  params.increment_counter();
    }

    mesh.remove_vertex(to_delete_v);

    params.increment_contracted_edges_counter();
}

bool Contraction_Simplifier::link_condition(int v0, int v1, VT &vt0, VT &vt1, ET &et, Mesh &mesh)
{
    //iset link_ab;
    //  iset link_e;
    //iset vts;

    //Update: Considering that the edge to be contracted should not be boundary edge
    //We can simplify the link condition check while checking if e is boundary edge
    if (et.first == -1 || et.second == -1)
        return false;
    int counter = 0;
    {
        iset vv0, vv1;
      
        for (int i = 0; i < vt0.size(); i++)
        {
            if (!mesh.is_triangle_removed(vt0[i]))
            {
                Triangle t = mesh.get_triangle(vt0[i]);
                int v0_id = t.vertex_index(v0);
                int v01_id = t.TV((v0_id + 1) % 3);
                int v02_id = t.TV((v0_id + 2) % 3);
                vv0.insert(v01_id);
                vv0.insert(v02_id);
            }
        }

        for (int i = 0; i < vt1.size(); i++)
        {
            if (!mesh.is_triangle_removed(vt1[i]))
            {
                Triangle t = mesh.get_triangle(vt1[i]);
                int v1_id = t.vertex_index(v1);

                vv1.insert(t.TV((v1_id + 1) % 3));
                vv1.insert(t.TV((v1_id + 2) % 3));
            }
        }
        //  cout<<v1<<"'s VV size: "<<vv1.size()<<endl;
        for (iset_iter it = vv1.begin(); it != vv1.end(); it++)
        {
            if (vv0.find(*it) != vv0.end())
            {
                //   link_ab.insert(*it);
                counter++;
            }
        }

    }
  
    return counter <= 2;
}

//// Function below is no longer used 
// bool Contraction_Simplifier::link_condition(int v0, int v1, VT &vt0, VT &vt1, ET &et, Node_V &n, VV &vv_locks, Mesh &mesh)
// {

//     //Update: Considering that the edge to be contracted should not be boundary edge
//     //We can simplify the link condition check while checking if e is boundary edge
//     if (et.first == -1 || et.second == -1)
//         return false;
//     int counter = 0;
//     //#pragma omp critical
//     {
//         iset vv0, vv1;
//         //vts.insert(vts.end(),vt0.begin(),vt0.end());
//         VT vts = vt0;
//         unify_vectors(vts, vt1);
//         // for (int i = 0; i < vts.size(); i++)
//         // {
//         //     omp_set_lock(&(t_locks[vts[i] - 1]));
//         // }

//         for (int i = 0; i < vt0.size(); i++)
//         {
//             if (!mesh.is_triangle_removed(vt0[i]))
//             {
//                 Triangle t = mesh.get_triangle(vt0[i]);
//                 int v0_id = t.vertex_index(v0);
//                 int v01_id = t.TV((v0_id + 1) % 3);
//                 int v02_id = t.TV((v0_id + 2) % 3);
//                 if (!n.indexes_vertex(v01_id))
//                     vv_locks.insert(v01_id);

//                 if (!n.indexes_vertex(v02_id))
//                     vv_locks.insert(v02_id);

//                 vv0.insert(v01_id);
//                 vv0.insert(v02_id);
//             }
//         }

//         for (int i = 0; i < vt1.size(); i++)
//         {
//             if (!mesh.is_triangle_removed(vt1[i]))
//             {
//                 Triangle t = mesh.get_triangle(vt1[i]);
//                 int v1_id = t.vertex_index(v1);
//                 int v11_id = t.TV((v1_id + 1) % 3);
//                 int v12_id = t.TV((v1_id + 2) % 3);
//                 if (!n.indexes_vertex(v11_id))
//                     vv_locks.insert(v11_id);
//                 if (!n.indexes_vertex(v12_id))
//                     vv_locks.insert(v12_id);
//                 vv1.insert(v11_id);
//                 vv1.insert(v12_id);
//             }
//         }

//         //  cout<<v1<<"'s VV size: "<<vv1.size()<<endl;
//         for (iset_iter it = vv1.begin(); it != vv1.end(); it++)
//         {
//             if (vv0.find(*it) != vv0.end())
//             {
//                 //   link_ab.insert(*it);
//                 counter++;
//             }
//         }

//         // for (int i = 0; i < vts.size(); i++)
//         // {
//         //     omp_unset_lock(&(t_locks[vts[i] - 1]));
//         // }

//         // for (iset_iter it = vv_locks.begin(); it != vv_locks.end(); it++)
//         // {
//         //     omp_set_lock(&(v_locks[*it - 1]));
//         // }

//         // for (int i=0;i<vt1.size();i++){
//         //     omp_unset_lock(&(t_locks[vt1[i]]));
//         // }
//         // if(et.first!=-1){
//         //     if(!mesh.is_triangle_removed(et.first)){
//         //     Triangle t1= mesh.get_triangle(et.first);
//         //     for(int i=0;i<3;i++){
//         //         if(t1.TV(i)!=v0&&t1.TV(i)!=v1)
//         //          link_e.insert(t1.TV(i));
//         //     }
//         //     }
//         // }
//         // if(et.second!=-1){
//         //      if(!mesh.is_triangle_removed(et.second)){
//         //     Triangle t2= mesh.get_triangle(et.second);
//         //     for(int i=0;i<3;i++){
//         //         if(t2.TV(i)!=v0&&t2.TV(i)!=v1)
//         //          link_e.insert(t2.TV(i));
//         //     }
//         // }
//         // }
//     }
//     return counter <= 2;
// }

bool Contraction_Simplifier::link_condition(int v0, int v1, VT &vt0, VT &vt1, ET &et, Node_V &n, Node_V &v_block, VV &vv_locks, Mesh &mesh)
{

    //Update: Considering that the edge to be contracted should not be boundary edge
    //We can simplify the link condition check while checking if e is boundary edge
    if (et.first == -1 || et.second == -1)
        return false;
    int counter = 0;
  

        iset vv0, vv1;
        //vts.insert(vts.end(),vt0.begin(),vt0.end());
        VT vts = vt0;

        unify_vectors(vts, vt1);
        // for (int i = 0; i < vts.size(); i++)
        // {
        //     omp_set_lock(&(t_locks[vts[i] - 1]));
        // }

        for (int i = 0; i < vt0.size(); i++)
        {
            // if (!mesh.is_triangle_removed(vt0[i]))
            // {
                Triangle t = mesh.get_triangle(vt0[i]);
                int v0_id = t.vertex_index(v0);
                int v01_id = t.TV((v0_id + 1) % 3);
                int v02_id = t.TV((v0_id + 2) % 3);

                // if (n.get_v_start() != v_block.get_v_start()&&!n.indexes_vertex(v01_id) && !v_block.indexes_vertex(v01_id))
                //     vv_locks.insert(v01_id);

                // if (n.get_v_start() != v_block.get_v_start()&&!n.indexes_vertex(v02_id) && !v_block.indexes_vertex(v02_id))
                //     vv_locks.insert(v02_id);

                vv0.insert(v01_id);
                vv0.insert(v02_id);
            // }
        }

        for (int i = 0; i < vt1.size(); i++)
        {
            // if (!mesh.is_triangle_removed(vt1[i]))
            // {
                Triangle t = mesh.get_triangle(vt1[i]);
                int v1_id = t.vertex_index(v1);
                int v11_id = t.TV((v1_id + 1) % 3);
                int v12_id = t.TV((v1_id + 2) % 3);
                if (n.get_v_start() != v_block.get_v_start()&&!n.indexes_vertex(v11_id) && !v_block.indexes_vertex(v11_id))
                    vv_locks.insert(v11_id);
                if (n.get_v_start() != v_block.get_v_start()&&!n.indexes_vertex(v12_id) && !v_block.indexes_vertex(v12_id))
                    vv_locks.insert(v12_id);
                vv1.insert(v11_id);
                vv1.insert(v12_id);
            // }
        }

        //  cout<<v1<<"'s VV size: "<<vv1.size()<<endl;
        for (iset_iter it = vv1.begin(); it != vv1.end(); it++)
        {
            if (vv0.find(*it) != vv0.end())
            {
                //   link_ab.insert(*it);
                counter++;
            }
        }

        // for (iset_iter it = vv_locks.begin(); it != vv_locks.end(); it++)
        // {
        //     omp_set_lock(&(v_locks[*it - 1]));
        // }
    


    return counter <= 2;
}

void Contraction_Simplifier::update_parallel(const ivect &current_e, VT &vt, VT &difference, Node_V &n, Node_V &v_block, edge_queue &edges,
                                             Mesh &mesh, contraction_parameters &params, map<vector<int>, double> &updated_edges)
{
    set<ivect> e_set; /// we insert the new edges first in this set to avoid duplicate insertions in the queue

    /// before updating the triangle, we check
    /// if the leaf block indexing e[0] does not contain the current triangle we have to add it
    /// NOTA: there is one possible case.. as leaf block n already indexes the triangle in e[1]
    /// NOTA2: we have just to check that n and v_block are different (as if they are equal the edge is internal in n)
    //cout<<"Update the relations"<<endl;
    if (params.is_QEM())
    {
        initialQuadric[current_e[0]] += initialQuadric[current_e[1]]; // Update the QEM of e[0]
        set<ivect> v1_related_set;
        for (int i = 0; i < vt.size(); i++)
        { //Store VE(v1)
            ivect new_e(2, 0);
            Triangle &t = mesh.get_triangle(vt[i]);
            int v1_pos = t.vertex_index(current_e[0]);
            for (int j = 0; j < t.vertices_num(); j++)
            {
                if (j != v1_pos)
                {
                    t.TE(j, new_e);
                    //     if(n.indexes_vertex(new_e[1])) /// we process an edge only if it has all the extrema already processed
                    v1_related_set.insert(new_e);
                }
            }
        }
        for (auto it = v1_related_set.begin(); it != v1_related_set.end(); ++it)
        {
            //Calculate updated edge costs of VE(v1)
            ivect new_edge(2, -1);
            Vertex &v1 = mesh.get_vertex((*it)[0]);
            Vertex &v2 = mesh.get_vertex((*it)[1]);
            int new_vertex_pos = -1;
            double error = compute_error((*it)[0], (*it)[1], mesh, new_vertex_pos);
            //   cout<<"[DEBUG] calculated error: "<<error<<endl;
            //assert(new_vertex_pos != -1);

            if (updated_edges.find(*it) != updated_edges.end())
                updated_edges[*it] = error;
            else
            {
                pair<ivect, double> updated_edge(*it, error);
                // #pragma omp critical
                updated_edges.insert(updated_edge);
            }

            if(params.output_stats()){
                #pragma omp critical 
                {
                    edge_costs_output[make_pair((*it)[0], (*it)[1])] = error;
            
                }
            }

            if (new_vertex_pos == 1)
                {
                    new_edge = {(*it)[1], (*it)[0]};
                }
            else
             new_edge = {(*it)[0], (*it)[1]};
            bool error_condition;
            if (params.no_limit())
                error_condition = true;
            else
                error_condition = (error + Zero < params.get_maximum_limit());


            if (error_condition && n.indexes_vertex(new_edge[1]))
            {
             // cout<<"["<<e[0]-1<<","<<e[1]-1<<"]  Error will be introduced: "<<error<<endl;
                edges.push(new Geom_Edge(new_edge, error));
            }
        }
    }

    for (ivect_iter it = difference.begin(); it != difference.end(); ++it)
    {
        int pos = -1;
        /// we have to add the new edges in the queue
 

        Triangle &t = mesh.get_triangle(*it);
        if (n.get_v_start() != v_block.get_v_start())
        {
            if (/*e[0] < e[1] &&*/ !v_block.index_triangle(*it, mesh, current_e[1]))
            {
                v_block.add_triangle(*it);
            }
            // else if (e[0] > e[1] && !n.index_triangle(*it, mesh, e[1]))
            // {
            //     n.add_triangle(*it);
            // }
        }

            /// then we update the triangle changing e[1] with e[0]
            pos = t.vertex_index(current_e[1]);
            t.setTV_keep_border(pos, current_e[0]);
            ivect new_e;
            new_e.assign(2, 0);
            for (int i = 0; i < 3; i++)
            {
                if (i != pos)
                {
                    t.TE(i, new_e); //t.TE(new_e,pos,i); need to check
                    e_set.insert(new_e);
                }
            }

    }

    /// we push the new "unique" edges in the queue
    for (auto it = e_set.begin(); it != e_set.end(); ++it)
    {
   
        ivect new_edge(2, -1);
        // if (n.indexes_vertex(e[1]))
        if (params.is_QEM())
        {
            int new_vertex_pos = -1;
            double error = compute_error((*it)[0], (*it)[1], mesh, new_vertex_pos);
            //    cout<<"[DEBUG] calculated error: "<<error<<endl;
          //  assert(new_vertex_pos != -1);
            if (updated_edges.find(*it) != updated_edges.end()) // updated_edges keep all the edges that have new cost values
                updated_edges[*it] = error;
            else
            {
                pair<ivect, double> updated_edge(*it, error);
                // #pragma omp critical
                updated_edges.insert(updated_edge);
            }

            if(params.output_stats()){
                #pragma omp critical 
                {
                    edge_costs_output[make_pair((*it)[0], (*it)[1])] = error;
            
                }
            }


           if (new_vertex_pos == 1)
             {
               new_edge = {(*it)[1], (*it)[0]};
            }
            else
            new_edge = {(*it)[0], (*it)[1]};
            bool error_condition;
            if (params.no_limit())
                error_condition = true;
            else
                error_condition = (error + Zero <params.get_maximum_limit() );

            if (error_condition && n.indexes_vertex(new_edge[1]))
            {
                 // cout<<"["<<e[0]-1<<","<<e[1]-1<<"]  Error will be introduced: "<<error<<endl;
                edges.push(new Geom_Edge(new_edge, error));
            }
        }
        else
        {
            double value;
            Vertex &v1 = mesh.get_vertex((*it)[0]);
            Vertex &v2 = mesh.get_vertex((*it)[1]);
            dvect dif = {v1.get_x() - v2.get_x(), v1.get_y() - v2.get_y(), v1.get_z() - v2.get_z()};
            value = sqrt(dif[0] * dif[0] + dif[1] * dif[1] + dif[2] * dif[2]);
            new_edge = {(*it)[0], (*it)[1]};
            bool error_condition;
            if (params.no_limit())
                error_condition = true;
            else
                error_condition = (value + Zero < params.get_maximum_limit());

            if (error_condition && n.indexes_vertex(new_edge[1]))
            {
                //     cout<<"["<<e[0]<<","<<e[1]<<"]  Edge length: "<<value<<endl;
                edges.push(new Geom_Edge(new_edge, value));
            }
        }
    }
    // for(ivect_iter it=t_locks_id.begin();it!=t_locks_id.end();it++){
    //     omp_unset_lock(&(t_locks[*it]));
    // }
    /// finally we update the VT relation of e[0]
    unify_vectors(vt, difference);
}

void Contraction_Simplifier::simplify_parallel(PRT_Tree &tree, Mesh &mesh, cli_parameters &cli)
{
    // cerr << "[NOTICED] Cache size: " << cli.cache_size << endl;
    // LRU_Cache<int, leaf_VT> cache(cli.cache_size); // the key is v_start while the value are the VT relations
    contraction_parameters params;
    if (cli.contract_all_edges == true)
        params.contract_all_possible_edges();
    params.set_maximum_limit(cli.maximum_limit);
   // omp_set_num_threads(cli.num_of_threads);
    // Set edge selection criteria
    if (cli.QEM_based)
        params.queue_criterion_QEM();
    else
    {
        params.queue_criterion_length();
    }
    // Set to be parallel mode
    params.parallel_compute();
    Timer time;
    int simplification_round;
    int round = 1;
    time.start();
    cout << "Number of threads used in the simplification:" << omp_get_max_threads() << endl;
    //  const int t_num = mesh.get_triangles_num();
   // const int v_num = mesh.get_vertices_num();
    const int l_num = tree.get_leaves_number();

    //  t_locks.resize(t_num);
  //  v_locks.resize(v_num);
    l_locks.resize(l_num);

// #pragma omp parallel for
//     for (int i = 0; i < t_num; i++)
//         omp_init_lock(&(t_locks[i]));
//     cout << "Initialize t_locks" << endl;
// #pragma omp parallel for
//     for (int i = 0; i < v_num; i++)
//         omp_init_lock(&(v_locks[i]));
//     cout << "Initialize v_locks" << endl;

#pragma omp parallel for
    for (int i = 0; i < l_num; i++)
        omp_init_lock(&(l_locks[i]));
    cout << "Initialize l_locks" << endl;

    time.stop();
    time.print_elapsed_time("[TIME] Initialization of locks: ");

    if (params.is_QEM())
    {
        time.start();
      //  trianglePlane = vector<dvect>(mesh.get_triangles_num(), dvect(4, 0));
        initialQuadric = vector<Matrix>(mesh.get_vertices_num() + 1, Matrix(0.0));
        cout << "=========Calculate triangle plane========" << endl;
   //     compute_triangle_plane(mesh, trianglePlane);
        cout << "=========Calculate initial QEM========" << endl;
    //    compute_initial_QEM_parallel(tree, mesh, trianglePlane);
        compute_initial_plane_and_QEM_parallel(tree, mesh);

        time.stop();
        time.print_elapsed_time("[TIME] Calculating initial QEM: ");
    }

    time.start();
    while (1)
    {
        simplification_round = params.get_contracted_edges_num(); //checked edges
        /// HERE YOU NEED TO DEFINE A PROCEDURE FOR SIMPLIFY THE TIN BY USING THE SPATIAL INDEX
        this->simplify_compute_parallel(mesh, tree.get_subdivision(), params, tree);

        cout << "Num of edges enqueued:" << params.get_sum_edge_queue_sizes() << endl;
        // PARTIAL SIMPLIFICATION STATS
        cerr << "=== end-of-round " << round << ") --> contracted edges: ";
        cerr << params.get_contracted_edges_num() - simplification_round << endl;
        round++;

        if (cli.debug_mode)
        {
            time.stop();
            time.print_elapsed_time("   [TIME] executing a simplification round: ");
            time.start();
            //  cerr << "   [RAM] peak for executing a simplification round: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
        }
        if (simplification_round == params.get_contracted_edges_num())
            break;
        
        Contraction_Simplifier::preprocess(tree, mesh, cli);

    }

    // #pragma omp parallel for
    //     for (int i = 0; i < t_num; i++)
    //         omp_destroy_lock(&(t_locks[i]));

// #pragma omp parallel for
//     for (int i = 0; i < v_num; i++)
//         omp_destroy_lock(&(v_locks[i]));

#pragma omp parallel for
    for (int i = 0; i < l_num; i++)
        omp_destroy_lock(&(l_locks[i]));

    time.stop();

    ///// Clear all the auxiliary data structures.
    //v_locks.clear();
    //vector<omp_lock_t>().swap(v_locks);
    vector<omp_lock_t>().swap(l_locks);
    vector<int>().swap(v_in_leaf);
    lists_leafs().swap(conflict_leafs);
    // l_locks.clear();
    // vector<Matrix>().swap(initialQuadric);
    // vector<dvect>().swap(trianglePlane);
    // map<vector<int>, double>().swap(updated_edges);
    if (!cli.debug_mode)
        time.print_elapsed_time("[TIME] Edge contraction simplification: ");
    cerr << "[MEMORY] peak for Simplification: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
    //else
    //params.print_simplification_partial_timings();
    //params.print_simplification_counters();
    /// finally we have to update/compress the mesh and the tree
    time.start();
    Contraction_Simplifier::update_mesh_and_tree(tree, mesh, params, cli);
    time.stop();
    time.print_elapsed_time("[TIME] Mesh and tree updating: ");

    cerr << "[MEMORY] peak for mesh and tree updating: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
}

void Contraction_Simplifier::simplify(PRT_Tree &tree, Mesh &mesh, cli_parameters &cli)
{

    //cerr<<"==Homology preserving simplification - weak-link condition=="<<endl;

    cerr << "[NOTICED] Cache size: " << cli.cache_size << endl;
    LRU_Cache<int, leaf_VT> cache(cli.cache_size); // the key is v_start while the value are the VT relations
    contraction_parameters params;
    if (cli.contract_all_edges == true)
        params.contract_all_possible_edges();
    params.set_maximum_limit(cli.maximum_limit);
    if (cli.QEM_based)
        params.queue_criterion_QEM();
    else
    {
        params.queue_criterion_length();
    }

    Timer time;
    int simplification_round;
    int round = 1;
    if (params.is_QEM())
    {

       // trianglePlane = vector<dvect>(mesh.get_triangles_num(), dvect(4, 0));
        initialQuadric = vector<Matrix>(mesh.get_vertices_num() + 1, Matrix(0.0));
        cout << "=========Calculate triangle plane========" << endl;
      //  compute_triangle_plane(mesh, trianglePlane);
        cout << "=========Calculate initial QEM========" << endl;
      //  compute_initial_QEM(mesh, trianglePlane);
        compute_plane_and_QEM(tree.get_root(),mesh,tree.get_subdivision(),tree);
    }
    time.start();
    while (1)
    {
        simplification_round = params.get_contracted_edges_num(); //checked edges

        /// HERE YOU NEED TO DEFINE A PROCEDURE FOR SIMPLIFY THE TIN BY USING THE SPATIAL INDEX
        this->simplify_compute(tree.get_root(), mesh, cache, tree.get_subdivision(), params, tree);

        cout << "Num of edges enqueued:" << params.get_sum_edge_queue_sizes() << endl;
        // PARTIAL SIMPLIFICATION STATS
        cerr << "=== end-of-round " << round << ") --> contracted edges: ";
        cerr << params.get_contracted_edges_num() - simplification_round << endl;
        round++;

        if (cli.debug_mode)
        {
            time.stop();
            time.print_elapsed_time("   [TIME] executing a simplification round: ");
            time.start();
            //  cerr << "   [RAM] peak for executing a simplification round: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
        }

        if (simplification_round == params.get_contracted_edges_num())
            break;

        cache.reset();
        //    if(round==2)
        //      break;
    }
    time.stop();
    if (!cli.debug_mode)
        time.print_elapsed_time("[TIME] Edge contraction simplification: ");
    //else
    //params.print_simplification_partial_timings();
    //params.print_simplification_counters();

    //  cerr << "[RAM peak] for contracting a simplicial complex: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;

    cerr << "[MEMORY] peak for Simplification: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;

    /// finally we have to update/compress the mesh and the tree
    Contraction_Simplifier::update_mesh_and_tree(tree, mesh, params, cli);
    cerr << "[MEMORY] peak for mesh and tree updating: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
}

void Contraction_Simplifier::simplify_compute(Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, Spatial_Subdivision &division, contraction_parameters &params, PRT_Tree &tree)
{
    if (n.is_leaf())
    {
        if (params.is_QEM() == true)
            simplify_leaf_QEM(n, mesh, cache, params, tree);
        else
            simplify_leaf(n, mesh, cache, params, tree);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            if (n.get_son(i) != NULL)
            {
                simplify_compute(*n.get_son(i), mesh, cache, division, params, tree);
            }
        }
    }
}

void Contraction_Simplifier::simplify_compute_parallel(Mesh &mesh, Spatial_Subdivision &division, contraction_parameters &params, PRT_Tree &tree)
{

    //First part for internal edges
    // #pragma omp parallel for // schedule(dynamic,1)
    //     for (unsigned i = 0; i < tree.get_leaves_number(); i++)
    //     {

    //         Node_V *leaf = tree.get_leaf(i);

    //         simplify_leaf(*leaf, mesh, cache, params, tree);

    //     }
    //Second part for cross edges

    //boost::dynamic_bitset<> conflict_nodes(tree.get_leaves_number());
    ivect nodes_status(tree.get_leaves_number(), 0);
    bool processed = false;
    int processed_node = 0;
    do
    {
        //  processed = false;

#pragma omp parallel for reduction(+ \
                                   : processed_node) // schedule(dynamic,1)
        for (unsigned i = 0; i < tree.get_leaves_number(); i++)
        {
            Node_V *leaf = tree.get_leaf(i);
            if (nodes_status[i] == -1)
                continue;
            if (!leaf->indexes_vertices())
            {
                processed_node++;
                nodes_status[i] = -1;
                continue;
            }
            //check the array of conflict_nodes
            // if nodes_status[i]==1, then continue
            //   cout << "Current leaf node:" << i << " on thread " << omp_get_thread_num() << endl;
            omp_set_lock(&(l_locks[i]));
            if (nodes_status[i] != 0)
            {
                omp_unset_lock(&(l_locks[i]));
                //   cout<<"Node "<<i<<" is skipped."<<endl;
                continue;
            }
            // else we set the conflict nodes of leaf[i] to be 1 in nodes_status
            else
            {
                nodes_status[i] = 2;
                omp_unset_lock(&(l_locks[i]));
                iset conflicts = conflict_leafs[i];
                // Check if the conflict nodes were set to 1 already
                // bool shared_conflicts = false;
                bool cannot_process = false;
                for (iset_iter it = conflicts.begin(); it != conflicts.end(); it++)
                {
                    //  cout<<"set leaf node:"<<*it<<" on thread "<<omp_get_thread_num()<<endl;
                    omp_set_lock(&(l_locks[*it])); // leafs should be locked when being checked and updated.
                    int status = 0;
                    //   #pragma omp atomic read
                    status = nodes_status[*it];
                    if (status == 1 || status == 2) //it is conflicted with another node being processed
                    {
                        // cout<<"conflict node id:"<<*it<<" with "<<i<<" on thread "<<omp_get_thread_num()<<endl;
                        omp_unset_lock(&(l_locks[*it]));
                        //it++; No need, the current one should not be change to 0
                        for (iset_iter it2 = conflicts.begin(); it2 != it; it2++)
                        {
                            //     cout<<"unset leaf node:"<<*it2<<" on thread "<<omp_get_thread_num()<<endl;
                            omp_set_lock(&(l_locks[*it2]));
                            if (nodes_status[*it2] == 1)
                                nodes_status[*it2] = 0;
                            omp_unset_lock(&(l_locks[*it2]));
                        }
                        //omp_unset_lock(&(l_locks[*it]));
                        // unset the locks that have been set
                        cannot_process = true;
                        //cout<<"neighbor conflict"<<endl;

                        break;
                    }
                    else if (status == 0)
                    {
                        nodes_status[*it] = 1;
                    }
                    omp_unset_lock(&(l_locks[*it]));
                }
                if (cannot_process == true)
                {
                    omp_set_lock(&(l_locks[i]));
                    nodes_status[i] = 0;
                    omp_unset_lock(&(l_locks[i]));
                    continue;
                }

                // Node_V *leaf = tree.get_leaf(i);
                //  cout<<"Node "<<i<<" will be processed."<<endl;
                // omp_set_lock(&(l_locks[i]));
                // // set nodes_status[i]=2 when node is being processed
                // nodes_status[i] = 2;
                // omp_unset_lock(&(l_locks[i]));

                //   processed = true;
                processed_node = processed_node + 1;
                // cout << "Start simplification" << endl;
                if (params.is_QEM() == true)
                    simplify_leaf_cross_QEM(*leaf, i, mesh, params, tree);
                else
                    simplify_leaf_cross(*leaf, i, mesh, params, tree);

                // cout << "Finish simplification" << endl;
                //set nodes_status[i]=-1 after processing

                omp_set_lock(&(l_locks[i]));
                nodes_status[i] = -1;
                omp_unset_lock(&(l_locks[i]));
                //cout<<"unset leaf node:"<<i<<" on thread "<<omp_get_thread_num()<<endl;

                for (iset_iter it = conflicts.begin(); it != conflicts.end(); it++)
                {
                    //  cout<<" set leaf lock "<<*it<<" on thread "<<omp_get_thread_num()<<endl;
                    omp_set_lock(&(l_locks[*it]));
                    int status = 0;
                    //   #pragma omp atomic read
                    status = nodes_status[*it];
                    if (status == 1)
                    {
                        nodes_status[*it] = 0;
                    }
                    //cout<<"unset leaf node:"<<*it<<" on thread "<<omp_get_thread_num()<<endl;
                    omp_unset_lock(&(l_locks[*it]));
                }
            }
        }
        // cout << "finished one for loop" << endl;
        cerr << "[MEMORY] peak for a simplification round: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
        cerr << "Number of processed nodes:" << processed_node << endl;
        //} while (processed == true);
    } while (processed_node != tree.get_leaves_number());
}

void Contraction_Simplifier::simplify_leaf_QEM(Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params, PRT_Tree &tree)
{

    if (!n.indexes_vertices())
        return;

    itype v_start = n.get_v_start();
    itype v_end = n.get_v_end();
    itype v_range = v_end - v_start;

    //cout<<"Simplification in leaf."<<endl;
    // leaf_VT local_vts(v_range,VT());
    boost::dynamic_bitset<> is_v_border(v_end-v_start);
    // n.get_VT(local_vts,mesh);
    leaf_VT &local_vts = get_VTS(n, mesh, cache, tree, params, is_v_border);
    // Create a priority queue of candidate edges
    edge_queue edges;
    find_candidate_edges_QEM(n, mesh, local_vts, edges, params);
    map<vector<int>, double> updated_edges;
    int edge_num = edges.size();
    int edges_contracted_leaf = 0;
    //cout << "Edge number:" << edges.size() << endl;
    params.add_edge_queue_size(edges.size());
    while (!edges.empty())
    {
        Geom_Edge *current = edges.top();
        ivect e = current->edge;
        //    cout<<"Start contraction."<<endl;
        //  cout<<"Edge error:"<<current->val<<endl;

        edges.pop();

        if (mesh.is_vertex_removed(e[0]) || mesh.is_vertex_removed(e[1]))
        {
            //       cout<<"skip current edge."<<endl;
            //   cout<<"[DEBUG] edge not complete: "<<e[0]<<", "<<e[1]<<endl;
            //   cout<<"Vertex removed"<<endl;
            delete current;
            // if(edges_contracted_leaf>edge_num*0.2)
            // cout<<"Num of deleted in a leaf"<<edges_contracted_leaf<<"; 20% of the queue num:"<<edge_num*0.2<<endl;
            continue;
        }

        ivect sorted_e = e;
        sort(sorted_e.begin(), sorted_e.end());
        auto it = updated_edges.find(sorted_e);
        if (it != updated_edges.end())
        {
            //int tmp=-1;
            // double error = compute_error(e[0],e[1],mesh,tmp);
            if (fabs(it->second - current->val) > Zero)
            {
                // cout<<"skip current edge."<<endl;
                // cout<<"[DEBUG] edge: "<<sorted_e[0]<<", "<<sorted_e[1]<<"; updated error: "<<it->second<<"old error: "<<current->val<<endl;
                delete current;
                continue;
            }
        }

        ET et(-1, -1);
        VT *vt0 = NULL, *vt1 = NULL;
        Node_V *outer_v_block = NULL;
        bool v1_is_border=false, v2_is_border=false;

        //get_edge_relations(e, et, vt0, vt1, outer_v_block, n, mesh, local_vts, cache, params, tree);
        get_edge_relations(e,et,vt0,vt1,v1_is_border,v2_is_border,outer_v_block,n,mesh,local_vts,is_v_border,cache,params,tree);

        if (link_condition(e[0], e[1], *vt0, *vt1, et, mesh)&&valid_boundary_condition(e[0],e[1],*vt0,*vt1,et,v1_is_border,v2_is_border,mesh)&&not_fold_over(e[0], e[1], *vt0, *vt1, et, mesh))
        {
            contract_edge(e, et, *vt0, *vt1, *outer_v_block, edges, n, mesh, cache, params,updated_edges);
            edges_contracted_leaf++;
            // break;
        }
        //cout << "Number of edges remaining:" << edges.size() << endl;

       delete current; 
    }

    //update the cache if the vts have been stored already.
    // if(cache.find(v_start) != cache.end()){
    //     cache.update(v_start,local_vts);
    // }
}

void Contraction_Simplifier::simplify_leaf(Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params, PRT_Tree &tree)
{

    if (!n.indexes_vertices())
        return;

    itype v_start = n.get_v_start();
    itype v_end = n.get_v_end();
    itype v_range = v_end - v_start;

    //cout<<"Simplification in leaf."<<endl;
    // leaf_VT local_vts(v_range,VT());
    // n.get_VT(local_vts,mesh);
    boost::dynamic_bitset<> is_v_border(v_end-v_start);
    ////// THIS method should be used when cache is used in parallel computing
    leaf_VT &local_vts = get_VTS(n, mesh, cache, tree, params, is_v_border);
    //// Here we use the simple get_VT to check if that is the cause of slow computing.

    // leaf_VT local_vts(v_range,VT());
    // n.get_VT(local_vts,mesh);

    // Create a priority queue of candidate edges
    edge_queue edges;

    find_candidate_edges(n, mesh, local_vts, edges, params);
    map<vector<int>, double> updated_edges;
    int edge_num = edges.size();
    int edges_contracted_leaf = 0;
    // cout<<"Edge number:"<<edges.size()<<endl;
    // cout<<"Number of threads used: "<<omp_get_num_threads()<<endl;
    // cout<<"Current thread id: "<<omp_get_thread_num()<<endl;

    params.add_edge_queue_size(edges.size());
    while (!edges.empty())
    {
        Geom_Edge *current = edges.top();
        ivect e = current->edge;
        //    cout<<"Start contraction."<<endl;
        //  cout<<"Edge Length:"<<current->val<<endl;

        edges.pop();
        // NOTE: Should be atomic when we have cross edges.
        if (mesh.is_vertex_removed(e[0]) || mesh.is_vertex_removed(e[1]))
        {

            //   cout<<"Vertex removed"<<endl;
            // cout<<"skip current edge"<<endl;
            delete current;
            // if(edges_contracted_leaf>edge_num*0.2)
            // cout<<"Num of deleted in a leaf"<<edges_contracted_leaf<<"; 20% of the queue num:"<<edge_num*0.2<<endl;
            continue;
        }

        ET et(-1, -1);
        VT *vt0 = NULL, *vt1 = NULL;
        Node_V *outer_v_block = NULL;
        bool v1_is_border=false, v2_is_border=false;
        
        //get_edge_relations(e, et, vt0, vt1, outer_v_block, n, mesh, local_vts, cache, params, tree);
        
        get_edge_relations(e,et,vt0,vt1,v1_is_border,v2_is_border,outer_v_block,n,mesh,local_vts,is_v_border,cache,params,tree);

        if (link_condition(e[0], e[1], *vt0, *vt1, et, mesh)&&valid_boundary_condition(e[0],e[1],*vt0,*vt1,et,v1_is_border,v2_is_border,mesh)&&not_fold_over(e[0], e[1], *vt0, *vt1, et, mesh))
        {
            contract_edge(e, et, *vt0, *vt1, *outer_v_block, edges, n, mesh, cache, params,updated_edges);
            edges_contracted_leaf++;
            // break;
        }
        delete current;
    }

}

void Contraction_Simplifier::simplify_leaf_cross(Node_V &n, int n_id, Mesh &mesh, contraction_parameters &params, PRT_Tree &tree)
{

    if (!n.indexes_vertices())
        return;

    itype v_start = n.get_v_start();
    itype v_end = n.get_v_end();
    itype v_range = v_end - v_start;

    //cout<<"Simplification in leaf."<<endl;
    // leaf_VT local_vts(v_range,VT());
    // n.get_VT(local_vts,mesh);
    boost::dynamic_bitset<> is_v_border(v_end - v_start);
    leaf_VT local_vts(v_range, VT());
    n.get_VT_and_border(local_vts, is_v_border, mesh);

    // Create a priority queue of candidate edges
    edge_queue edges;
    // if (params.is_parallel())
    // {
    //     find_candidate_edges_parallel(n, mesh, local_vts, edges, params, true);
    // }
    // else
    find_candidate_edges(n, mesh, local_vts, edges, params);
    map<vector<int>, double> updated_edges;

    int edge_num = edges.size();
    int edges_contracted_leaf = 0;
    // cout << "Edge number:" << edges.size() << endl;
    // cout<<"Number of threads used: "<<omp_get_num_threads()<<endl;
    // cout<<"Current thread id: "<<omp_get_thread_num()<<endl;

    params.add_edge_queue_size(edges.size());
    map<int, leaf_VT> local_cache;
    while (!edges.empty())
    {
        Geom_Edge *current = edges.top();
        ivect e = current->edge;
        //    cout<<"Start contraction."<<endl;
        //  cout<<"Edge Length:"<<current->val<<endl;

        edges.pop();

        if (mesh.is_vertex_removed(e[0]) || mesh.is_vertex_removed(e[1]))
        {

            delete current;
            continue;
        }

        ET et(-1, -1);
        VT *vt0 = NULL, *vt1 = NULL;
        Node_V *outer_v_block = NULL;
        bool v1_is_border = false, v2_is_border = false;
        get_edge_relations(e, et, vt0, vt1, v1_is_border, v2_is_border, outer_v_block, n, mesh, local_vts, is_v_border, local_cache, params, tree);

        // if(params.is_parallel()){
        VV vv_locks;
        if (link_condition(e[0], e[1], *vt0, *vt1, et, n, *outer_v_block, vv_locks, mesh) && valid_boundary_condition(e[0], e[1], *vt0, *vt1, et, v1_is_border, v2_is_border, mesh) && not_fold_over(e[0], e[1], *vt0, *vt1, et, mesh))
        {
            contract_edge(e, et, *vt0, *vt1, *outer_v_block, edges, n, mesh, params, updated_edges);
            edges_contracted_leaf++;
            // break;

            // A new step for cross edge case
            // Check possible new conflict nodes by checking the vv_locks
            // vv_locks stores all the vertices in the VV(v0) & VV(v1) that are not contained by n or outer_v_block
            update_conflict_nodes(vv_locks, n_id, tree);
        }
        // for (iset_iter it = vv_locks.begin(); it != vv_locks.end(); it++)
        // {
        //     omp_unset_lock(&(v_locks[*it - 1]));
        // }
        // }
        // else
        // {
        // if (link_condition(e[0], e[1], *vt0, *vt1, et, mesh))
        // {
        //     contract_edge(e, et, *vt0, *vt1, *outer_v_block, edges, n, mesh, cache, params);
        //     edges_contracted_leaf++;
        //     // break;
        // }
        // }
        delete current;
    }

    // leaf_VV vvs;
    // n.get_VV(vvs,mesh);
}

void Contraction_Simplifier::simplify_leaf_cross_QEM(Node_V &n, int n_id, Mesh &mesh, contraction_parameters &params, PRT_Tree &tree)
{

    if (!n.indexes_vertices())
        return;

    itype v_start = n.get_v_start();
    itype v_end = n.get_v_end();
    itype v_range = v_end - v_start;

    //cout<<"Simplification in leaf."<<endl;
    boost::dynamic_bitset<> is_v_border(v_end - v_start);
    leaf_VT local_vts(v_range, VT());
    n.get_VT_and_border(local_vts, is_v_border, mesh);

    // Create a priority queue of candidate edges
    edge_queue edges;

    find_candidate_edges_QEM(n, mesh, local_vts, edges, params);
    map<vector<int>, double> updated_edges;
    int edge_num = edges.size();
    int edges_contracted_leaf = 0;
    // cout << "Edge number:" << edges.size() << endl;
    // cout<<"Number of threads used: "<<omp_get_num_threads()<<endl;
    // cout<<"Current thread id: "<<omp_get_thread_num()<<endl;

    params.add_edge_queue_size(edges.size());
    map<int, leaf_VT> local_cache;
    while (!edges.empty())
    {
        Geom_Edge *current = edges.top();
        ivect e = current->edge;
        //    cout<<"Start contraction."<<endl;
        //   cout<<"Edge Length:"<<current->val<<endl;

        edges.pop();

        if (mesh.is_vertex_removed(e[0]) || mesh.is_vertex_removed(e[1]))
        {

            delete current;
            continue;
        }
        ivect sorted_e = e;
        sort(sorted_e.begin(), sorted_e.end());
        auto it = updated_edges.find(sorted_e);
        if (it != updated_edges.end())
        {
            //int tmp=-1;
            // double error = compute_error(e[0],e[1],mesh,tmp);
            if (fabs(it->second - current->val) > Zero)
            {
                delete current;
                continue;
            }
        }
        ET et(-1, -1);
        VT *vt0 = NULL, *vt1 = NULL;
        Node_V *outer_v_block = NULL;
        bool v1_is_border = false, v2_is_border = false;

        get_edge_relations(e, et, vt0, vt1, v1_is_border, v2_is_border, outer_v_block, n, mesh, local_vts, is_v_border, local_cache, params, tree);

        // if(params.is_parallel()){
        VV vv_locks;
        if (link_condition(e[0], e[1], *vt0, *vt1, et, n, *outer_v_block, vv_locks, mesh) && valid_boundary_condition(e[0], e[1], *vt0, *vt1, et, v1_is_border, v2_is_border, mesh) && not_fold_over(e[0], e[1], *vt0, *vt1, et, mesh))
        {
            contract_edge(e, et, *vt0, *vt1, *outer_v_block, edges, n, mesh, params, updated_edges);
            edges_contracted_leaf++;
            // break;
            n_id = v_in_leaf[e[0]];
            // A new step for cross edge case
            // Check possible new conflict nodes by checking the vv_locks
            // vv_locks stores all the vertices in the VV(v0) & VV(v1) that are not contained by n or outer_v_block
            update_conflict_nodes(vv_locks, n_id, tree);
        }
        // for (iset_iter it = vv_locks.begin(); it != vv_locks.end(); it++)
        // {
        //     omp_unset_lock(&(v_locks[*it - 1]));
        // }
        delete current;
       // delete vt0, vt1, outer_v_block;
    }

    // leaf_VV vvs;
    // n.get_VV(vvs,mesh);
}

void Contraction_Simplifier::get_edge_relations(ivect &e, ET &et, VT *&vt0, VT *&vt1, Node_V *&outer_v_block, Node_V &n, Mesh &mesh, leaf_VT &vts, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params, PRT_Tree &tree)
{

    //cout<<"[NOTICE]get edge relation"<<endl;
    outer_v_block = NULL;
    /// inverted order as I only need the block indexing v1

    if (e[1] > e[0])
    {
        vt1 = Contraction_Simplifier::get_VT(e[1], n, mesh, vts, cache, tree, outer_v_block, params);
        vt0 = Contraction_Simplifier::get_VT(e[0], n, mesh, vts, cache, tree, outer_v_block, params);
    }
    else
    {

        vt0 = Contraction_Simplifier::get_VT(e[0], n, mesh, vts, cache, tree, outer_v_block, params);
        vt1 = Contraction_Simplifier::get_VT(e[1], n, mesh, vts, cache, tree, outer_v_block, params);
    }

    Contraction_Simplifier::get_ET(e, et, n, mesh, vts);
}

//A new version of get_edge_relations for parallel computing without global cache
void Contraction_Simplifier::get_edge_relations(ivect &e, ET &et, VT *&vt0, VT *&vt1, bool &v1_is_border, bool &v2_is_border, Node_V *&outer_v_block, Node_V &n, Mesh &mesh, leaf_VT &vts, boost::dynamic_bitset<> is_border_edge, map<int, leaf_VT> &cache, contraction_parameters &params, PRT_Tree &tree)
{

    //cout<<"[NOTICE]get edge relation"<<endl;
    outer_v_block = NULL;

    // Using local cache
  //  if (e[1] > e[0])
  //  {

        vt1 = Contraction_Simplifier::get_VT(e[1], n, mesh, vts, cache, tree, outer_v_block, params);
        vt0 = Contraction_Simplifier::get_VT(e[0], n, mesh, vts, cache, tree, outer_v_block, params);
        v2_is_border = is_border_edge[e[1] - n.get_v_start()];
       
        if (n.indexes_vertex(e[0]))
        {
            v1_is_border = is_border_edge[e[0] - n.get_v_start()];
        }
        else
        {
            for (auto it = vt0->begin(); it != vt0->end(); it++)
            {
                Triangle &t = mesh.get_triangle(*it);
                int v_pos = t.vertex_index(e[0]);
                for (int v1 = 1; v1 < t.vertices_num(); v1++)
                {
                    if (t.is_border_edge((v1 + v_pos) % t.vertices_num()))
                    {
                        v1_is_border = true;
                        break;
                    }
                }
            }
        }
   // }
    // else
    // {
 
    //     vt0 = get_VT(e[0], n, mesh, vts, cache, tree, outer_v_block, params);
    //     vt1 = get_VT(e[1], n, mesh, vts, cache, tree, outer_v_block, params);

    //     v1_is_border = is_border_edge[e[0] - n.get_v_start()];
    //     if (n.indexes_vertex(e[1]))
    //     {
    //         v2_is_border = is_border_edge[e[1] - n.get_v_start()];
    //     }
    //     else
    //     {
    //         for (auto it = vt1->begin(); it != vt1->end(); it++)
    //         {
    //             Triangle &t = mesh.get_triangle(*it);
    //             int v_pos = t.vertex_index(e[1]);
    //             for (int v1 = 1; v1 < t.vertices_num(); v1++)
    //             {
    //                 if (t.is_border_edge((v1 + v_pos) % t.vertices_num()))
    //                 {
    //                     v2_is_border = true;
    //                     break;
    //                 }
    //             }
    //         }
    //     }
    // }

    Contraction_Simplifier::get_ET(e, et, n, mesh, vts);
}

void Contraction_Simplifier::get_edge_relations(ivect &e, ET &et, VT *&vt0, VT *&vt1,bool& v1_is_border, bool& v2_is_border, Node_V *&outer_v_block,
                                                Node_V &n, Mesh &mesh, leaf_VT &vts,boost::dynamic_bitset<> is_border_edge, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params,PRT_Tree &tree){
    //cout<<"[NOTICE]get edge relation"<<endl;
    outer_v_block = NULL;
    /// inverted order as I only need the block indexing v1
    // if(e[1]>e[0]){
    vt1 = get_VT(e[1],n,mesh,vts,cache,tree,outer_v_block,params);
    vt0 = get_VT(e[0],n,mesh,vts,cache,tree,outer_v_block,params);
   
    v2_is_border=is_border_edge[e[1]-n.get_v_start()];
    if(n.indexes_vertex(e[0])){
        v1_is_border=is_border_edge[e[0]-n.get_v_start()];
    }
    else{
        for(auto it=vt0->begin();it!=vt0->end();it++){
            Triangle& t=mesh.get_triangle(*it);
            int v_pos=t.vertex_index(e[0]);
            for(int v1=1;v1<t.vertices_num();v1++){
                if(t.is_border_edge((v1+v_pos)%t.vertices_num()))
                {
                    v1_is_border=true;
                    break;
                }
            }

        }
    }
   
    Contraction_Simplifier::get_ET(e,et,n,mesh,vts);                                              


}



VT *Contraction_Simplifier::get_VT(int v_id, Node_V &n, Mesh &mesh, leaf_VT &vts, LRU_Cache<int, leaf_VT> &cache,
                                   PRT_Tree &tree, Node_V *&v_block, contraction_parameters &params)
{
    int local_index;
    bool debug = false;
    if (n.indexes_vertex(v_id))
    {
        if (debug)
            cout << "[get_VT] " << v_id << " -> LOCAL VERTEX " << n << endl;

        local_index = v_id - n.get_v_start();
        VT *vt = &(vts[local_index]);
        Contraction_Simplifier::clean_coboundary(*vt, mesh);
        v_block = &n;
        return vt;
    }
    else
    {
        if (debug) // if v is external
            cout << "[get_VT] " << v_id << " -> EXTERNAL VERTEX " << n << endl;

        tree.get_leaf_indexing_vertex(tree.get_root(), v_id, v_block); // the vertex should be protected by a lock in the future
        local_index = v_id - v_block->get_v_start();

        LRU_Cache<int, leaf_VT>::mapIt it_c = cache.end();

        {
            it_c = cache.find(v_block->get_v_start()); //First check in the cache
            if (it_c == cache.end())                   //if not in the cache
            {
                if (debug)
                    cout << "    -> LEAF BLOCK OUTSIDE CACHE - REGEN " << *v_block << endl;

                leaf_VT lVT;
                v_block->get_VT(lVT, mesh);
                it_c = cache.insert(v_block->get_v_start(), lVT);
            }
            else
            {
                if (debug)
                {
                    cout << "    -> LEAF BLOCK IN CACHE - CLEAN " << *v_block << endl;
                }

                // if(debug/* || v_id ==2355*/)
                //     cout<<"num_elem_in_vt: "<<get_num_elements_in_container_of_containers((it_c->second)[local_index])<<endl;
                //                print_container_of_containers_content("VTop(2355) ",(it_c->second)[local_index]);
            }
            Contraction_Simplifier::clean_coboundary((it_c->second)[local_index], mesh);

            if (debug)
                cout << "[NOTICE]cleaned coboundary" << endl;
        }

        return &(it_c->second)[local_index];
    }
}

// New version of get_VT() without using global cache
VT *Contraction_Simplifier::get_VT(int v_id, Node_V &n, Mesh &mesh, leaf_VT &vts, map<int, leaf_VT> &cache,
                                   PRT_Tree &tree, Node_V *&v_block, contraction_parameters &params)
{
    int local_index;
    bool debug = false;

    if (n.indexes_vertex(v_id))
    {
        if (debug)
            cout << "[get_VT] " << v_id << " -> LOCAL VERTEX " << n << endl;

        local_index = v_id - n.get_v_start();
        VT *vt = &(vts[local_index]);
        Contraction_Simplifier::clean_coboundary(*vt, mesh);
        v_block = &n;
        return vt;
    }
    else
    {
        if (debug) // if v is external
            cout << "[get_VT] " << v_id << " -> EXTERNAL VERTEX " << n << endl;

        //UPDATE: we can use v_in_leaf array to directly locate the node containing it
        //tree.get_leaf_indexing_vertex(tree.get_root(), v_id, v_block); // the vertex should be protected by a lock in the future
        int n_index = v_in_leaf[v_id];
        v_block = tree.get_leaf(n_index);
        local_index = v_id - v_block->get_v_start();
        //LRU_Cache<int, leaf_VT>::mapIt it_c = cache.end();
        map<int, leaf_VT>::iterator it_c = cache.end();
        VT *vt = NULL;
        it_c = cache.find(v_block->get_v_start()); //First check in the cache
        if (it_c == cache.end())                   //if not in the cache
        {
            if (debug)
                cout << "    -> LEAF BLOCK OUTSIDE CACHE - REGEN " << *v_block << endl;

            leaf_VT lVT;
            v_block->get_VT(lVT, mesh);
            cache.insert({v_block->get_v_start(), lVT});
            vt = &(cache[v_block->get_v_start()][local_index]);
        }
        else
        {
            if (debug)
            {
                cout << "    -> LEAF BLOCK IN CACHE - CLEAN " << *v_block << endl;
            }
            vt = &((it_c->second)[local_index]);
        }
        Contraction_Simplifier::clean_coboundary(*vt, mesh);

        if (debug)
            cout << "[NOTICE]cleaned coboundary" << endl;
        return vt;
    }
}

leaf_VT &Contraction_Simplifier::get_VTS(Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache,
                                         PRT_Tree &tree, contraction_parameters &params, boost::dynamic_bitset<>& is_v_border)
{
    int local_index;
 
    
    
    LRU_Cache<int, leaf_VT>::mapIt it_c = cache.end();

    
        it_c = cache.find(n.get_v_start()); //First check in the cache
        if (it_c == cache.end())            //if not in the cache
        {
            leaf_VT local_vts;
            n.get_VT_and_border(local_vts, is_v_border, mesh);
            it_c = cache.insert(n.get_v_start(), local_vts);
        }
        else{
            leaf_VT local_vts = it_c->second;
            for(int i=0;i<local_vts.size();i++){

            VT vt = local_vts[i];
            for(int j=0;j<vt.size();j++){
            if(!is_v_border[i])
              {
                Triangle t= mesh.get_triangle(vt[j]);
                int v = t.vertex_index(n.get_v_start()+i);
                for(int v1=1; v1<t.vertices_num(); v1++)
            {
                if(t.is_border_edge((v1+v)%t.vertices_num()))
                {
                    is_v_border[i] = true;
                    break;
                }
            }
            }
            }   
        }
        }
        
    return it_c->second;
}

void Contraction_Simplifier::update_mesh_and_tree(PRT_Tree &tree, Mesh &mesh, contraction_parameters &params, cli_parameters &cli)
{
    Timer time;

    ///  UPDATE OF MESH AND TREE
    ivect new_v_positions;
    ivect new_t_positions;
    ivect surviving_vertices;

    time.start();
    //    cerr<<"[TREE] compact vertices lists"<<endl;
    tree.compact_vertices_lists(tree.get_root(), mesh, surviving_vertices);
    time.stop();
    time.print_elapsed_time("[TIME] Compact tree vertices lists: ");
    cerr << "[MEMORY] peak for compacting tree vertices lists: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;

    //    print_container_content("surviving vertices: ",surviving_vertices);
    //    mesh.print_mesh(cout);
    //    int a; cin>>a;

    time.start();
    //    cerr<<"[MESH] compact"<<endl;
    Mesh_Updater mu;
    cout << "number of surviving vertices:" << surviving_vertices.size() << endl;
    mu.clean_vertices_array(mesh, new_v_positions, surviving_vertices);
    cout << "number of deleted triangles:" << params.get_counter() << endl;
    /// NEW: the update_and_compact procedure check internally if we have removed all the top d-simplices
    bool all_deleted = mu.update_and_clean_triangles_arrays(mesh, new_v_positions, new_t_positions, params.get_counter());
    time.stop();
    time.print_elapsed_time("[TIME] Compact and update mesh: ");
    cerr << "[MEMORY] peak for compacting and updating the mesh: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;

    // cerr<<"[STAT] mesh "<<endl;
    // cerr<<"------------"<<endl;
    // mesh.print_mesh_stats(cerr);
    // cerr<<"------------"<<endl;

   
    //    cerr<<"[TREE] update indices in the tree"<<endl;
    int index_counter = 1;
    time.start();
    tree.update_vertex_index(tree.get_root(), new_v_positions, index_counter);
    time.stop();
    time.print_elapsed_time("[TIME] Update tree structure (merging blocks): ");

    time.start();
    tree.reinsert_triangles();
    time.stop();
    time.print_elapsed_time("[TIME] Update tree (triangles): ");
    //    Reindexer r;
    //    r.reorganize_index_and_mesh(tree,mesh,false);

    //cerr << "[RAM peak] for updating the mesh and the tree: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
}

void Contraction_Simplifier::find_candidate_edges_QEM(Node_V &n, Mesh &mesh, leaf_VT &vts, edge_queue &edges, contraction_parameters &params)
{
    map<ivect, coord_type> edge_map;
    ivect e;

    for (RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const &t_id = itPair.first;
        if (mesh.is_triangle_removed(*t_id))
        {
            //  cout<<"triangle removed"<<endl;
            continue;
        }

        Triangle &t = mesh.get_triangle(*t_id);
        // if(!n.indexes_triangle_vertices(t))
        //     continue;

        for (int i = 0; i < t.vertices_num(); i++)
        {
            t.TE(i, e);
            if((!n.indexes_vertex(e[0]))&&(!n.indexes_vertex(e[1])))
                 continue;
            int new_vertex_pos = -1;
            double error = compute_error(e[0], e[1], mesh, new_vertex_pos);

            if(params.output_stats()&&n.indexes_vertex(e[1])){
                #pragma omp critical
                    edge_costs_output[make_pair(e[0],e[1])] = error;
            }

         //   assert(new_vertex_pos != -1);
            if (new_vertex_pos == 1)
                {
                    int tmp = e[1];
                    e[1] = e[0];
                    e[0] = tmp;
            }  // Ensure that e[1] is always the vertex to be removed, and 
               // it should be contained by the current node
            if (n.indexes_vertex(e[1])) // e (v1,v2) is a candidate edge if at least v2 is in n
            {


                map<ivect, coord_type>::iterator it = edge_map.find(e);
                // cout<<e[0]<<" and "<<e[1]<<endl;
                if (it == edge_map.end())
                {

                    //  Edge e((*it)[0],(*it)[1]);
                    edge_map[e] = error;
                    //Edge edge_obj(e[0],e[1]);
                    if (params.no_limit())
                    {
                        //cout<<"["<<e[0]-1<<","<<e[1]-1<<"]  Error will be introduced: "<<error<<endl;

                        edges.push(new Geom_Edge(e, error));
                    }
                    else
                    {
                        if ((error - params.get_maximum_limit() )< Zero)
                        {
                            //      cout<<"["<<e[0]-1<<","<<e[1]-1<<"]  Error will be introduced: "<<error<<endl;
                            edges.push(new Geom_Edge(e, error));
                            //     cout<<"ENQUEUE"<<endl;
                        }
                    }
                }
            }
        }
    }
    // cout << "**** [Number] " << edges.size() << " edges enqueued. Start simplification.****" << endl;
    // cout << "number of remaining triangles:" << t_count << endl;

    // cout<<"======NEXT NODE======"<<endl;
}

void Contraction_Simplifier::compute_initial_QEM(Mesh &mesh, vector<dvect> &planes)
{
#pragma omp parallel for
    for (int i = 1; i <= mesh.get_triangles_num(); i++)
    {
        /* faces are triangles */
        //  #pragma omp parallel for
        for (int j = 0; j < 3; j++)
        {
            double *a = &(planes[i - 1][0]);
            //
            Matrix tmp = Matrix(a);
            //#pragma omp critical
           // omp_set_lock(&(v_locks[mesh.get_triangle(i).TV(j) - 1]));
            initialQuadric[mesh.get_triangle(i).TV(j)] += tmp;
           // omp_unset_lock(&(v_locks[mesh.get_triangle(i).TV(j) - 1]));
        }
    }
}

void Contraction_Simplifier::compute_triangle_plane(Mesh &mesh, vector<dvect> &trPl)
{

#pragma omp parallel for
    for (int i = 1; i <= mesh.get_triangles_num(); i++)
    {
        double coords[3][3];
        for (int v = 0; v < 3; v++)
        {
            Vertex v1 = mesh.get_vertex(mesh.get_triangle(i).TV(v));
            coords[0][v] = v1.get_x();
            coords[1][v] = v1.get_y();
            coords[2][v] = v1.get_z();
        }

        double a, b, c, m;

        a = (coords[1][1] - coords[1][0]) * (coords[2][2] - coords[2][0]) - (coords[2][1] - coords[2][0]) * (coords[1][2] - coords[1][0]);

        b = (coords[2][1] - coords[2][0]) * (coords[0][2] - coords[0][0]) - (coords[0][1] - coords[0][0]) * (coords[2][2] - coords[2][0]);

        c = (coords[0][1] - coords[0][0]) * (coords[1][2] - coords[1][0]) - (coords[1][1] - coords[1][0]) * (coords[0][2] - coords[0][0]);

        double tmp = a * a + b * b + c * c;
        m = sqrt(tmp);

        a = a / m;
        b = b / m;
        c = c / m;

        trPl[i - 1][0] = a;
        trPl[i - 1][1] = b;
        trPl[i - 1][2] = c;
        trPl[i - 1][3] = -1 * (a * coords[0][0] + b * coords[1][0] + c * coords[2][0]);
    }
}

double Contraction_Simplifier::compute_error(int v1, int v2, Mesh &mesh, int &new_vertex_pos)
{
    double min_error;
    Matrix q_bar;
    Matrix q_delta;

    /* computer quadric of virtual vertex vf */
    //cout<<"v1: "<<v1<<" and v2: "<<v2<<endl;
    q_bar = initialQuadric[v1] + initialQuadric[v2];

    q_delta = Matrix(q_bar[0], q_bar[1], q_bar[2], q_bar[3],
                     q_bar[4], q_bar[5], q_bar[6], q_bar[7],
                     q_bar[8], q_bar[9], q_bar[10], q_bar[11],
                     0, 0, 0, 1);

    Vertex vertex_1 = mesh.get_vertex(v1);
    double vx1 = vertex_1.get_x();
    double vy1 = vertex_1.get_y();
    double vz1 = vertex_1.get_z();

    Vertex vertex_2 = mesh.get_vertex(v2);
    double vx2 = vertex_2.get_x();
    double vy2 = vertex_2.get_y();
    double vz2 = vertex_2.get_z();

    double error1 = vertex_error(q_bar, vx1, vy1, vz1);
    double error2 = vertex_error(q_bar, vx2, vy2, vz2);

    min_error = std::min(error1, error2);

    if (fabs(error1 - min_error) < Zero)
    {
        new_vertex_pos = 0;
        min_error = error1;
    }
    else
    {
        new_vertex_pos = 1;
        min_error = error2;
    }

    //  min_error = vertex_error(q_bar, new_vertex[0], new_vertex[1], new_vertex[2]);
    if (min_error < Zero)
        min_error = 0;
    return min_error;
}

void Contraction_Simplifier::preprocess(PRT_Tree &tree, Mesh &mesh, cli_parameters &cli)
{
    if(cli.num_of_threads != 1){
    conflict_leafs.assign(tree.get_leaves_number(), iset());
    }
    //map<int, ivect> nodes_of_t;
    v_in_leaf.assign(mesh.get_vertices_num() + 1, -1);
    vector<ivect> cross_triangles(tree.get_leaves_number(), ivect());
    // Can be parallel
  //  cout << "preprocessing" << endl;
#pragma omp parallel for
    for (int i = 0; i < tree.get_leaves_number(); i++)
    {
        Node_V *n = tree.get_leaf(i);
        if (!n->indexes_vertices())
            continue;
        // UPDATE: check triangles and vertices in the current node in the same step
        // Visit all triangles in the current node
        // update nodes_of_t array
        // Visit all vertices in the current leaf node
        // update v_in_leaf array
        for (RunIteratorPair itPair = n->make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const &t_id = itPair.first;

            if (mesh.is_triangle_removed(*t_id))
                continue;
            Triangle t = mesh.get_triangle(*t_id);

            int num_v_indexed = 0;
            for (int j = 0; j < t.vertices_num(); j++)
            {
                if (n->indexes_vertex(t.TV(j)))
                {
                    num_v_indexed++;
                    //  cout<<t.TV(j)<<endl;
                    v_in_leaf[t.TV(j)] = i;
                }
            }
            if (num_v_indexed > 0 && num_v_indexed < 3)
            {
                cross_triangles[i].push_back(*t_id);
                //    omp_set_lock(&(t_locks[*t_id - 1]));
                //  #pragma omp critical
                //      nodes_of_t[*t_id].push_back(i);
                //  omp_unset_lock(&(t_locks[*t_id - 1]));
            }

            // if (n->partial_indexes_triangle_vertices(mesh.get_triangle(*t_id)))
            // {
            //     nodes_of_t[*t_id].push_back(i);
            // }
        }

        // for (auto it = n->v_array_begin_iterator(); it != n->v_array_end_iterator(); it++)
        // {
        //     v_in_leaf[*it] = i;
        // }
    }
if(cli.num_of_threads != 1){
#pragma omp parallel for
    for (int i = 0; i < tree.get_leaves_number(); i++)
    {
        Node_V *n = tree.get_leaf(i);
        if (!n->indexes_vertices())
            continue;
        for (auto it = cross_triangles[i].begin(); it != cross_triangles[i].end(); it++)
        {
            Triangle t = mesh.get_triangle(*it);
            for (int j = 0; j < t.vertices_num(); j++)
            {
                if (!n->indexes_vertex(t.TV(j)))
                {
                    conflict_leafs[i].insert(v_in_leaf[t.TV(j)]);
                    //  cout<<t.TV(j)<<endl;
                }
            }
        }
    }
}




    // cout<<"[DEBUG] print the conflict list"<<endl;
    // for(int i=0;i<conflict_leafs.size();i++){
    //     cout<<"The conflict list of node "<<i<<" :";
    //     for(iset_iter it=conflict_leafs[i].begin();it!=conflict_leafs[i].end();it++){
    //         cout<<*it<<", ";
    //     }
    //     cout<<endl;
    // }
}

void Contraction_Simplifier::update_conflict_nodes(VV &vv_locks, int n_id, PRT_Tree &tree)
{

    for (auto it = vv_locks.begin(); it != vv_locks.end(); it++)
    {
        int indexed_node_id;
        //indexed_node_id = tree.get_leaf_indexing_vertex_index(tree.get_root(), vv_locks[i]);
        indexed_node_id = v_in_leaf[*it];
        if (indexed_node_id == -1)
        {
            cout << "SHOULD NOT HAPPEN" << endl;
            int a;
            cin >> a;
        }
        if (conflict_leafs[n_id].find(indexed_node_id) == conflict_leafs[n_id].end())
        {
         //   cout<<"node "<<n_id<<" is now conflicted with "<<indexed_node_id<<endl;
            conflict_leafs[n_id].insert(indexed_node_id);
            conflict_leafs[indexed_node_id].insert(n_id);
        }
    }
}

void Contraction_Simplifier::generate_v_in_leaf(PRT_Tree &tree, Mesh &mesh, cli_parameters &cli)
{
    itype num_leaf = tree.get_leaves_number();
    for (int i = 0; i < num_leaf; i++)
    {
        Node_V *leaf = tree.get_leaf(i);
        for (auto it = leaf->v_array_begin_iterator(); it != leaf->v_array_end_iterator(); it++)
        {
            v_in_leaf[*it] = i;
        }
    }
}

bool Contraction_Simplifier::valid_boundary_condition(int v1, int v2, VT &vt1, VT &vt2, ET &et, bool v1_is_border, bool v2_is_border, Mesh &mesh)
{

    if (v1_is_border || v2_is_border)
    {
        //cout<<"border edge"<<endl;
        return false;
    }
    if (vt1.size() < 4 || vt2.size() < 4)
    {
        // cout<<"less than 4 triangles"<<endl;
        return false;
    }

    return true;
}

bool Contraction_Simplifier::not_fold_over(int v1, int v2, VT &vt1, VT &vt2, ET &et, Mesh &mesh)
{


    VT vt2_sub_et = vt2;
    ivect et_vec;
    et_vec.push_back(et.first);
    if (et.second != -1)
        et_vec.push_back(et.second);
    difference_of_vectors(vt2_sub_et, et_vec);
    for (int i = 0; i < vt2_sub_et.size(); i++)
    {
        Triangle t = mesh.get_triangle(vt2_sub_et[i]);
        int v2_pos = t.vertex_index(v2);
        ivect e;
        t.TE(v2_pos, e);
        bool same_side = Geometry_Wrapper::same_side_of_edge(v1, v2, e[0], e[1], mesh);

        if (!same_side)
        {
            return false;
        }
    }

    return true;
}

void Contraction_Simplifier::update_QEM(ivect &surviving_vertices, Mesh &mesh)
{
    vector<Matrix> newQuadric(surviving_vertices.size() + 1, Matrix(0.0));
    int v_counter = 1;
    for (ivect_iter it = surviving_vertices.begin(); it != surviving_vertices.end(); ++it)
    {
        // mesh.add_vertex(old_list[*it-1]);
        newQuadric[v_counter] = initialQuadric[*it];
        v_counter++;
    }
    initialQuadric = newQuadric;
    vector<Matrix>().swap(newQuadric);
}

void Contraction_Simplifier::compute_leaf_QEM(Node_V &n, Mesh &mesh, leaf_VT &vts)
{

    // itype v_begin = *(n.v_array_begin_iterator());
    // if(!initialQuadric[v_begin].is_zero())
    // return;
    for (auto it = n.v_array_begin_iterator(); it != n.v_array_end_iterator(); it++)
    {
        itype local_index = *it - n.get_v_start();
        for (auto it_t = vts[local_index].begin(); it_t != vts[local_index].end(); it_t++)
        {
            double *a = &(trianglePlane[(*it_t) - 1][0]);
            initialQuadric[*it] += Matrix(a);
        }
    }
}
void Contraction_Simplifier::compute_initial_QEM_parallel(PRT_Tree &tree, Mesh &mesh, vector<dvect> &planes)
{

#pragma omp parallel for
    for (int i = 0; i < tree.get_leaves_number(); i++)
    {
        Node_V *n = tree.get_leaf(i);
        if (!n->indexes_vertices())
            continue;

        for (RunIteratorPair itPair = n->make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const &t_id = itPair.first;

            Triangle t = mesh.get_triangle(*t_id);

            for (int j = 0; j < t.vertices_num(); j++)
            {
                if (n->indexes_vertex(t.TV(j)))
                {
                    double *a = &(trianglePlane[(*t_id) - 1][0]);
                    initialQuadric[t.TV(j)] += Matrix(a);
                }
            }
        }

        // itype v_start = n->get_v_start();
        // itype v_end = n->get_v_end();
        // itype v_range = v_end - v_start;
        //  leaf_VT local_vts(v_range, VT());
        //  n->get_VT(local_vts,mesh);
        //  compute_leaf_QEM(*n,mesh,local_vts);
    }
}

void Contraction_Simplifier::compute_initial_plane_and_QEM_parallel(PRT_Tree &tree, Mesh &mesh)
{

#pragma omp parallel for
    for (int i = 0; i < tree.get_leaves_number(); i++)
    {
        Node_V *n = tree.get_leaf(i);
        if (!n->indexes_vertices())
            continue;

        for (RunIteratorPair itPair = n->make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const &t_id = itPair.first;

            Triangle t = mesh.get_triangle(*t_id);
            if (!n->indexes_triangle_vertices(t))
                continue;
            else
            {
                dvect triangle_plane = dvect(4, 0);
                double coords[3][3];
                for (int v = 0; v < 3; v++)
                {
                    Vertex v1 = mesh.get_vertex(t.TV(v));
                    coords[0][v] = v1.get_x();
                    coords[1][v] = v1.get_y();
                    coords[2][v] = v1.get_z();
                }
                double a, b, c, m;
                a = (coords[1][1] - coords[1][0]) * (coords[2][2] - coords[2][0]) - (coords[2][1] - coords[2][0]) * (coords[1][2] - coords[1][0]);
                b = (coords[2][1] - coords[2][0]) * (coords[0][2] - coords[0][0]) - (coords[0][1] - coords[0][0]) * (coords[2][2] - coords[2][0]);
                c = (coords[0][1] - coords[0][0]) * (coords[1][2] - coords[1][0]) - (coords[1][1] - coords[1][0]) * (coords[0][2] - coords[0][0]);
                double tmp = a * a + b * b + c * c;
                m = sqrt(tmp);
                a = a / m;
                b = b / m;
                c = c / m;
                triangle_plane[0] = a;
                triangle_plane[1] = b;
                triangle_plane[2] = c;
                triangle_plane[3] = -1 * (a * coords[0][0] + b * coords[1][0] + c * coords[2][0]);

                for (int j = 0; j < t.vertices_num(); j++)
                {
                    if (n->indexes_vertex(t.TV(j))) // Only add the plane error to vertices indexed by the current node
                    {                               // to avoid duplicates.
                        double *plane = &(triangle_plane[0]);
                        initialQuadric[t.TV(j)] += Matrix(plane);
                    }
                }
            }
        }

        // itype v_start = n->get_v_start();
        // itype v_end = n->get_v_end();
        // itype v_range = v_end - v_start;
        //  leaf_VT local_vts(v_range, VT());
        //  n->get_VT(local_vts,mesh);
        //  compute_leaf_QEM(*n,mesh,local_vts);
    }
}

//// Another version of parallel QEM computation, need to compare the performance
// void Contraction_Simplifier::compute_initial_plane_and_QEM_parallel(PRT_Tree &tree, Mesh &mesh)
// {

// #pragma omp parallel for
//     for (int i = 0; i < tree.get_leaves_number(); i++)
//     {
//         Node_V *n = tree.get_leaf(i);
//         if (!n->indexes_vertices())
//             continue;

//         itype v_start = n->get_v_start();
//         itype v_end = n->get_v_end();
//         itype v_range = v_end - v_start;
//         leaf_VT local_vts(v_range, VT());
//         n->get_VT(local_vts, mesh);
//         for (RunIteratorPair itPair = n->make_v_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
//         {
//             RunIterator const &v_id = itPair.first;
//             VT vt = local_vts[*v_id-v_start];
//             for (auto it = vt.begin(); it != vt.end(); it++)
//             {
//                 Triangle t = mesh.get_triangle(*it);

//                 dvect triangle_plane = dvect(4, 0);
//                 double coords[3][3];
//                 for (int v = 0; v < 3; v++)
//                 {
//                     Vertex v1 = mesh.get_vertex(t.TV(v));
//                     coords[0][v] = v1.get_x();
//                     coords[1][v] = v1.get_y();
//                     coords[2][v] = v1.get_z();
//                 }
//                 double a, b, c, m;
//                 a = (coords[1][1] - coords[1][0]) * (coords[2][2] - coords[2][0]) - (coords[2][1] - coords[2][0]) * (coords[1][2] - coords[1][0]);
//                 b = (coords[2][1] - coords[2][0]) * (coords[0][2] - coords[0][0]) - (coords[0][1] - coords[0][0]) * (coords[2][2] - coords[2][0]);
//                 c = (coords[0][1] - coords[0][0]) * (coords[1][2] - coords[1][0]) - (coords[1][1] - coords[1][0]) * (coords[0][2] - coords[0][0]);
//                 double tmp = a * a + b * b + c * c;
//                 m = sqrt(tmp);
//                 a = a / m;
//                 b = b / m;
//                 c = c / m;
//                 triangle_plane[0] = a;
//                 triangle_plane[1] = b;
//                 triangle_plane[2] = c;
//                 triangle_plane[3] = -1 * (a * coords[0][0] + b * coords[1][0] + c * coords[2][0]);

//                 double *plane = &(triangle_plane[0]);
//                 initialQuadric[*v_id] += Matrix(plane);
//             }

//         }
//     }
// }

void Contraction_Simplifier::error_range(PRT_Tree &tree, Mesh &mesh, cli_parameters &cli, itype num_bin){

    initialQuadric = vector<Matrix>(mesh.get_vertices_num() + 1, Matrix(0.0));
    compute_initial_plane_and_QEM_parallel(tree, mesh);
    dvect edge_costs;
    coord_type min=INFINITY, max=-INFINITY;
    ivect hist_counter (num_bin,0);
    
    //#pragma omp parallel for  
    for (unsigned i = 0; i < tree.get_leaves_number(); i++)
        {
                Node_V *leaf = tree.get_leaf(i);
            //   coord_type n_min=INFINITY, n_max=-INFINITY;
                error_range_leaf(*leaf, mesh,edge_costs, min, max);

        }
    cout<<"min: "<<min<<" ,max: "<<max<<endl;

    //// Below is part for computing histogram
    // coord_type gap = (max-min)/num_bin;
    // for(int i=0;i< edge_costs.size();i++){

    //     itype bin_num = (edge_costs[i]-min)/gap;
    //     if(bin_num==hist_counter.size())
    //     {
    //         hist_counter[bin_num-1]++;
    //     }
    //     else
    //     hist_counter[bin_num]++;
    // }
    // cout<<"Print the count of each subrange"<<endl;
    // for(int i=0; i<hist_counter.size();i++){
    //     cout<<"bin "<<i<<": "<<hist_counter[i]<<endl;
    // }

    //// Below is for computing quartiles
    sort(edge_costs.begin(),edge_costs.end());

    auto it_second_half = edge_costs.cbegin() + edge_costs.size() / 2;
    auto it_first_half = it_second_half;
    if ((edge_costs.size() % 2) == 0) --it_first_half;

    double q1 = median(edge_costs.begin(), it_first_half);
    double q2 = median(edge_costs.begin(), edge_costs.end());
    double q3 = median(it_second_half, edge_costs.end());

    cout<<"first quartile: "<<q1<<" second quartile: "<<q2<<" third quartile: "<<q3<<endl;
    // itype size = edge_costs.size();
    // int mid = size/2;
    // double median;
    // median = size % 2 == 0 ? (edge_costs[mid] + edge_costs[mid-1])/2 : edge_costs[mid];


}


void Contraction_Simplifier::error_range_leaf(Node_V &n, Mesh &mesh, dvect& edge_costs, coord_type& min, coord_type& max){
    if (!n.indexes_vertices())
        return;

    map<ivect, coord_type> edge_map;
    ivect e;
    for (RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const &t_id = itPair.first;
        // if (mesh.is_triangle_removed(*t_id))
        // {
        //     //  cout<<"triangle removed"<<endl;
        //     continue;
        // }
        Triangle &t = mesh.get_triangle(*t_id);

        for (int i = 0; i < t.vertices_num(); i++)
        {
            t.TE(i, e);
            if(!n.indexes_vertex(e[1]))
                 continue;
            int new_vertex_pos = -1;
            double error = compute_error(e[0], e[1], mesh, new_vertex_pos);
            if(error>max)
            max = error;
            if(error<min)
            min = error;
            if (n.indexes_vertex(e[1])) // e (v1,v2) is a candidate edge if at least v2 is in n
            {
                map<ivect, coord_type>::iterator it = edge_map.find(e);
                 
                if (it == edge_map.end())
                {
   //                 cout<<e[0]<<" and "<<e[1]<<"'s error: "<< error<<endl;
                    edge_map[e] = error;
                    edge_costs.push_back(error);
                }
            }       
        }
    }
}

// void Contraction_Simplifier::error_hist_leaf(Node_V &n, Mesh &mesh,ivect& count, coord_type min, coord_type max){

// coord_type gap = (max-min)/count.size();


// }
void Contraction_Simplifier::compute_plane_and_QEM(Node_V &n, Mesh &mesh,Spatial_Subdivision &division, PRT_Tree &tree){

if (n.is_leaf())
    {
            QEM_leaf(n, mesh, tree);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            if (n.get_son(i) != NULL)
            {
                compute_plane_and_QEM(*n.get_son(i), mesh,  division, tree);
            }
        }
    }

}
void Contraction_Simplifier::QEM_leaf(Node_V &n, Mesh &mesh, PRT_Tree &tree){

        if (!n.indexes_vertices())
            return;

        for (RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const &t_id = itPair.first;

            Triangle t = mesh.get_triangle(*t_id);
            if (!n.indexes_triangle_vertices(t))
                continue;
            else
            {
                dvect triangle_plane = dvect(4, 0);
                double coords[3][3];
                for (int v = 0; v < 3; v++)
                {
                    Vertex v1 = mesh.get_vertex(t.TV(v));
                    coords[0][v] = v1.get_x();
                    coords[1][v] = v1.get_y();
                    coords[2][v] = v1.get_z();
                }
                double a, b, c, m;
                a = (coords[1][1] - coords[1][0]) * (coords[2][2] - coords[2][0]) - (coords[2][1] - coords[2][0]) * (coords[1][2] - coords[1][0]);
                b = (coords[2][1] - coords[2][0]) * (coords[0][2] - coords[0][0]) - (coords[0][1] - coords[0][0]) * (coords[2][2] - coords[2][0]);
                c = (coords[0][1] - coords[0][0]) * (coords[1][2] - coords[1][0]) - (coords[1][1] - coords[1][0]) * (coords[0][2] - coords[0][0]);
                double tmp = a * a + b * b + c * c;
                m = sqrt(tmp);
                a = a / m;
                b = b / m;
                c = c / m;
                triangle_plane[0] = a;
                triangle_plane[1] = b;
                triangle_plane[2] = c;
                triangle_plane[3] = -1 * (a * coords[0][0] + b * coords[1][0] + c * coords[2][0]);

                for (int j = 0; j < t.vertices_num(); j++)
                {
                    if (n.indexes_vertex(t.TV(j))) // Only add the plane error to vertices indexed by the current node
                    {                               // to avoid duplicates.
                        double *plane = &(triangle_plane[0]);
                        initialQuadric[t.TV(j)] += Matrix(plane);
                    }
                }
            }
        }

}


void Contraction_Simplifier::check_delaunay(PRT_Tree &tree, Mesh &mesh){


    vector<bool> tri_delaunay(mesh.get_triangles_num(),true);
    int t_num = mesh.get_triangles_num();

      #pragma omp parallel for
        for (unsigned i = 0; i < tree.get_leaves_number(); i++)
        {
            Node_V *leaf = tree.get_leaf(i);
            check_delaunay_leaf(tri_delaunay,*leaf, mesh, tree);

        }

    int count_invalid = 0;
    for(int j=0; j<tri_delaunay.size();j++){

        if(tri_delaunay[j]==false)
          count_invalid++;
    }
    cout<<"Invalid triangle num: "<<count_invalid<<endl;
    cout<<"Invalid percentage: "<<100*count_invalid/t_num<<"%"<<endl;


}

void Contraction_Simplifier::check_delaunay_leaf(vector<bool>& tri_delaunay, Node_V &n, Mesh &mesh,  PRT_Tree &tree){
    if (!n.indexes_vertices())
        return;

    itype v_start = n.get_v_start();
    itype v_end = n.get_v_end();
    itype v_range = v_end - v_start;



   for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        
        RunIterator const& t_id = itPair.first;
        Triangle& tri = mesh.get_triangle(*t_id);
        if(n.indexes_vertex(tri.minindex())){
            Vertex v0=mesh.get_vertex(tri.TV(0));
            Vertex v1=mesh.get_vertex(tri.TV(1));
            Vertex v2=mesh.get_vertex(tri.TV(2));
            coord_type x_circle;
            coord_type y_circle;
            
            double A = v1.get_x() - v0.get_x();
            double B = v1.get_y() - v0.get_y();
            double C =v2.get_x()- v0.get_x();
            double D = v2.get_y() - v0.get_y();
            double E= A*(v0.get_x()+v1.get_x())+B*(v0.get_y()+v1.get_y());
            double F= C*(v0.get_x()+v2.get_x())+D*(v0.get_y()+v2.get_y());
            double G=2.0*((A*(v2.get_y()-v1.get_y()))-(B*(v2.get_x()-v1.get_x())));
            x_circle = ( (D*E) - (B*F) ) / G;
            y_circle = ( (A*F) - (C*E) ) / G;
    
            coord_type radius=sqrt((v0.get_x()-x_circle)*(v0.get_x()-x_circle)+(v0.get_y()-y_circle)*(v0.get_y()-y_circle));
            Point center(x_circle, y_circle);
            bool is_delaunay = true;
            check_delaunay_triangle(is_delaunay,center,radius,tree.get_root(),mesh,tree.get_mesh().get_domain(),0,tree.get_subdivision());
            if(!is_delaunay)
                tri_delaunay[*t_id-1] = false;

        }

    }


}

void Contraction_Simplifier::check_delaunay_triangle(bool& is_delaunay, Point& center, coord_type radius, Node_V &n, Mesh &mesh, Box &n_dom, int level, Spatial_Subdivision &division){

    Point max(center.get_x()+radius, center.get_y()+radius);
    Point min(center.get_x()-radius, center.get_y()-radius);
    Box bounding_box(min,max);
    if(!n_dom.intersects(bounding_box)||is_delaunay==false)
        return;
    
    if (n.is_leaf())
    {
           if(n.indexes_vertices()&&Geometry_Wrapper::point_in_circle_range(center, radius,n.get_v_start(),n.get_v_end(),mesh)){
                is_delaunay = false;
                return;
           }
                
    }
    else
    {

        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(n_dom,level,i);
            int son_level = level +1;
            if (n.get_son(i) != NULL)
            {
             check_delaunay_triangle(is_delaunay, center,radius,*n.get_son(i),mesh,son_dom,son_level,division);
            }
        }
    }

}

void Contraction_Simplifier::compute_compactness(PRT_Tree &tree, Mesh &mesh){

    dvect meshCompactness(mesh.get_triangles_num(), 0);
    #pragma omp parallel for
    for (unsigned i = 0; i < tree.get_leaves_number(); i++)
    {
            Node_V *leaf = tree.get_leaf(i);

            if (!leaf->indexes_vertices())
                continue;


        for(RunIteratorPair itPair = leaf->make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
            {
                RunIterator const& t_id = itPair.first;
                Triangle& tri = mesh.get_triangle(*t_id);
                if(leaf->indexes_vertex(tri.minindex())){
                    meshCompactness[*t_id-1] = Geometry_Wrapper::triangle_compactness(*t_id,mesh);

                }

            }
        }

    coord_type max=-1, min = 1.1, avg=0;
    for(int i=0; i < meshCompactness.size(); i++){
        max = std::max(meshCompactness[i],max);
        min = std::min(meshCompactness[i],min);
        avg += meshCompactness[i];
    }

    avg = avg/meshCompactness.size();
    cout<<"Triangle compactness: Max: "<<max<<", Min: "<<min<<", Avg: "<<avg<<endl;
    //TODO: write the vector to file
}

vector<pair<coord_type, int>>& Contraction_Simplifier::get_edge_costs(bool contracted){
    if(contracted)
        return this->contracted_costs;
    else 
    {
        return this->skipped_costs;
    }

}
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

void Forman_Gradient_Simplifier::exec_local_topological_simplification(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, Spatial_Subdivision &division,
                                                                       OpType operation, int cache_size, coord_type persistence)
{
    /// reset the statistical counters
    reset_extraction_critica_counters();

    mig_cache cache = mig_cache(cache_size);

    this->local_topological_simplification(n,mesh,gradient,division,n,operation,cache,persistence);
    print_arc_nums();
  
}

void Forman_Gradient_Simplifier::local_topological_simplification(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, Spatial_Subdivision &division,
                                                                  Node_V &root, OpType operation, mig_cache &cache, coord_type persistence)
{
    if (n.is_leaf())
    {
//        cout<<n<<endl;
        /// if there are no vertices in the leaf we have nothing to do..
        if(!n.indexes_vertices())
            return;

        Timer time;
        IG local_ig;
        local_VTstar_ET local_rels;

        if(operation == TIME_VERBOSE)
            time.start();
        Forman_Gradient_Topological_Relations::get_VTstar_ET(local_rels,n,mesh,gradient);
//        cout<<" get_vtstar_et"<<endl;
        if(operation == TIME_VERBOSE)
        {
            time.stop();
            rels_gathering += time.get_elapsed_time();
            time.start();
        }
        this->get_local_MIG(n,mesh,gradient,local_ig,local_rels,root,division,cache);
        //  print_arc_nums();
       //local_ig.print_saddles();
       // Error at saddle vertices (45.25,11) (45.25,11.25) Should be removed, but becomes two vertices. 
       //cout<<"Coordinates:"<< mesh.get_vertex()
//        cout<<" get_local_mig"<<flush;
        if(operation == TIME_VERBOSE)
        {
            time.stop();
            new_paths += time.get_elapsed_time();

            /// get the MIG statistics before simplifying it!
            stats.set_local_MIG_stats(local_ig);
            stats.set_cache_stats(cache);
        }
//        cout<<" local_saddles: "<<local_ig.saddle_number()<<flush;

        if(local_ig.saddle_number() > 0)        /// I have a local MIG representation
        {
            if(operation == TIME_VERBOSE)
                time.start();
            this->local_topological_simplification_leaf(n,mesh,gradient,local_ig,local_rels, cache,root,division,persistence);
            if(operation == TIME_VERBOSE)
            {
                time.stop();
                u_paths += time.get_elapsed_time();
            }
        }
//        cout<<" simplified"<<flush;

        if(operation == TIME_VERBOSE)
            time.start();
        //aggiungo alla cache la struttura della foglia corrente
        if(mesh.get_vertices_num() > n.get_v_end())
        {
            cache.addItem(n.get_v_end()+n.get_v_start(),local_rels);
        }
        if(operation == TIME_VERBOSE)
        {
            time.stop();
            cache_t += time.get_elapsed_time();
        }

//        cout<<" cached"<<flush;

        local_ig.clear();

//        cout<<endl;
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                this->local_topological_simplification(*n.get_son(i),mesh,gradient,division,root,operation,cache,persistence);
            }
        }
    }
}

void Forman_Gradient_Simplifier::get_local_MIG(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, IG &ig, local_VTstar_ET &local_rels, Node_V &root,
                                         Spatial_Subdivision &division, mig_cache &cache)
{
    for(leaf_ET::iterator it_e=local_rels.begin_ETs(); it_e!=local_rels.end_ETs(); ++it_e)
    {
        //each edge is visited only once.. precisely in the leaf block containing the extreme with maximum index
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

            /// visit the descending 1 cells ///
            /// we visit the two extremes with two distinct visits!
            explore_desc1cell_local_mig(n, make_pair(it_e->first[0],v1_pos), mesh, gradient, ig, local_rels, node, root, division, cache);
            explore_desc1cell_local_mig(n, make_pair(it_e->first[1],v2_pos), mesh, gradient, ig, local_rels, node, root, division, cache);

            /// visit the ascending 1 cells ///
            pair<itype,short> coppia = make_pair(it_e->second.first,e_pos_first);
            explore_asc1cell_local_mig(n, coppia, node, it_e->second.first, local_rels.get_ETs(), mesh, gradient, ig, root, division, cache);
            if(it_e->second.second != -1)
            {
                pair<itype,short> coppia2 = make_pair(it_e->second.second,e_pos_second);
                explore_asc1cell_local_mig(n, coppia2, node, it_e->second.second, local_rels.get_ETs(), mesh, gradient, ig, root, division, cache);
            }
        }
    }
}

void Forman_Gradient_Simplifier::explore_desc1cell_local_mig(Node_V &n, const pair<itype, short> &v_pair, Mesh &mesh, Forman_Gradient &gradient, IG &ig,
                                                       local_VTstar_ET &local_rels, iNode *saddle_node, Node_V &root, Spatial_Subdivision &division, mig_cache &cache)
{
    /// v_pair contains the vertex index from which start the navigation
    /// and the position of the saddle vertex in the first triangle of the ET(saddle)
    queue<pair<itype, short> > coda;
    coda.push(v_pair);
    pair<itype, short> current_pair;

    itype last_v = v_pair.first; /// last vertex visited

    while(!coda.empty())
    {
        current_pair = coda.front();
        coda.pop();

        /// we need to expand fully the navigation here
        itype vtstar = Forman_Gradient_Topological_Relations::get_VTstar(n,current_pair.first,local_rels.get_VTstars(),cache.get_vtstar_cache(),root,division,mesh,gradient);
        if(!gradient.is_vertex_critical(current_pair.first,vtstar,mesh)) //*it is paired
        {
            Triangle &t = mesh.get_triangle(vtstar);
            short v1i = gradient.convert_compressed_to_expand(vtstar).get_vertex_pair(t.vertex_index(current_pair.first));
            itype v1 = t.TV(v1i);
            if(v1 != current_pair.first) /// this should not happen
            {
                coda.push(make_pair(v1,current_pair.second));
                last_v = current_pair.first;
            }
        }
        else
        {
            pair_minimum_saddle(saddle_node,current_pair.second,current_pair.first,last_v,ig);
        }
    }
}

void Forman_Gradient_Simplifier::explore_asc1cell_local_mig(Node_V &n, const pair<itype, short> &tf, iNode *saddle_node, itype first_t_id, leaf_ET &local_ef,
                                                      Mesh &mesh, Forman_Gradient &gradient, IG &ig, Node_V &root, Spatial_Subdivision &division, mig_cache &cache)
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

        /// we need to expand fully the navigation here
        short f1 = gradient.convert_compressed_to_expand(current_pair.first).get_face_pair();

        if(f1 != current_pair.second)
        {
            if(f1 > -1)
            {
                t.TE(f1,edge);
                pair<itype,itype> et = Forman_Gradient_Topological_Relations::get_ET(n,edge,local_ef,cache.get_et_cache(),root,division,mesh);

                itype next = (et.first == current_pair.first) ? et.second : et.first;
                if(next != -1)
                {
                    f1 = mesh.get_triangle(next).edge_index(edge);
                    coda.push(make_pair(next,f1));
                    last_t = current_pair.first;
                }

            }
            else if (gradient.is_triangle_critical(current_pair.first))
            {
         
                pair_saddle_maximum(saddle_node,first_t_id,current_pair.first,last_t,ig);
            }
        }
    }
}

///when we call this function we must have built the local relations and the local MIG
void Forman_Gradient_Simplifier::local_topological_simplification_leaf(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, IG& ig, local_VTstar_ET &local_rels, mig_cache &cache,
                                                                 Node_V& root, Spatial_Subdivision &division, coord_type persistence)
{
//    cout<<n<<endl;
    priority_arcs_queue queue;
    this->build_persistence_queue(queue,ig,mesh,persistence); /// add a subset of the arcs of the local MIG
//    cout<<"built persistence_queue"<<endl;

    if(max_priority_queue_size < queue.size())
        max_priority_queue_size = queue.size();

//    cout<<"start simplification"<<endl;
    simplify(queue,n,mesh,gradient,ig,local_rels, cache,root,division,persistence);
}

void Forman_Gradient_Simplifier::simplify(priority_arcs_queue &queue, Node_V &n, Mesh &mesh, Forman_Gradient &gradient, IG& ig, local_VTstar_ET &local_rels,
                                    mig_cache &cache, Node_V& root, Spatial_Subdivision &division, coord_type persistence)
{
    Topo_Sempl sempl;
    iNode* saddle=NULL;
    nNode* extrema=NULL;
    int lvl0=0,lvl1=0,lvl2=0;
    while(!queue.empty())
    {
        sempl = queue.top();
        queue.pop();

        if(sempl.arc->getLabel() != 1 )
        {
            lvl2++;
            if(sempl.arc->getLabel() == -1)
                delete sempl.arc;
            continue;
        }

        if(sempl.lvl == 0)
        {
            lvl0++;
            saddle  = (iNode*)sempl.arc->getNode_j();
//FOR DEBUG:
        itype sad_vertex= saddle->get_critical_index();
       // itype sad_vertex = node->get_critical_index();
           // cout<<"Coordinates of saddle: "<<mesh.get_vertex(sad_vertex).get_x()<<", "<<mesh.get_vertex(sad_vertex).get_y()<<endl;
            extrema = (nNode*)sempl.arc->getNode_i();

            if(saddle->getArcs(true).size() != 2) continue;
           // cout<<"Persistence value:"<<sempl.val<<endl;
            contraction(extrema, saddle, queue, ig, mesh, gradient, local_rels, cache, n, root, division, persistence);

            refined_topo++;
        }
        else
        {
            lvl1++;
            saddle  = (iNode*)sempl.arc->getNode_i();

           //  itype sad_vertex= saddle->get_critical_index();
       // itype sad_vertex = node->get_critical_index();
           // cout<<"Coordinates of saddle: "<<mesh.get_vertex(sad_vertex).get_x()<<", "<<mesh.get_vertex(sad_vertex).get_y()<<endl;


            extrema = (nNode*)sempl.arc->getNode_j();
           // itype t_id= sempl.arc->getSimplexj();
            //cout<<"Ending triangle "<<t_id<<endl;
            if(saddle->getArcs(false).size() != 2) continue;
            removal(extrema,saddle, queue, ig, mesh, gradient, local_rels, cache, n, root, division, persistence);

            refined_topo++;
        }

        if(max_priority_queue_size < queue.size())
            max_priority_queue_size = queue.size();
    }

//cout<<"Level 0:"<<lvl0<<", Level 1:"<<lvl1<<", Level 2:"<<lvl2<<endl;
}

/// involves a minimum and a saddle
void Forman_Gradient_Simplifier::contraction(nNode *extrema, iNode *saddle, priority_arcs_queue &q, IG &ig, Mesh &mesh, Forman_Gradient &gradient,
                                       local_VTstar_ET &local_rels, mig_cache &cache, Node_V &n, Node_V& root, Spatial_Subdivision &division,
                                       coord_type persistence)
{
    vector<Arc*> arcs = saddle->get_vector_Arcs(true);

    pair<itype,itype> critical_edge_tetra=saddle->get_edge_id();
    ivect critical_edge;

    Triangle &first = mesh.get_triangle(critical_edge_tetra.first);

    if(critical_edge_tetra.second<0)
    {
        first.TE(-critical_edge_tetra.second-1,critical_edge);
    }
    else
    {
        Triangle &second = mesh.get_triangle(critical_edge_tetra.second);
        for(int i=0; i<3; i++)
        {
            if(!(second.has_vertex(first.TV(i))))
            {
                first.TE(i,critical_edge);
                break;
            }
        }
    }

    nNode* other_extrema;
    itype vertex,next_vertex;
    itype ending_path_simplex;
    //SETUP per dopo
    if(arcs[0]->getNode_i() == extrema){
        other_extrema = (nNode*)arcs[1]->getNode_i();
        next_vertex = first.TV(arcs[0]->getSimplexj());
        vertex = first.TV(arcs[1]->getSimplexj());
        ending_path_simplex = arcs[1]->getSimplexi();
    }
    else{
        other_extrema = (nNode*)arcs[0]->getNode_i();
        next_vertex = first.TV(arcs[1]->getSimplexj());
        vertex = first.TV(arcs[0]->getSimplexj());
        ending_path_simplex = arcs[0]->getSimplexi();
    }

    /// ----------------- MODIFICHE SUL GRADIENTE ----------------- ///
    contraction_update_gradient(vertex, extrema->get_critical_index(), next_vertex, critical_edge, mesh, gradient, local_rels, cache, n, root, division);

    /// ----------------- MODIFICHE SUL MIG ----------------- ///
    arcs[0]->setLabel(-1);
    arcs[1]->setLabel(-1);

    /// remove maxima-saddle arcs
    remove_saddle_arcs(saddle,false,ig);

    /// remove extrema arcs
    remove_extreme_arcs(extrema, saddle, other_extrema, ending_path_simplex, true, ig, q, mesh, persistence);

    /// remove minima-saddle arcs
    remove_saddle_arcs(saddle,true,ig);

    /// remove the saddle and the minimum

   // cout<<"CONTRACTION SADDLE: "<<*saddle<<endl;
    ig.remove_saddle(saddle->get_edge_id(),saddle);
    ig.remove_minimum(extrema->get_critical_index(),extrema);
//    int a; cin>>a;

//    ivect edge = {11376, 11377};
//    ET et = {22741, 22728};
//    cout<<"IS THE TARGET EDGE CRITICAL? "<<gradient.is_edge_critical(edge,et,mesh)<<endl;
}

///involves a saddle and a maximum
void Forman_Gradient_Simplifier::removal(nNode *extrema, iNode *saddle, priority_arcs_queue &q, IG &ig, Mesh &mesh, Forman_Gradient &gradient,
                                   local_VTstar_ET &local_rels, mig_cache &cache, Node_V &n, Node_V &root, Spatial_Subdivision &division, coord_type persistence)
{
    vector<Arc*> arcs = saddle->get_vector_Arcs(false);
    nNode* other_extrema;
    itype triangle;
    itype ending_path_simplex;

    if(arcs[0]->getNode_j() == extrema){
        other_extrema = (nNode*)arcs[1]->getNode_j();
        triangle = arcs[0]->getSimplexi();
        ending_path_simplex = arcs[1]->getSimplexj();
    }
    else{
        other_extrema = (nNode*)arcs[0]->getNode_j();
        triangle = arcs[1]->getSimplexi();
        ending_path_simplex = arcs[0]->getSimplexj();
    }

    /// ----------------- MODIFICHE SUL GRADIENTE ----------------- ///
    removal_update_gradient(triangle, saddle, mesh, gradient, local_rels.get_ETs(), cache, n, root, division);

    /// ----------------- MODIFICHE SUL MIG ----------------- ///
    arcs[0]->setLabel(-1);
    arcs[1]->setLabel(-1);

    /// remove minima-saddle arcs
    remove_saddle_arcs(saddle,true,ig);

    /// remove extrema arcs
    remove_extreme_arcs(extrema, saddle, other_extrema, ending_path_simplex, false, ig, q, mesh, persistence);

    /// remove maxima-saddle arcs
    remove_saddle_arcs(saddle,false,ig);

    /// remove the saddle and the maximum
   // cout<<"REMOVAL SADDLE: "<<*saddle<<endl;
    ig.remove_saddle(saddle->get_edge_id(),saddle);
    ig.remove_maximum(extrema->get_critical_index(),extrema);
}

void Forman_Gradient_Simplifier::contraction_update_gradient(itype vertex, itype ex_minimum, itype next_vertex, const ivect &critical_edge, Mesh &mesh,
                                                             Forman_Gradient &gradient, local_VTstar_ET &local_rels, mig_cache &cache,
                                                             Node_V &n, Node_V &root, Spatial_Subdivision &division)
{
    ivect old_edge = critical_edge;
    ivect edge;
    pair<itype,itype> old_ef = Forman_Gradient_Topological_Relations::get_ET(n,old_edge,local_rels.get_ETs(),cache.get_et_cache(),root,division,mesh);
    while(next_vertex != ex_minimum)
    {
        //trovo il nuovo edge;
        itype vtstar = Forman_Gradient_Topological_Relations::get_VTstar(n,next_vertex,local_rels.get_VTstars(),cache.get_vtstar_cache(),root,division,mesh,gradient);
        Triangle &t = mesh.get_triangle(vtstar);
        short v1i = gradient.convert_compressed_to_expand(vtstar).get_vertex_pair(t.vertex_index(next_vertex));

        /// no checking vli != -1 because it must be paired
        itype v1 = t.TV(v1i);

        edge = { min(next_vertex,v1), max(next_vertex,v1) };
        vertex=next_vertex;
        next_vertex = v1;

        //azzero la sua adiacenza;
        pair<itype,itype> ef = Forman_Gradient_Topological_Relations::get_ET(n,edge,local_rels.get_ETs(),cache.get_et_cache(),root,division,mesh);
        gradient.free_VE(vertex,next_vertex,ef,mesh);

        //accoppio il vecchio edge;
        if(old_edge[0] == vertex)
        {
            gradient.set_VE(old_edge[0], old_edge[1], old_ef, mesh);
            Forman_Gradient_Topological_Relations::set_VTstar(n,old_edge[0],old_ef.first,local_rels.get_VTstars(),cache.get_vtstar_cache(),root,division);
        }
        else{
            gradient.set_VE(old_edge[1], old_edge[0], old_ef, mesh);
            Forman_Gradient_Topological_Relations::set_VTstar(n,old_edge[1],old_ef.first,local_rels.get_VTstars(),cache.get_vtstar_cache(),root,division);
        }
        /// save for the next step
        old_edge = edge;
        old_ef = ef;
    }
    //assert(next_vertex == ex_minimum);
    gradient.set_VE(next_vertex, vertex, old_ef, mesh);
    Forman_Gradient_Topological_Relations::set_VTstar(n,next_vertex,old_ef.first,local_rels.get_VTstars(),cache.get_vtstar_cache(),root,division);
}

void Forman_Gradient_Simplifier::removal_update_gradient(itype triangle, iNode *saddle, Mesh &mesh, Forman_Gradient &gradient, leaf_ET &local_ef, mig_cache &cache,
                                                         Node_V &n, Node_V& root, Spatial_Subdivision &division)
{
    ivect old_edge;
    ivect edge;
    pair<itype,itype> twotriangles= saddle->get_edge_id();
    Triangle &tri = mesh.get_triangle(twotriangles.first);

    if(twotriangles.second >= 0)
    {
        ivect e;
        Triangle &tri2 = mesh.get_triangle(twotriangles.second);
        for(int j=0; j<tri2.vertices_num();j++)
        {
            tri.TE(j,e);
            if(tri2.has_edge(e))
            {
                old_edge = e;
                break;
            }
        }
    }
    else
    {
        tri.TE(-twotriangles.second -1, old_edge);
    }

    /// ----------------- GRADIENTE UPDATES ----------------- ///
    itype next_triangle=triangle;

    while(!gradient.is_triangle_critical(next_triangle))
    {
        Triangle &tri = mesh.get_triangle(next_triangle);
        int i = get_paired_edge(gradient,tri,next_triangle,edge);

        triangle=next_triangle;

        pair<itype,itype> et = Forman_Gradient_Topological_Relations::get_ET(n,edge,local_ef,cache.get_et_cache(),root,division,mesh);
        
        next_triangle = (et.first == triangle) ? et.second : et.first;
        int vertex=mesh.get_triangle(triangle).TV(i);

        gradient.free_ET(i,triangle);

        i = tri.edge_index(old_edge);
        vertex=mesh.get_triangle(triangle).TV(i);
        gradient.set_ET(i,triangle);

        old_edge = edge;

      
  
    }

    Triangle &tri_next = mesh.get_triangle(next_triangle);    
    int i = tri_next.edge_index(old_edge);

    gradient.set_ET(i,next_triangle);
}

/// push into the priority queue only those arc below the avg
void Forman_Gradient_Simplifier::build_persistence_queue(priority_arcs_queue &q, IG &ig, Mesh &mesh, coord_type persistence)
{
    for(int i=0; i<2; i++)
    {
        set<Arc*>& arcs = ig.getLevelArcs(i);
        for(set<Arc*>::iterator it = arcs.begin(); it != arcs.end(); it++)
        {
            if((*it)->getLabel() == 1)
            {
                itype index_i=(*it)->getNode_i()->get_critical_index();
                Vertex &v1 = mesh.get_vertex(index_i);
                itype index_j;
                coord_type val;
                Topo_Sempl ts;
                
                if(i==1)
                {
                    pair<itype,itype> critical_edge_tetra=((iNode*)(*it)->getNode_i())->get_edge_id();
                    ivect critical_edge;
                    Triangle &first = mesh.get_triangle(critical_edge_tetra.first);
                    if(critical_edge_tetra.second<0)
                    {
                        first.TE(-critical_edge_tetra.second-1,critical_edge);
                    }
                    else
                    {
                        Triangle &second = mesh.get_triangle(critical_edge_tetra.second);
                        for(int i=0; i<3; i++)
                        {
                            if(!(second.has_vertex(first.TV(i))))
                            {
                                first.TE(i,critical_edge);
                                break;
                            }
                        }
                    }

                     ivect filt_ex;
                        Triangle t=mesh.get_triangle((*it)->getNode_j()->get_critical_index());

                        for (int i=0;i<3;i++)
                            filt_ex.push_back(filtration[t.TV(i)-1]);
                    sort(filt_ex.begin(), filt_ex.end(), greater<int>()); 
                   // index_j=get_max_elevation_vertex(mesh.get_triangle((*it)->getNode_j()->get_critical_index()));
                   Vertex &v2=mesh.get_vertex(index_j);
                    val = fabs(v1.get_z() - v2.get_z());
                     ts= Topo_Sempl(*it, val, i,filtration[critical_edge[0]-1],filtration[critical_edge[1]-1],filt_ex);
                }
                else{
                    index_j=(*it)->getNode_j()->get_critical_index();

                    pair<itype,itype> critical_edge_tetra=((iNode*)(*it)->getNode_j())->get_edge_id();
                    ivect critical_edge;
                    Triangle &first = mesh.get_triangle(critical_edge_tetra.first);
                    if(critical_edge_tetra.second<0)
                    {
                        first.TE(-critical_edge_tetra.second-1,critical_edge);
                    }
                    else
                    {
                        Triangle &second = mesh.get_triangle(critical_edge_tetra.second);
                        for(int i=0; i<3; i++)
                        {
                            if(!(second.has_vertex(first.TV(i))))
                            {
                                first.TE(i,critical_edge);
                                break;
                            }
                        }
                    }
            
                    Vertex &v2 = mesh.get_vertex(index_j);
                    val = fabs(v1.get_z() - v2.get_z());
                    ivect filt_ex;
                    filt_ex.push_back(filtration[index_i-1]);
                    ts= Topo_Sempl(*it, val, i,filtration[critical_edge[0]-1],filtration[critical_edge[1]-1],filt_ex);


                }
                
         
                if(val <= persistence) /// NEW <= instead of <
                {
                   // cout<<"SADDLE is "<< *saddle<<endl;
                   // cout<<"persistence value is:"<<val<<endl;
                    
                    q.push(ts);
                }
            }
        }
    }
}


void Forman_Gradient_Simplifier::remove_extreme_arcs(nNode* extrema, iNode* saddle, nNode *other_extrema, itype ending_path_simplex,
                                                     bool is_minimum, IG &ig, priority_arcs_queue &q, Mesh &mesh, coord_type persistence)
{
    iNode* node_saddle1=NULL;
    itype starting_path_simplex;
    coord_type val;

    set<Arc*>::iterator it=extrema->begin();
    while(it!=extrema->end())
    {
        if(is_minimum)
        {
            node_saddle1 = ((iNode*)(*it)->getNode_j());
            starting_path_simplex = (*it)->getSimplexj();
        }
        else
        {
            node_saddle1 = ((iNode*)(*it)->getNode_i());
            starting_path_simplex = (*it)->getSimplexi();
        }

        (*it)->setLabel(-1);

        ig.removeArc(!is_minimum,*it);
        node_saddle1->removeArc(is_minimum, *it);
        extrema->removeArc(it);

        if(node_saddle1 != saddle)
        {
            Arc* existing_arc = ig.already_connected(other_extrema,node_saddle1);
            if(existing_arc==NULL)
            {
                Arc* arco = NULL;

                if(is_minimum)
                    arco = ig.addArc(other_extrema, ending_path_simplex, node_saddle1, starting_path_simplex, 0);
                else
                    arco = ig.addArc(node_saddle1, starting_path_simplex, other_extrema, ending_path_simplex, 1);

                //cout<<"addArc: "<<*arco<<endl;

                if(arco->getLabel() == 1)
                {
                    int index_ex;
                    ivect critical_edge;
                    ivect filt_ex;
                    /// check how to compute val..
                    if(is_minimum)
                       {
                        index_ex=arco->getNode_i()->get_critical_index();
                        filt_ex.push_back(filtration[index_ex-1]);
                        itype index_j=arco->getNode_j()->get_critical_index();
                        pair<itype,itype> critical_edge_tetra=((iNode*)arco->getNode_j())->get_edge_id();
                        
                        Triangle &first = mesh.get_triangle(critical_edge_tetra.first);
                        if(critical_edge_tetra.second<0)
                        {
                            first.TE(-critical_edge_tetra.second-1,critical_edge);
                        }
                        else
                        {
                            Triangle &second = mesh.get_triangle(critical_edge_tetra.second);
                            for(int i=0; i<3; i++)
                            {
                                if(!(second.has_vertex(first.TV(i))))
                                {
                                    first.TE(i,critical_edge);
                                    break;
                                }
                            }
                        }

                           index_j=arco->getNode_j()->get_critical_index();
                            val = fabs(mesh.get_vertex(index_ex).get_z() -
                                   mesh.get_vertex(index_j).get_z());
                       }
                    else
                      { 
                         itype index_i=arco->getNode_i()->get_critical_index();

                        pair<itype,itype> critical_edge_tetra=((iNode*)arco->getNode_i())->get_edge_id();
                        
                        Triangle &first = mesh.get_triangle(critical_edge_tetra.first);
                        if(critical_edge_tetra.second<0)
                        {
                            first.TE(-critical_edge_tetra.second-1,critical_edge);
                        }
                        else
                        {
                            Triangle &second = mesh.get_triangle(critical_edge_tetra.second);
                            for(int i=0; i<3; i++)
                            {
                                if(!(second.has_vertex(first.TV(i))))
                                {
                                    first.TE(i,critical_edge);
                                    break;
                                }
                            }
                        }
                            Triangle t=mesh.get_triangle(arco->getNode_j()->get_critical_index());
                          index_ex=get_max_elevation_vertex(t);
                           val = fabs(mesh.get_vertex(index_i).get_z() -  mesh.get_vertex(index_ex).get_z());
                      
                          for (int i=0;i<3;i++)
                            filt_ex.push_back(filtration[t.TV(i)-1]);
                        sort(filt_ex.begin(), filt_ex.end(), greater<int>()); 
                      
                      }
                    /// we push in queue if:
                    /// we simplify using a percentage based critarion
                    /// or if the persistence is below the average
                    //---------COMMENTED SINCE IT IS NOT INCLUDED IN IA VERSION---------
                    if(val <= persistence) /// NEW <= instead of <
                    {
                        /// this must be enable if we simplify only topologically!! (the same holds in the removal function)
                        Topo_Sempl ts = Topo_Sempl(arco, val, !is_minimum,filtration[critical_edge[0]-1],filtration[critical_edge[1]-1],filt_ex);
                        q.push(ts);
                    }
                }
            }
            else
            {
                existing_arc->setLabel(2);
            }

        }

        /// set again to begin
        it=extrema->begin();
    }
}

coord_type Forman_Gradient_Simplifier::get_average_persistence_value(IG &ig, Mesh &mesh)
{
    coord_type avg = 0.0;

    for(int i=0; i<2; i++)
    {
        set<Arc*>& arcs = ig.getLevelArcs(i);
        for(set<Arc*>::iterator it = arcs.begin(); it != arcs.end(); it++)
        {
            if((*it)->getLabel() == 1)
            {
                Vertex &v1 = mesh.get_vertex((*it)->getNode_i()->get_critical_index());

                if(i==1)
                {
                    Vertex &v2 = mesh.get_vertex(get_max_elevation_vertex(mesh.get_triangle((*it)->getNode_j()->get_critical_index())));
                    avg += fabs(v1.get_z() - v2.get_z());
                }
                else{
                    Vertex &v2 = mesh.get_vertex((*it)->getNode_j()->get_critical_index());
                    avg += fabs(v1.get_z() - v2.get_z());
                }

            }
        }
    }

    avg /= ig.get_arcs_number();

    /// for debug: check how many are below the average
    cerr<<"total arc number: "<<ig.get_arcs_number()<<endl;
    int arc_below_counter = 0;
    for(int i=0; i<2; i++)
    {
        set<Arc*>& arcs = ig.getLevelArcs(i);
        for(set<Arc*>::iterator it = arcs.begin(); it != arcs.end(); it++)
        {
            if((*it)->getLabel() == 1)
            {
                Vertex &v1 = mesh.get_vertex((*it)->getNode_i()->get_critical_index());

                if(i==1)
                {
                    Vertex &v2 = mesh.get_vertex(get_max_elevation_vertex(mesh.get_triangle((*it)->getNode_j()->get_critical_index())));
                    if(avg > fabs(v1.get_z() - v2.get_z()))
                        arc_below_counter++;
                }
                else{
                    Vertex &v2 = mesh.get_vertex((*it)->getNode_j()->get_critical_index());
                    if(avg > fabs(v1.get_z() - v2.get_z()))
                        arc_below_counter++;
                }
            }
        }
    }
    cerr<<"arcs below the average: "<<arc_below_counter<<endl;
    cerr<<"average_persistence: "<<avg<<endl;

    return avg;
}

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

#include "forman_gradient_topological_relations.h"
#include "utilities/sorting.h"

pair<itype,itype> Forman_Gradient_Topological_Relations::get_ET(Node_V &n, ivect &e, leaf_ET &local_et, et_cache &cache, Node_V &root, Spatial_Subdivision &division, Mesh &mesh)
{
    leaf_ET::iterator it;
    if(n.indexes_vertex(e[1]))
    {
        it = local_et.find(e);
        if(it != local_et.end())
        {
            return it->second;
        }
    }
    else
    {
        //the face is completely outside the current leaf..
        //then I search in cache for the ft of its maximum vertex
        itype key = -1;
        Node_V& other_n = Dangling_Paths_Handler::find_leaf(e[1],root,division,key);

        et_cache::mapIt iter = cache.find(key);
        if(iter == cache.end())
        {
            //unlucky.. I have to compute the vt lists and save it in cache
            leaf_ET inside;
            other_n.get_ET(inside,mesh);
            cache.insert(key,inside);
            it = inside.find(e);

            if(it != inside.end())
            {
                return it->second;
            }
        }
        else
        {
            it = iter->second.find(e);
            if(it != iter->second.end())
            {
                return it->second;
            }
        }
    }

    //// debug purpose only --> executed only if something goes wrong
    pair<itype,itype> ret;
    cout<<"[get_FT] non dovrei mai arrivarci"<<endl;
    cout<<"failing edge: "<<e[0]<<" "<<e[1]<<endl;
    cout<<"current leaf: "<<n<<endl;
    if(n.indexes_vertex(e[1]))
    {
        int counter = 0;
        cout<<"indexed"<<endl;
        cout<<"brute force algorithm"<<endl;
        for(itype i=1; i<=mesh.get_triangles_num(); i++)
        {
            Triangle &tri = mesh.get_triangle(i);
            if(tri.has_edge(e))
            {
                cout<<i<<" is incident into the edge "<<endl;
                counter++;
                if(counter == 2)
                    break;
            }
        }
        cout<<"end of brute force"<<endl;

    }
    else
    {
        itype key = -1;
        Node_V& other_n = Dangling_Paths_Handler::find_leaf(e[1],root,division,key);

        cout<<"key: "<<key<<endl;

        et_cache::mapIt iter = cache.find(key);
        if(iter == cache.end())
        {
            cout<<"not cached"<<endl;
            //unlucky.. I have to compute the vt lists and save it in cache
            leaf_ET inside;
            other_n.get_ET(inside,mesh);
            cache.insert(key,inside);
            it = inside.find(e);

            if(it != inside.end())
            {
                cout<<"it's here"<<endl;
            }
            else
                cout<<"nothing found.. it seems that the simplex does not exists.."<<endl;
        }
        else
        {
            cout<<"cached"<<endl;
            it = iter->second.find(e);
            if(it != iter->second.end())
            {
                cout<<"it's here"<<endl;
            }
            else
            {
                cout<<"nothing found.. it seems that the simplex does not exists.."<<endl;
                cout<<other_n<<endl;
                leaf_ET inside;
                other_n.get_ET(inside,mesh);
                for(leaf_ET::iterator it=inside.begin(); it!=inside.end(); ++it)
                {
                    cout<<it->first[0]<<" "<<it->first[1]<<" -> "<<it->second.first<<" "<<it->second.second<<endl;
                }
            }
        }
    }
    int a; cin>>a;
    return ret;
    /////////////// --------------------------------------- /////////////
}

pair<itype,itype> Forman_Gradient_Topological_Relations::get_ET(Node_T &n, itype v_start, itype v_end, ivect &e, leaf_ET &local_et, et_cache &cache,
                                                                Node_T &root, Spatial_Subdivision &division, Mesh &mesh)
{
    leaf_ET::iterator it;
    if(n.indexes_vertex(e[1],v_start,v_end))
    {
        it = local_et.find(e);
        if(it != local_et.end())
        {
            return it->second;
        }
    }
    else
    {
        //the face is completely outside the current leaf..
        //then I search in cache for the ft of its maximum vertex
        itype other_v_start, other_v_end;
        Node_T& other_n = Dangling_Paths_Handler::find_leaf(mesh.get_vertex(e[1]),root,mesh.get_domain(),0,mesh,division,other_v_start,other_v_end);
        itype key = other_v_start+other_v_end;
//        Node_T& other_n = Dangling_Paths_Handler::find_leaf(e[1],root,division,key);

        et_cache::mapIt iter = cache.find(key);
        if(iter == cache.end())
        {
            //unlucky.. I have to compute the vt lists and save it in cache
            leaf_ET inside;
            other_n.get_ET(inside,other_v_start, other_v_end,mesh);
            cache.insert(key,inside);
            it = inside.find(e);

            if(it != inside.end())
            {
                return it->second;
            }
        }
        else
        {
            it = iter->second.find(e);
            if(it != iter->second.end())
            {
                return it->second;
            }
        }
    }
}

iset& Forman_Gradient_Topological_Relations::get_VV(Node_V &n, itype v, leaf_VV &all_vv, asc3_cache &cache, Node_V &root, Spatial_Subdivision &division, Mesh &mesh)
{
    itype v_pos = -1;
    if(n.indexes_vertex(v))
    {
        v_pos = v - n.get_v_start();
        return all_vv[v_pos];
    }

    //then I search in cache for the vt of this vertex
    itype key = -1;
    Node_V& other_n = Dangling_Paths_Handler::find_leaf(v,root,division,key);
    v_pos = v - other_n.get_v_start();

    vv_cache::mapIt it = cache.find_vv_lists(key);
    if(it == cache.end_vv_cache())
    {
        //unlucky.. I have to compute the vt lists and save it in cache
        leaf_VV inside;
        other_n.get_VV(inside,mesh/*,forman_gradient*/);
        vv_cache::mapIt ret = cache.add_vv_lists(key,inside);
        return ret->second[v_pos];
    }
    else
    {
        return it->second[v_pos];
    }
}

itype Forman_Gradient_Topological_Relations::get_VTstar(Node_V &n, itype v, leaf_VTstar &all_vtstar, vtstar_cache &cache, Node_V &root,
                                                          Spatial_Subdivision &division, Mesh &mesh, Forman_Gradient &gradient)
{
    itype v_pos = -1;
    if(n.indexes_vertex(v))
    {
        v_pos = v - n.get_v_start();
        return all_vtstar[v_pos];
    }

    //then I search in cache for the vt of this vertex
    itype key = -1;
    Node_V& other_n = Dangling_Paths_Handler::find_leaf(v,root,division,key);
    v_pos = v - other_n.get_v_start();

    vtstar_cache::mapIt it = cache.find(key);
    if(it == cache.end())
    {
        //unlucky.. I have to compute the vt lists and save it in cache
        leaf_VTstar inside;
        Forman_Gradient_Topological_Relations::get_VTstar(inside,other_n,mesh,gradient);
        vtstar_cache::mapIt ret = cache.insert(key,inside);
        return ret->second[v_pos];
    }
    else
    {
        return it->second[v_pos];
    }
}

itype Forman_Gradient_Topological_Relations::get_VTstar(Node_T &n, itype v_start, itype v_end, itype v, leaf_VTstar &all_vtstar, vtstar_cache &cache, Node_T &root,
                                                          Spatial_Subdivision &division, Mesh &mesh, Forman_Gradient &gradient)
{
    itype v_pos = -1;
    if(n.indexes_vertex(v_start,v_end,v))
    {
        v_pos = v - v_start;
        return all_vtstar[v_pos];
    }

    //then I search in cache for the vt of this vertex
    itype other_v_start, other_v_end;
    Node_T& other_n = Dangling_Paths_Handler::find_leaf(mesh.get_vertex(v),root,mesh.get_domain(),0,mesh,division,other_v_start,other_v_end);
    v_pos = v - other_v_start;
    itype key = other_v_start+other_v_end;

    vtstar_cache::mapIt it = cache.find(key);
    if(it == cache.end())
    {
        //unlucky.. I have to compute the vt lists and save it in cache
        leaf_VTstar inside;
        Forman_Gradient_Topological_Relations::get_VTstar(inside,other_n,other_v_start,other_v_end,mesh,gradient);
        vtstar_cache::mapIt ret = cache.insert(key,inside);
        return ret->second[v_pos];
    }
    else
    {
        return it->second[v_pos];
    }
}


/// NOTA: I update a VT* of a vertex only if it is locally indexed by the leaf or if the leaf that index it is cached..
/// I do not need to rebuild.. because, once needed, it will be consistent
void Forman_Gradient_Topological_Relations::set_VTstar(Node_V &n, itype v, itype new_vt_star, leaf_VTstar &all_vtstar, vtstar_cache &cache,
                                                       Node_V &root, Spatial_Subdivision &division)
{
    itype v_pos = -1;
    if(n.indexes_vertex(v))
    {
        v_pos = v - n.get_v_start();
        all_vtstar[v_pos] = new_vt_star;
    }

    //then I search in cache for the vt of this vertex
    itype key = -1;
    Node_V& other_n = Dangling_Paths_Handler::find_leaf(v,root,division,key);
    v_pos = v - other_n.get_v_start();

    vtstar_cache::mapIt it = cache.find(key);
    if(it != cache.end())
    {
        it->second[v_pos] = new_vt_star;
    }
}

void Forman_Gradient_Topological_Relations::get_VTstar(leaf_VTstar &vtstars, Node_V &n, Mesh& mesh, Forman_Gradient &gradient)
{
    vtstars.assign(n.get_v_end()-n.get_v_start(),-1);

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);

        for(int v=0; v<t.vertices_num(); v++)
        {
            if(n.indexes_vertex(t.TV(v)))
                Forman_Gradient_Topological_Relations::check_VTstar(t.TV(v),*t_id,t,vtstars,n.get_v_start(),mesh,gradient);
        }
    }
}

void Forman_Gradient_Topological_Relations::get_VTstar(leaf_VTstar &vtstars, Node_T &n, itype v_start, itype v_end, Mesh& mesh, Forman_Gradient &gradient)
{
    vtstars.assign(v_end-v_start,-1);

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);

        for(int v=0; v<t.vertices_num(); v++)
        {
            if(n.indexes_vertex(v_start,v_end,t.TV(v)))
                Forman_Gradient_Topological_Relations::check_VTstar(t.TV(v),*t_id,t,vtstars,v_start,mesh,gradient);
        }
    }
}

void Forman_Gradient_Topological_Relations::get_VTstar_ETstar(desc1rels &all_rels, Node_V &n, Mesh& mesh, Forman_Gradient &gradient)
{
    all_rels.init((n.get_v_end()-n.get_v_start()));
    ivect e;
    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);
        for(int v=0; v<t.vertices_num(); v++)
        {
            if(n.indexes_vertex(t.TV(v)))
                Forman_Gradient_Topological_Relations::check_VTstar(t.TV(v),*t_id,t,all_rels.vtstars,n.get_v_start(),mesh,gradient);
            t.TE(v,e);
            if(n.indexes_vertex(e[1]))
                Forman_Gradient_Topological_Relations::check_ETstar(e,*t_id,t,all_rels.etstars,mesh,gradient);
        }
    }
}

void Forman_Gradient_Topological_Relations::get_VTstar_ET(local_VTstar_ET &all_rels, Node_V &n, Mesh& mesh, Forman_Gradient &gradient)
{
    all_rels.init((n.get_v_end()-n.get_v_start()));
    ivect e;
    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);
        for(int v=0; v<t.vertices_num(); v++)
        {
            //if vertex is in the leaf node, check VT*
            if(n.indexes_vertex(t.TV(v)))
                Forman_Gradient_Topological_Relations::check_VTstar(t.TV(v),*t_id,t,all_rels.get_VTstars(),n.get_v_start(),mesh,gradient);
            //the first 4 edge are checked into the vertices loop
            t.TE(v,e);////corresponding edge
            if(n.indexes_vertex(e[1])) //if e[1] is in the current node, extract ET.
            {
                leaf_ET::iterator it = all_rels.find_ET(e);
                if(it != all_rels.end_ETs())
                {
                    ET &inside = it->second;
                    inside.second = *t_id;

                    /// force the ordering <max,min>
                    if(inside.first < inside.second)
                    {
                        itype tmp = inside.first;
                        inside.first = inside.second;
                        inside.second = tmp;
                    }
                }
                else
                {
                    ET new_entry = make_pair(*t_id,-1);
                    all_rels.add_ET(e,new_entry);
                }
            }
        }
    }
}


void Forman_Gradient_Topological_Relations::get_VTstar_ET(local_VTstar_ET &all_rels, Node_V &n, Mesh& mesh, Forman_Gradient &gradient, mig_cache& cache)
{
    int key = n.get_v_start() + n.get_v_end();
    auto ef_cache = cache.find_ef_lists(key);
    all_rels.init((n.get_v_end() - n.get_v_start()));
    bool has_ef_in_cache = false, has_vt_in_cache = false;
    if(ef_cache != cache.end_et_cache())
    {
        has_ef_in_cache = true;
        all_rels.set_ETs(ef_cache->second);
    }
    auto vt_cache = cache.find_vtstar_lists(key);
    if(vt_cache != cache.end_vtstar_cache()){
        has_vt_in_cache = true;
        all_rels.set_VTstars(vt_cache->second);
    }
    if(has_ef_in_cache && has_vt_in_cache)
        return;
    ivect e;
    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);
        for(int v=0; v < t.vertices_num(); v++)
        {
            //if vertex is in the leaf node, check VT*
            if(!has_vt_in_cache && n.indexes_vertex(t.TV(v)))
                Forman_Gradient_Topological_Relations::check_VTstar(t.TV(v),*t_id,t,all_rels.get_VTstars(),n.get_v_start(),mesh,gradient);
            //the first 4 edge are checked into the vertices loop
            t.TE(v,e);////corresponding edge
            if(!has_ef_in_cache && n.indexes_vertex(e[1])) //if e[1] is in the current node, extract ET.
            {
                leaf_ET::iterator it = all_rels.find_ET(e);
                if(it != all_rels.end_ETs())
                {
                    ET &inside = it->second;
                    inside.second = *t_id;
                    /// force the ordering <max,min>
                    if(inside.first < inside.second)
                    {
                        itype tmp = inside.first;
                        inside.first = inside.second;
                        inside.second = tmp;
                    }
                }
                else
                {
                    ET new_entry = make_pair(*t_id,-1);
                    all_rels.add_ET(e,new_entry);
                }
            }
        }
    }
    if(!has_ef_in_cache)
        cache.add_ef_lists(key,all_rels.get_ETs());
    if(!has_vt_in_cache)
        cache.add_vtstar_lists(key,all_rels.get_VTstars());
}

void Forman_Gradient_Topological_Relations::get_VTstar_ET(local_VTstar_ET &all_rels, Node_T &n, itype v_start, itype v_end, Mesh& mesh, Forman_Gradient &gradient)
{
    all_rels.init((v_end-v_start));
    ivect e;
    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);

        for(int v=0; v<t.vertices_num(); v++)
        {
            if(n.indexes_vertex(v_start,v_end,t.TV(v)))
                Forman_Gradient_Topological_Relations::check_VTstar(t.TV(v),*t_id,t,all_rels.get_VTstars(),v_start,mesh,gradient);
            //the first 4 edge are checked into the vertices loop
            t.TE(v,e);
            if(n.indexes_vertex(v_start,v_end,e[1]))
            {
                leaf_ET::iterator it = all_rels.find_ET(e);
                if(it != all_rels.end_ETs())
                {
                    ET &inside = it->second;
                    inside.second = *t_id;

                    /// force the ordering <max,min>
                    if(inside.first < inside.second)
                    {
                        itype tmp = inside.first;
                        inside.first = inside.second;
                        inside.second = tmp;
                    }
                }
                else
                {
                    ET new_entry = make_pair(*t_id,-1);
                    all_rels.add_ET(e,new_entry);
                }
            }
        }
    }
}

void Forman_Gradient_Topological_Relations::get_VTstar_VV(asc2rels &all_rels, Node_V &n, Mesh& mesh, Forman_Gradient &gradient)
{
    all_rels.init((n.get_v_end()-n.get_v_start()));

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);

        for(int v=0; v<t.vertices_num(); v++)
        {
            itype v_id = t.TV(v);
            if(n.indexes_vertex(t.TV(v)))
            {
                //init delle VT
                itype v_pos = v_id - n.get_v_start();
                for(int j=1; j<t.vertices_num(); j++)
                    all_rels.vvs[v_pos].insert(abs(t.TV((v+j)%t.vertices_num())));

                Forman_Gradient_Topological_Relations::check_VTstar(v_id,*t_id,t,all_rels.vtstars,n.get_v_start(),mesh,gradient);
            }
        }
    }
}

bool Forman_Gradient_Topological_Relations::check_VTstar(itype v_id, itype t_id, Triangle& t, leaf_VTstar &vtstars, itype v_start, Mesh &mesh, Forman_Gradient &gradient)
{
    itype v_pos = v_id - v_start;

    if(vtstars[v_pos] == -1)
    {
        //first one is simply set
        vtstars[v_pos] = t_id;
        return false;
    }

    TriGradient t_grad = gradient.convert_compressed_to_expand(t_id);
    //if in the current top simplex the vertex is unpaired.. nothing to do
    if(t_grad.is_vertex_unpaired(t.vertex_index(v_id)))
        return true;

    Triangle &first_t = mesh.get_triangle(vtstars[v_pos]);
    TriGradient t_grad_first = gradient.convert_compressed_to_expand(vtstars[v_pos]);

    //if we arrive here the current top simplex has the vertex paired
    //I have only to check if the first one has the vertex unpaired
    if(t_grad_first.is_vertex_unpaired(first_t.vertex_index(v_id)))
    {
        vtstars[v_pos] = t_id;
        return true;
    }

    return false;
}

bool Forman_Gradient_Topological_Relations::check_ETstar(ivect &e, itype t_id, Triangle& t, leaf_ETstar &etstars, Mesh &mesh, Forman_Gradient &gradient)
{
    leaf_ETstar::iterator it = etstars.find(e);
    if(it == etstars.end())
    {
        etstars.insert(make_pair(e,t_id));
    }
    else
    {
        TriGradient t_grad = gradient.convert_compressed_to_expand(t_id);
        if(t_grad.is_edge_unpaired(t.vertex_index(e[0]),t.vertex_index(e[1]))) //se non e' paired stop..
            return true;

        Triangle &current_etstar = mesh.get_triangle(it->second);
        TriGradient t_grad_first = gradient.convert_compressed_to_expand(it->second);

        if(t_grad_first.is_edge_unpaired(current_etstar.vertex_index(e[0]),current_etstar.vertex_index(e[1])))
        {
            etstars[e] = t_id;
            return true;
        }
    }
    return false;
}

itype Forman_Gradient_Topological_Relations::get_triangle_id(ivect max_tri, VT &vt, Mesh &mesh)
{
    sort(max_tri.begin(),max_tri.end());
    for(auto t_id : vt)
    {
        Triangle tri = mesh.get_triangle(t_id);
        ivect tri_vect;
        tri.convert_to_vec(tri_vect);
        sort(tri_vect.begin(),tri_vect.end());
        if(tri_vect == max_tri)
            return t_id;
    }
    // it should never reach this point if we give the correct VT relation
    return -1;
}

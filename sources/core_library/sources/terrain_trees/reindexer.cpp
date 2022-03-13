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

#include "reindexer.h"

void Reindexer::reindex_tree_and_mesh(PMRT_Tree& tree, bool save_v_indices, ivect &original_vertex_indices,
                                      bool save_t_indices, ivect &original_triangle_indices)
{
    coherent_indices.assign(tree.get_mesh().get_vertices_num(),-1);
    reindex_vertices(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_subdivision(),tree.get_mesh(),save_v_indices,original_vertex_indices);
    update_triangles_boundaries(tree.get_mesh());
    update_mesh_vertices(tree.get_mesh());
    reset();

    coherent_indices.assign(tree.get_mesh().get_triangles_num(),-1);
    extract_leaf_tri_association(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_subdivision(),tree.get_mesh(),tree.get_root());
    get_triangles_reordered_indexes(save_t_indices,original_triangle_indices);
    compress_tree_representation_top(tree.get_root());
    update_mesh_triangles(tree.get_mesh());

    reset();
    return;
}

void Reindexer::reindex_vertices(Node_T& n, Box &domain, int level, Spatial_Subdivision& decomposition, Mesh &mesh,
                                 bool save_v_indices, ivect &original_vertex_indices)
{
    if (n.is_leaf())
    {
        iset contained_vertices;
        //we recollect the vertices geometrically contained by the leaf block
        for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& t_id = itPair.first;
            Triangle& t = mesh.get_triangle(*t_id);
            for(int j=0;j<t.vertices_num();j++)
            {
                int real_v_index = t.TV(j);
                if(domain.contains(mesh.get_vertex(real_v_index),mesh.get_domain().get_max()))
                {                    
                    contained_vertices.insert(real_v_index);
                }
            }
        }
        //set the coherent position indices for those vertices
        for(iset_iter it=contained_vertices.begin(); it!=contained_vertices.end();it++)
        {
            coherent_indices[*it-1] = indices_counter;

            if(save_v_indices)
                original_vertex_indices[indices_counter - 1] = *it;

            indices_counter++;
        }
    }
    else
    {
        for (int i = 0; i < decomposition.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                Box son_dom = decomposition.compute_domain(domain,level,i);
                int son_level = level +1;
                reindex_vertices(*n.get_son(i),son_dom,son_level,decomposition,mesh,save_v_indices,original_vertex_indices);
            }
        }
    }
}

void Reindexer::reindex_vertices(Node_V &n, Spatial_Subdivision &decomposition, bool save_v_indices, ivect &original_vertex_indices)
{
    if (n.is_leaf())
    {
        //set the coherent position indices for the indexed vertices
        if(n.get_real_v_array_size()>0)
        {
            itype start = indices_counter;
            for(RunIteratorPair itPair = n.make_v_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
            {
                RunIterator const& v_id = itPair.first;
                coherent_indices[*v_id-1] = indices_counter;

                if(save_v_indices)
                    original_vertex_indices[indices_counter - 1] = *v_id;

                indices_counter++;
            }
            itype end = indices_counter;
            n.clear_v_array();
            n.set_v_range(start,end);
        }
    }
    else
    {
        itype start = indices_counter;
        for (int i = 0; i < decomposition.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                reindex_vertices(*n.get_son(i),decomposition,save_v_indices,original_vertex_indices);
            }
        }
        itype end = indices_counter;
        n.set_v_range(start,end);
    }
}

void Reindexer::update_triangles_boundaries(Mesh &mesh)
{
    for(itype i=1;i<=mesh.get_triangles_num();i++)
    {
        Triangle& t = mesh.get_triangle(i);
        for(int j=0;j<t.vertices_num();j++)
            t.setTV(j,coherent_indices[t.TV(j)-1]);
    }
}

void Reindexer::update_mesh_vertices(Mesh &mesh/*, bool save_original_indices, ivect &original_vertex_indices*/)
{
    // update global vertices array
    //     -1 identifies an already set vertex
    // ** the reordering is done inline the global vertex array ** ///
    for(itype i=1;i<=mesh.get_vertices_num();i++)
    {
        itype j = i -1; // -1 is needed to avoid array lookup error

        if(coherent_indices[j] == i || coherent_indices[j] < 0)
        {
            // mark the last entry visited...
            coherent_indices[j] = -1;
            continue;
        }

        while(coherent_indices[j] != i)
        {
            mesh.vertices_swap(i,coherent_indices[j]);
            itype j_prime = coherent_indices[j] -1; // -1 is needed to avoid array lookup error
            coherent_indices[j] = -1;
            j = j_prime;
        }
        // mark the last entry visited...
        coherent_indices[j] = -1;
    }
}

void Reindexer::extract_leaf_tri_association(Node_T& n, Box& dom, int level, Spatial_Subdivision& decomposition, Mesh &mesh, Node_T &root)
{
    if (n.is_leaf())
    {
        extract_leaf_tri_association_leaf(n,dom,decomposition,mesh,root);
    }
    else
    {
        for (int i = 0; i < decomposition.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                Box son_dom = decomposition.compute_domain(dom,level,i);
                int son_level = level +1;
                extract_leaf_tri_association(*n.get_son(i),son_dom,son_level,decomposition,mesh,root);
            }
        }
    }
}

void Reindexer::extract_leaf_tri_association_leaf(Node_T& n, Box& dom, Spatial_Subdivision& decomposition, Mesh &mesh, Node_T& root)
{
    itype v_start;
    itype v_end;

    n.get_v_range(v_start,v_end,dom,mesh); // we need to gather the vertices range..

    if(v_start == v_end) //no internal vertices..
        return;

    ivect internal_t_key; internal_t_key.push_back(v_start);

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle &t = mesh.get_triangle(*t_id);
        itype min_v = t.minindex();
        if(n.indexes_vertex(v_start,v_end,min_v))
        {
            ivect key = internal_t_key;
            pair<utype,utype> value = make_pair(this->leaf_tuples_array.size(),1);

            if(!n.completely_indexes_triangle_vertices_dom(t,dom,mesh))
            {
                ivect outer_v_ids;

                for(int i=0; i<t.vertices_num(); i++)
                {
                    if(!n.indexes_vertex(v_start,v_end,t.TV(i)))
                        outer_v_ids.push_back(t.TV(i));
                }

                get_tri_tuple_key(root,mesh.get_domain(),0,decomposition,mesh,outer_v_ids,key);
            }

            pair<map<ivect,pair<utype,utype> >::iterator,bool> ret = this->leaf_tuples_array.insert(make_pair(key,value));
            if(ret.second) // inserted
            {
                this->coherent_indices[*t_id-1] = ret.first->second.first;
            }
            else
            {
                this->coherent_indices[*t_id-1] = ret.first->second.first;
                ret.first->second.second++;
            }
        }
    }
}

void Reindexer::get_tri_tuple_key(Node_T &n, Box& dom, int level, Spatial_Subdivision& decomposition, Mesh &mesh, ivect &v_ids, ivect &key)
{
    if (n.is_leaf())
    {
        itype v_start;
        itype v_end;

        n.get_v_range(v_start,v_end,dom,mesh); // we need to gather the vertices range..

        // when I'm entering a leaf the v of a top are for sure indexed..
        key.push_back(v_start);
        // clean-up phase for v_ids array
        for(ivect_iter it=v_ids.begin(); it!=v_ids.end();)
        {
            if(n.indexes_vertex(v_start,v_end,*it))
            {
                v_ids.erase(it);
            }
            else
                ++it;
        }
    }
    else
    {
        for (int i = 0; i < decomposition.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                Box son_dom = decomposition.compute_domain(dom,level,i);
                int son_level = level +1;

                for(ivect_iter itv=v_ids.begin(); itv!=v_ids.end(); ++itv)
                {
                    if(son_dom.contains(mesh.get_vertex(*itv),mesh.get_domain().get_max()))
                    {
                        get_tri_tuple_key(*n.get_son(i),son_dom,son_level,decomposition,mesh,v_ids,key);
                        break;
                    }
                }
            }
        }
    }
}

void Reindexer::update_mesh_triangles(Mesh &mesh/*, bool save_original_indices, ivect &original_triangle_indices*/)
{
    for(itype i=1; i<=mesh.get_triangles_num(); i++)
    {
        itype j = i -1; // -1 is needed to avoid array lookup error

        if(coherent_indices[j] == i || coherent_indices[j] < 0)
        {
            // mark the last entry visited...
            coherent_indices[j] = -1;
            continue;
        }

        while(coherent_indices[j] != i && coherent_indices[j] > 0)
        {
            mesh.triangles_swap(i,coherent_indices[j]);
            int j_prime = coherent_indices[j] -1; // -1 is needed to avoid array lookup error
            coherent_indices[j] = -1;
            j = j_prime;
        }

        // mark the last entry visited...
        coherent_indices[j] = -1;
    }
}

void Reindexer::get_triangles_reordered_indexes(bool save_original_indices, ivect &original_triangle_indices)
{
    ivect I;
    I.assign(leaf_tuples_array.size(),0);
    utype counter = 1;

    //get the prefix sum of the leaf_tuple_array counts by iterating it
    //this will give us, for each grouped set of top cells, the initial index of this group in the reindexed array
    for(map<ivect,pair<utype,utype> >::iterator it=leaf_tuples_array.begin(); it!=leaf_tuples_array.end(); ++it)
    {
        I[it->second.first] = counter;
        counter += it->second.second;
    }
    leaf_tuples_array.clear();
    //we updated the values in coherent_indices setting the new index value for each triangle
    for(utype j=0; j<coherent_indices.size(); j++)
    {
        utype leaf_key = coherent_indices[j];
        coherent_indices[j] = I[leaf_key];

        if(save_original_indices)
            original_triangle_indices[I[leaf_key]-1] = j+1;

        I[leaf_key]++;
    }
}

/*
    This file is part of the Terrain Trees library.

    Author(s): Riccardo Fellegara (riccardo.fellegara@gmail.com)
               Federico Iuricich (federico.iuricich@gmail.com)

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

void Forman_Gradient_Features_Extractor::extract_descending_1cells(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, Spatial_Subdivision &division,
                                                                   Node_V &root, OpType operation, int cache_size)
{
    /// reset the statistical counters
    reset_extraction_critica_counters();

    leaves_1_desc_map dangling_paths = leaves_1_desc_map();
    vtstar_cache cache = vtstar_cache(cache_size);
    this->descending_1cells_extraction(n,mesh,gradient,division,root,operation,dangling_paths,cache);
    if(dangling_paths.size() != 0)
        cout<<"[ERROR] the procedure to not visit fully the dangling paths"<<endl;
    dangling_paths.clear();
    cache.reset();
    cout<<"one saddles visited: "<<found_1selle<<endl;
    found_1selle = 0;
}

void Forman_Gradient_Features_Extractor::descending_1cells_extraction(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, Spatial_Subdivision &division,  Node_V &root,
                                                                      OpType operation, leaves_1_desc_map &dangling_paths, vtstar_cache &cache)
{
    if (n.is_leaf())
    {
        this->descending_1cells_extraction_leaf(n,mesh,gradient,root,division,operation,dangling_paths,cache);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                this->descending_1cells_extraction(*n.get_son(i),mesh,gradient,division,root,operation,dangling_paths,cache);
            }
        }
    }
}

void Forman_Gradient_Features_Extractor::descending_1cells_extraction_leaf(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, Node_V &root, Spatial_Subdivision &division,
                                                                           OpType operation, leaves_1_desc_map &dangling_paths, vtstar_cache &cache)
{
    /// if there are no vertices in the leaf we have nothing to do..
    if(!n.indexes_vertices())
        return;

    Timer time = Timer();

    desc1rels all_rels;

    if(operation == TIME_VERBOSE)
        time.start();
    Forman_Gradient_Topological_Relations::get_VTstar_ETstar(all_rels,n,mesh,gradient);
    if(operation == TIME_VERBOSE)
    {
        time.stop();
        rels_gathering += time.get_elapsed_time();
        time.start();
    }

    get_new_descending_1cells(n,mesh,gradient,all_rels,root,division,operation,dangling_paths,cache);
    if(operation == TIME_VERBOSE)
    {
        time.stop();
        new_paths += time.get_elapsed_time();
        time.start();
    }

    get_dangling_descending_1cells(n,mesh,gradient,all_rels,root,division,operation,dangling_paths,cache);
    if(operation == TIME_VERBOSE)
    {
        time.stop();
        u_paths += time.get_elapsed_time();
        time.start();
    }
    //add the current vt to cache (in any leaves except the latest)
    if(mesh.get_vertices_num() > n.get_v_end())
        cache.insert(n.get_v_end()+n.get_v_start(),all_rels.vtstars);
    if(operation == TIME_VERBOSE)
    {
        time.stop();
        cache_t += time.get_elapsed_time();
    }

    //get stats
    if(operation == OUTPUT)
    {
        stats.set_cache_stats(cache);
        stats.set_1desc_dangling_paths_stats(dangling_paths);
    }
}

void Forman_Gradient_Features_Extractor::get_new_descending_1cells(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, desc1rels &all_rels, Node_V &root, Spatial_Subdivision &division,
                                                                   OpType operation, leaves_1_desc_map &dangling_paths, vtstar_cache &cache)
{
    for(leaf_ETstar::iterator it_e=all_rels.begin_etstars(); it_e!=all_rels.end_etstars(); ++it_e)
    {
        //each edge is visited once, precisely, inside the leaf block indexing the vertex, in its boundary, with the maximum position index
        if(gradient.is_edge_critical(it_e->first,it_e->second,mesh))
        {
            found_1selle++;
            Edge* crit_edge  = new Edge(it_e->first[0],it_e->first[1]);
            if(operation == OUTPUT)
            {
                edges_1celle[it_e->first].push_back(it_e->first[0]);
                // edges_1celle.insert(make_pair(it_e->first,it_e->first[0]));
            }
            get_one_descending_1cells(n, crit_edge, mesh, gradient, it_e->first[0], all_rels, root, division, operation, dangling_paths, cache);
        }
    }
}

void Forman_Gradient_Features_Extractor::get_dangling_descending_1cells(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, desc1rels &all_rels,
                                                                        Node_V &root, Spatial_Subdivision &division, OpType operation,
                                                                        leaves_1_desc_map &dangling_paths, vtstar_cache &cache)
{
    //search for path started in a previous leaf block
    int key = n.get_v_start()+n.get_v_end();

    leaves_1_desc_map::iterator it = dangling_paths.find(key);
    if(it != dangling_paths.end())
    {
        map_paths_pair &unfinished = it->second;

        for(map_paths_pair::iterator it_e = unfinished.begin(); it_e != unfinished.end(); ++it_e)
        {
            Edge* edge = new Edge(it_e->first.first,it_e->first.second);
            get_one_descending_1cells(n, edge, mesh, gradient, it_e->second, all_rels, root, division, operation, dangling_paths, cache);
        }
        unfinished.clear();
        dangling_paths.erase(it);
    }
}

void Forman_Gradient_Features_Extractor::get_one_descending_1cells(Node_V &n, Edge *edge, Mesh &mesh, Forman_Gradient &gradient, itype label,
                                                                   desc1rels &all_rels, Node_V &root, Spatial_Subdivision &division, OpType operation,
                                                                   leaves_1_desc_map &dangling_paths, vtstar_cache &cache)
{
    queue<Edge*> coda;
    coda.push(edge);

    while(!coda.empty())
    {
        Edge* e = coda.front();
        coda.pop();

        if(n.has_all_vertices_visited(*e))
        {
            for(int k=0;k<e->vertices_num();k++)
            {
                itype v = e->EV(k);
                itype vtstar = Forman_Gradient_Topological_Relations::get_VTstar(n,v,all_rels.vtstars,cache,root,division,mesh,gradient);

                if(!gradient.is_vertex_critical(v,vtstar,mesh)) //*it is paired
                {
                    Triangle &t = mesh.get_triangle(vtstar);
                    short v1i = gradient.convert_compressed_to_expand(vtstar).get_vertex_pair(t.vertex_index(v));
                    itype v1 = t.TV(v1i);
                    Edge* edg = new Edge(v,v1);
                    if(*edg != *e)
                    {
                        if(operation == OUTPUT)
                        {
                            ivect tmp = { edg->EV(0), edg->EV(1) };
                            edges_1celle[tmp].push_back(label);
                            // edges_1celle.insert(make_pair(tmp,label));
                        }
                        coda.push(edg);
                    }
                    else
                        delete edg;
                }
            }
        }
        else
        {
            Dangling_Paths_Handler::add_dangling_path(root,division,dangling_paths,e,label);
        }

        delete e;

        if(operation == OUTPUT && stats.max_queue_size < coda.size())
            stats.max_queue_size = coda.size();
    }
}

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

void Forman_Gradient_Features_Extractor::extract_ascending_2cells(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, Spatial_Subdivision &division,
                                                                  Node_V &root, OpType operation, int cache_size)
{
    /// reset the statistical counters
    reset_extraction_critica_counters();

    leaves_2_asc_map dangling_paths = leaves_2_asc_map();
    asc3_cache cache = asc3_cache(cache_size);
    this->ascending_2cells_extraction(n,mesh,gradient,division,root,operation,dangling_paths,cache);
    if(dangling_paths.size() != 0)
        cout<<"[ERROR] the procedure to not visit fully the dangling paths"<<endl;
    dangling_paths.clear();
    cache.reset();
    cout<<"min visited: "<<found_min<<endl;
    found_min = 0;
}

void Forman_Gradient_Features_Extractor::ascending_2cells_extraction(Node_V &n, Mesh &mesh,  Forman_Gradient &gradient,Spatial_Subdivision &division, Node_V &root,
                                                                     OpType operation, leaves_2_asc_map &dangling_paths, asc3_cache &cache)
{
    if (n.is_leaf())
    {
        this->ascending_2cells_extraction_leaf(n,mesh,gradient,root,division,operation,dangling_paths,cache);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                this->ascending_2cells_extraction(*n.get_son(i),mesh,gradient,division,root,operation,dangling_paths,cache);
            }
        }
    }
}

void Forman_Gradient_Features_Extractor::ascending_2cells_extraction_leaf(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, Node_V &root, Spatial_Subdivision &division, OpType operation,
                                                                          leaves_2_asc_map &dangling_paths, asc3_cache &cache)
{
    /// if there are no vertices in the leaf we have nothing to do..
    if(!n.indexes_vertices())
        return;

    Timer time = Timer();
    if(operation == TIME_VERBOSE)
        time.start();
    asc2rels all_rels;
    Forman_Gradient_Topological_Relations::get_VTstar_VV(all_rels,n,mesh,gradient);

    if(operation == TIME_VERBOSE)
    {
        time.stop();
        rels_gathering += time.get_elapsed_time();
        time.start();
    }

    get_new_ascending_2cells(n,mesh,gradient,all_rels,root,division,operation,dangling_paths,cache);
    if(operation == TIME_VERBOSE)
    {
        time.stop();
        new_paths += time.get_elapsed_time();
        time.start();
    }
    get_dangling_ascending_2cells(n,mesh,gradient,all_rels,root,division,operation,dangling_paths,cache);
    if(operation == TIME_VERBOSE)
    {
        time.stop();
        u_paths += time.get_elapsed_time();
        time.start();
    }
    //add the current vt to cache (in any leaves except the latest)
    if(mesh.get_vertices_num() > n.get_v_end())
        cache.addItem(n.get_v_end()+n.get_v_start(),all_rels);
    if(operation == TIME_VERBOSE)
    {
        time.stop();
        cache_t += time.get_elapsed_time();
    }

    //get stats
    if(operation == OUTPUT)
    {
        stats.set_cache_stats(cache);
        stats.set_2asc_dangling_paths_stats(dangling_paths);
    }
}

void Forman_Gradient_Features_Extractor::get_new_ascending_2cells(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, asc2rels &all_rels, Node_V &root, Spatial_Subdivision &division, OpType operation,
                                                                  leaves_2_asc_map &dangling_paths, asc3_cache &cache)
{
    for(itype v_id=n.get_v_start(); v_id<n.get_v_end(); v_id++)
    {
        itype v_pos = v_id - n.get_v_start();
        itype vtstar = all_rels.get_vtstar(v_pos);
        if(vtstar == -1){
            // This will lead to error in is_vertex_critical function if not handled. 
            // In other cases, the list of triangles or vv relations are traversed so the vertex must be connected. 
            cout << "Vertex " << v_id <<" is not in any triangle of the mesh." <<endl;
            continue; 
        }
        if(gradient.is_vertex_critical(v_id, vtstar, mesh))
        {
            found_min++;
            if(operation == OUTPUT)
            {
                manifold_2celle_asc.at(v_id-1) = v_id-1;
            }
            get_one_ascending_2cells(n,v_id,mesh,gradient,v_id,all_rels,root,division,operation, dangling_paths, cache);
        }
    }
}

void Forman_Gradient_Features_Extractor::get_dangling_ascending_2cells(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, asc2rels &all_rels, Node_V &root, Spatial_Subdivision &division,
                                                                       OpType operation, leaves_2_asc_map &dangling_paths, asc3_cache &cache)
{
    //quindi cerco path iniziati in altre foglie che continuano in quella corrente
    int key = n.get_v_start()+n.get_v_end();

    leaves_2_asc_map::iterator it = dangling_paths.find(key);
    if(it != dangling_paths.end())
    {
        map_paths_pair &unfinished = it->second;

        for(map_paths_pair::iterator it_e = unfinished.begin(); it_e != unfinished.end(); ++it_e)
        {
            itype vtstar = Forman_Gradient_Topological_Relations::get_VTstar(n,it_e->first.first,all_rels.vtstars,
                                                                               cache.get_vtstar_cache(),root,division,mesh,gradient);
            if(vtstar  == -1){
                cout << "Error getting vtstar when processing dangling cells." << endl;
            }
            if(!gradient.is_vertex_critical(it_e->first.first,vtstar,mesh)) //*it is paired
            {
                Triangle &t = mesh.get_triangle(vtstar);
                TriGradient grad = gradient.convert_compressed_to_expand(vtstar/*vv_inside.get_vtstar()*/);
                short v1i = grad.get_vertex_pair(t.vertex_index(it_e->first.first));
                itype v1 = t.TV(v1i);
                if(v1 == it_e->first.second)
                {
                    if(operation == OUTPUT)
                    {
                        manifold_2celle_asc[it_e->first.first-1] = it_e->second-1;
                    }
                    get_one_ascending_2cells(n, it_e->first.first, mesh, gradient, it_e->second, all_rels, root, division, operation, dangling_paths, cache);
                }
            }
        }
        unfinished.clear();

        dangling_paths.erase(it);
    }
}

void Forman_Gradient_Features_Extractor::get_one_ascending_2cells(Node_V &n, itype crit_v, Mesh &mesh, Forman_Gradient &gradient, itype label, asc2rels &all_rels,
                                                                  Node_V &root, Spatial_Subdivision &division, OpType operation, leaves_2_asc_map &dangling_paths, asc3_cache &cache)
{
    iqueue coda;

    coda.push(crit_v);

    while(!coda.empty())
    {
        itype current_v = coda.front();
        coda.pop();

        iset &vv = Forman_Gradient_Topological_Relations::get_VV(n,current_v,all_rels.vvs,cache,root,division,mesh);

        for(iset_iter it=vv.begin(); it!=vv.end(); ++it)
        {
            itype vv_vert = *it;
            if(n.visited_vertex(vv_vert))
            {
                itype vtstar = Forman_Gradient_Topological_Relations::get_VTstar(n,vv_vert,all_rels.vtstars,cache.get_vtstar_cache(),root,division,mesh,gradient);
                if(vtstar == -1){
                    cout << "Error when checking vtstar of vv_vert" << endl;
                }
                if(!gradient.is_vertex_critical(vv_vert,vtstar,mesh)) //*it is paired
                {
                    Triangle &t = mesh.get_triangle(vtstar);
                    short v1i = gradient.convert_compressed_to_expand(vtstar).get_vertex_pair(t.vertex_index(vv_vert));
                    if(v1i <= -1)
                    {
                        continue;
                    }
                    itype point_to = t.TV(v1i);
                    if(point_to == current_v)
                    {
                        if(operation == OUTPUT)
                        {
                            manifold_2celle_asc[vv_vert-1] = label-1;
                        }
                        coda.push(vv_vert);
                    }
                }
            }
            else
            {
                Dangling_Paths_Handler::add_dangling_path(root, division, dangling_paths, vv_vert, current_v, label);
            }
        }

        if(operation == OUTPUT && stats.max_queue_size < coda.size())
            stats.max_queue_size = coda.size();
    }
}

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

void Forman_Gradient_Features_Extractor::extract_ascending_1cells(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, Spatial_Subdivision &division,
                                                                  Node_V &root, OpType operation, int cache_size)
{
    /// reset the statistical counters
    reset_extraction_critica_counters();

    leaves_asc_map dangling_paths = leaves_asc_map();
    et_cache cache = et_cache(cache_size);
    this->ascending_1cells_extraction(n,mesh,gradient,division,root,operation,dangling_paths,cache);

    if(dangling_paths.size() != 0)
        cout<<"[ERROR] the procedure to not visit fully the dangling paths"<<endl;

    cache.reset();
    dangling_paths.clear();
    cout<<"one saddle visited: "<<found_1selle<<endl;
    found_1selle = 0;
}

void Forman_Gradient_Features_Extractor::ascending_1cells_extraction(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, Spatial_Subdivision &division, Node_V &root,
                                                                     OpType operation, leaves_asc_map &dangling_paths, et_cache &cache)
{
    if (n.is_leaf())
    {
        this->ascending_1cells_extraction_leaf(n,mesh,gradient,root,division,operation,dangling_paths,cache);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                this->ascending_1cells_extraction(*n.get_son(i),mesh,gradient,division,root,operation,dangling_paths,cache);
            }
        }
    }
}

void Forman_Gradient_Features_Extractor::ascending_1cells_extraction_leaf(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, Node_V &root, Spatial_Subdivision &division,
                                                                          OpType operation, leaves_asc_map &dangling_paths, et_cache &cache)
{
    /// if there are no vertices in the leaf we have nothing to do..
    if(!n.indexes_vertices())
        return;
    Timer time;

    leaf_ET local_ef;

    if(operation == TIME_VERBOSE)
        time.start();
    n.get_ET(local_ef,mesh);
    if(operation == TIME_VERBOSE)
    {
        time.stop();
        rels_gathering += time.get_elapsed_time();
        time.start();
    }
    get_new_ascending_1cells(n,mesh,gradient,local_ef,root,division,operation,dangling_paths,cache);
    if(operation == TIME_VERBOSE)
    {
        time.stop();
        new_paths += time.get_elapsed_time();
        time.start();
    }
    get_dangling_ascending_1cells(n,mesh,gradient,local_ef,root,division,operation,dangling_paths,cache);
    if(operation == TIME_VERBOSE)
    {
        time.stop();
        u_paths += time.get_elapsed_time();
        time.start();
    }
    //add the current vt to cache (in any leaves except the last)
    if(mesh.get_vertices_num() > n.get_v_end())
        cache.insert(n.get_v_end()+n.get_v_start(),local_ef);
    if(operation == TIME_VERBOSE)
    {
        time.stop();
        cache_t += time.get_elapsed_time();
    }

    //get stats
    if(operation == OUTPUT)
    {
        stats.set_cache_stats(cache);
        stats.set_1asc_dangling_paths_stats(dangling_paths);
    }
}

void Forman_Gradient_Features_Extractor::get_new_ascending_1cells(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, leaf_ET &local_ef, Node_V &root, Spatial_Subdivision &division,
                                                                  OpType operation, leaves_asc_map &dangling_paths, et_cache &cache)
{
    for(leaf_ET::iterator it_f=local_ef.begin(); it_f!=local_ef.end(); ++it_f)
    {
        if(!gradient.is_edge_critical(it_f->first,it_f->second,mesh))
            continue;

        /// if it_f is a critical edge visit the ascending 1 cells
        found_1selle++;

        Triangle& t_first = mesh.get_triangle(it_f->second.first);
        short f_pos_first = t_first.edge_index(it_f->first);
        short f_pos_second = -1;

        if(it_f->second.second != -1)
        {
            Triangle& t_second = mesh.get_triangle(it_f->second.second);
            f_pos_second = t_second.edge_index(it_f->first);
        }

        //salvo in strutture globali solo se i dati estratti mi servono per la visualizzazione
        if(operation == OUTPUT)
        {
            ivect f_vec;
            t_first.convert_to_vec(f_vec);
            triangles_2celle[f_vec].push_back(it_f->second.first);
            // triangles_2celle.insert(make_pair(f_vec, it_f->second.first));
        }

        pair<itype,short> coppia = make_pair(it_f->second.first,f_pos_first);
        get_one_ascending_1cells(n, coppia, mesh, gradient, it_f->second.first, local_ef, root, division,operation, dangling_paths, cache);
        if(it_f->second.second != -1)
        {
            pair<itype,short> coppia2 = make_pair(it_f->second.second,f_pos_second);
            get_one_ascending_1cells(n, coppia2, mesh, gradient, it_f->second.first, local_ef, root, division,operation, dangling_paths, cache);
        }
    }
}

void Forman_Gradient_Features_Extractor::get_dangling_ascending_1cells(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, leaf_ET &local_ef, Node_V &root, Spatial_Subdivision &division,
                                                                       OpType operation, leaves_asc_map &dangling_paths, et_cache &cache)
{
    //quindi cerco path iniziati in altre foglie che continuano in quella corrente
    itype key = n.get_v_start()+n.get_v_end();
    leaves_asc_map::iterator it = dangling_paths.find(key);
    if(it != dangling_paths.end())
    {
        triple_set &pairs = it->second;
        for(triple_set::iterator it_p = pairs.begin(); it_p != pairs.end();++it_p)
        {
            if(n.has_all_vertices_visited(mesh.get_triangle(it_p->get_first())))
                get_one_ascending_1cells(n, it_p->get_pair(), mesh, gradient, it_p->get_label(), local_ef, root, division,operation, dangling_paths, cache);            
        }
        pairs.clear();
        dangling_paths.erase(it);
    }
}

void Forman_Gradient_Features_Extractor::get_one_ascending_1cells(Node_V &n, const pair<itype, short> &tf, Mesh &mesh, Forman_Gradient &gradient,
                                                                  itype label, leaf_ET &local_ef, Node_V &root, Spatial_Subdivision &division,
                                                                  OpType operation, leaves_asc_map &dangling_paths, et_cache &cache)
{
    queue<pair<itype, short> > coda = queue<pair<itype, short> >();
    coda.push(tf);
    pair<itype, short> current_pair;

    ivect e;

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
                    if(operation == OUTPUT)
                    {
                        ivect f_vec;
                        t.convert_to_vec(f_vec);
                        triangles_2celle[f_vec].push_back(label);
                        // triangles_2celle.insert(make_pair(f_vec, label));
                    }

                    t.TE(f1,e);
                    pair<itype,itype> ft = Forman_Gradient_Topological_Relations::get_ET(n,e,local_ef,cache,root,division,mesh);

                    itype next = (ft.first == current_pair.first) ? ft.second : ft.first;
                    if(next != -1)
                    {
                        f1 = mesh.get_triangle(next).edge_index(e);
                        coda.push(make_pair(next,f1));
                    }
                }
                else if (operation == OUTPUT && gradient.is_triangle_critical(current_pair.first))
                {
                    ivect f_vec;
                    t.convert_to_vec(f_vec);
                    triangles_2celle[f_vec].push_back(label);
                    // triangles_2celle.insert(make_pair(f_vec, label));
                }
            }
        }
        else
        {
            itype max_id = t.maxindex();
            Dangling_Paths_Handler::add_dangling_path(root,division,dangling_paths,max_id,current_pair,label);
        }

        if(operation == OUTPUT && stats.max_queue_size < coda.size())
            stats.max_queue_size = coda.size();
    }
}

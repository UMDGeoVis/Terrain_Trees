#include "node_v.h"

void Node_V::get_VT(leaf_VT &all_vt, Mesh &mesh)
{
    // here we have a reindexed index thus, if there are no vertices indexed the array size is zero
    if(!indexes_vertices())
        return;

    itype v_start = get_v_start();
    itype v_end = get_v_end();

    all_vt.assign(v_end-v_start,ivect());

    
    for(RunIteratorPair itPair = make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);
        for(int v=0; v<t.vertices_num(); v++)
        {
            itype real_v_index = t.TV(v);
            //a vertex must be inside the leaf and inside the box (this avoids to insert the same vertex in different leaves)
            if (indexes_vertex(real_v_index))
                all_vt[real_v_index-v_start].push_back(*t_id);
        }
    }
}

void Node_V::get_VT_and_border(leaf_VT &all_vt, boost::dynamic_bitset<> &is_v_border,Mesh &mesh)
{
    // here we have a reindexed index thus, if there are no vertices indexed the array size is zero
    if(!indexes_vertices())
        return;

    itype v_start = get_v_start();
    itype v_end = get_v_end();

    all_vt.assign(v_end-v_start,ivect());
    
    for(RunIteratorPair itPair = make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        
        Triangle& t = mesh.get_triangle(*t_id);
        
         //cout<<"DEBUG:"<<"t_id:"<<*t_id<<" V:"<<t.TV(0)<<", "<<t.TV(1)<<", "<<t.TV(2)<<endl;
        
        
        for(int v=0; v<t.vertices_num(); v++)
        {
            itype real_v_index = t.TV(v);
          
            //a vertex must be inside the leaf and inside the box (this avoids to insert the same vertex in different leaves)
            if (indexes_vertex(real_v_index)){
                all_vt[real_v_index-v_start].push_back(*t_id);

                if(!is_v_border[real_v_index-v_start])
        {
            for(int v1=1; v1<t.vertices_num(); v1++)
            {
                if(t.is_border_edge((v1+v)%t.vertices_num()))
                {
                   // cout<<"[DEBUG] triangle id:"<<*t_id<<endl;
                    is_v_border[real_v_index-v_start] = true;
                   // cout<<"Vertex "<<real_v_index<<" is on the boundary."<<endl;
                    break;
                }
            }
        }
        }
        }
    }
}




void Node_V::get_VV(leaf_VV &all_vv,  Mesh& mesh)
{
    if(!this->indexes_vertices())
        return;

    all_vv.assign((this->get_v_end()-this->get_v_start()),iset());

    for(RunIteratorPair itPair = make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);
        for(int v=0; v<t.vertices_num(); v++)
        {
            itype v_id = t.TV(v);
            if(this->indexes_vertex(v_id))
            {
                //init VV
                itype v_pos = v_id - this->get_v_start();
                for(int j=1; j<t.vertices_num(); j++)
                    all_vv[v_pos].insert(t.TV((v+j)%t.vertices_num()));
            }
        }
    }
}


void Node_V::get_VV_vector(leaf_VV_vec &all_vv,  Mesh& mesh)
{
    if(!this->indexes_vertices())
        return;
    all_vv.assign(this->get_v_end()-this->get_v_start(),ivect());

    for(RunIteratorPair itPair = make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);
       
        for(int v=0; v<t.vertices_num(); v++)
        {
                                 
            itype v_id = t.TV(v);
            if(this->indexes_vertex(v_id))
            {
            itype v_pos = v_id - this->get_v_start();
            for(int v2=1;v2<t.vertices_num();v2++){// check if one of its adjacent vertex is on the border only when it is on the border. 
                 if(t.is_border_edge((v+v2)%t.vertices_num())) //so if one vertex index is negative, the other two are on the boundary
                 {
                     all_vv[v_pos].push_back(t.TV((v+t.vertices_num()-v2)%t.vertices_num()));    //only work for triangle
                 }
                    all_vv[v_pos].push_back(t.TV((v+v2)%t.vertices_num()));
             }      
                //init VV
                
              
            }
        }
    }
}

void Node_V::get_VV_VT(leaf_VV &all_vv, leaf_VT &all_vt, Mesh &mesh)
{
    if(!this->indexes_vertices())
        return;

    itype v_start = get_v_start();
    itype v_end = get_v_end();
    all_vt.assign(v_end-v_start,ivect());
    all_vv.assign(v_end-v_start,iset());

    for(RunIteratorPair itPair = make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);
        for(int v=0; v<t.vertices_num(); v++)
        {
            itype v_id = t.TV(v);
            if(this->indexes_vertex(v_id))
            {
                //init VV and VT
                itype v_pos = v_id - v_start;
                all_vt[v_pos].push_back(*t_id);
                for(int j=1; j<t.vertices_num(); j++)
                {
                    all_vv[v_pos].insert(t.TV((v+j)%t.vertices_num()));                    
                }

            }
        }
    }
}

void Node_V::get_ET(leaf_ET &ets,  Mesh &mesh)
{
    if(!this->indexes_vertices())
        return;

    ivect e;

    for(RunIteratorPair itPair = make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle& t = mesh.get_triangle(*t_id);
        for(int v=0; v<t.vertices_num(); v++)
        {
            //ET
            t.TE(v,e);

            // an edge is considered processed only if the current leaf block indexes the extreme with the higher position index
            // in this way each edge is processed once during each traversal of the tree
            if(indexes_vertex(e[1]))
            {
                leaf_ET::iterator it = ets.find(e);
                if(it != ets.end())
                {
                    pair<itype,itype> &inside = it->second;
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
                    /// we add to the local structure only the edges with an internal vertex
                    pair<itype,itype> new_entry = make_pair(*t_id,-1);
                    ets.insert(make_pair(e,new_entry));
                }
            }
        }
    }
}



void Node_V::update_vertex_indices(ivect &new_v_indices, itype &index_counter){

    ivect new_v_list;
    // itype start = index_counter;

    for(RunIteratorPair itPair = make_v_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& v_id = itPair.first;
        if(new_v_indices[*v_id-1] != -1){
        //    assert(index_counter == new_v_indices[*v_id-1]);
            new_v_list.push_back(new_v_indices[*v_id-1]);
            index_counter++;            
           // cout<<new_v_indices[*v_id-1]<<endl;
            }
    }
    this->clear_v_array();
    // itype end = index_counter;

    if(new_v_list.size()>0){    
        
    this->set_v_range(*new_v_list.begin(),*(new_v_list.end()-1)+1);    
    // cout<<start<<", "<<end<<endl;
    }

    // for(ivect_iter i=new_v_list.begin();i!=new_v_list.end();i++){
    //     this->add_vertex(*i);
    // }
}


void Node_V::compress_triangle_array( ivect &new_t_list){
    sort(new_t_list.begin(),new_t_list.end());
    int count=0;
    int start_t_id = new_t_list[0];

    if(new_t_list.size()==1)
    {
       this->add_triangle(start_t_id);
       return;
    }
    //otherwise visit new_t_list

    //now obtain the new encoding
    for(unsigned i=0; i<new_t_list.size(); i++)
    {
        if((i+1<new_t_list.size()) && (new_t_list[i]+1) == new_t_list[i+1]) //I have a consecutive range in the t_list
        {
            count++;
        }
        else //found a possible new range
        {
            if(count > 1) //more than two consecutive tetrahedra
            {
                //the range should be from start_t_id to start_t_id + count
                //example: 12 - 13 - 14 - 15
                //is encoded: -12 - 3
                //the first index is implicit
                this->add_triangle(-start_t_id);
                this->add_triangle(count);
            }
            else //less or equal to two
            {
                this->add_triangle(start_t_id);
                if(count==1)
                    this->add_triangle(start_t_id+count);
            }
            //re-init
            count = 0;
            start_t_id = new_t_list[i+1];
        }
    }
}


void Node_V::update_and_compress_triangles_arrays(ivect &new_t_positions, bool all_deleted)
{

        if(all_deleted)  
            this->clear_t_array();
        else
            this->update_and_compress_triangles_array(new_t_positions);
    
}

void Node_V::update_and_compress_triangles_array(ivect &new_indices)
{
    ivect t_list;
    //int d=1;
    for(RunIteratorPair itPair = this->make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        if(new_indices[*t_id-1] != -1) // the triangle still exists
        {
            t_list.push_back(new_indices[*t_id-1]);
        }
    }

    this->clear_t_array();
    if(t_list.size() > 0)
        this->compress_triangle_array(t_list);
}

void Node_V::compact_vertices_array(Mesh &mesh, ivect &surviving_vertices)
{
    ivect new_v_list;
    for(RunIteratorPair itPair = make_v_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& v_id = itPair.first;
        if(!mesh.is_vertex_removed(*v_id))
        {
            new_v_list.push_back(*v_id);
            surviving_vertices.push_back(*v_id);
        }
    }
    this->clear_v_array();
    this->set_v_array(new_v_list);
}
void Node_V::compact_triangle_array( Mesh &mesh)
{
    ivect t_list;

    for(RunIteratorPair itPair = this->make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        // NEW the top cell must be into the leaf block as well
        // (for reusing position indices)
        // NOTA -> this completely removes the run-lenght encoding
        if(!mesh.is_triangle_removed(*t_id) && this->indexes_triangle_vertices(mesh.get_triangle(*t_id)))
        {
            t_list.push_back(*t_id);
        }
    }

    this->clear_t_array();
    if(t_list.size() > 0)
        this->compress_triangle_array(t_list); // this sets again the SRE encoding
}
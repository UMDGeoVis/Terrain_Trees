#include "prt_tree.h"

PRT_Tree::PRT_Tree(int vertices_per_leaf, int sons_num) : Tree<Node_V>(sons_num)
{
    this->vertices_threshold = vertices_per_leaf;
    this->mesh = Mesh();
    this->root = Node_V();
}

void PRT_Tree::build_tree()
{
    for(itype i=1;i<=this->mesh.get_vertices_num(); i++)
    {
        this->add_vertex(this->root,this->mesh.get_domain(),0,i);
    }
    for(itype i=1;i<=this->mesh.get_triangles_num();i++)
    {
        this->add_triangle(this->root,this->mesh.get_domain(),0,i);
    }
}

void PRT_Tree::reinsert_triangles(){
    cout<<"total: "<<mesh.get_triangles_num()<<endl;
    for(itype i=1;i<=this->mesh.get_triangles_num();i++)
    {
       Node_V* ancestor;
       int anc_level = 0;
       Box block = find_lowest_common_ancestor(this->root, this->mesh.get_domain(), 0, i, ancestor, anc_level);

       this->add_triangle_new(*ancestor, block, anc_level, i);

    }
}

int PRT_Tree::visit_and_unify(Node_V &n, Mesh &mesh){
//// The unify function can be simpler. We have internal range and we can just delete sons if the internal range is less than the threshold. 
    int vertex_counter = 0;
    if(n.is_leaf())
    { 
        
        vertex_counter = n.get_real_v_array_size();//this->count_indexed_vertices(n,mesh);
        // cout<< "vertex counter (leaf):"<<vertex_counter<<endl;
    }
    else
    {
        for(Node_V::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL)
            {
                int local_num_v = this->visit_and_unify(**it,mesh);
                // cout<<**it<<endl;
                if(local_num_v == 0)
                {
                    // auto it = n.get_son(i);
                    // delete *it;
                    // *it = NULL;
                }
                else
                    vertex_counter += local_num_v;
            }
        }

        
        if(vertex_counter <= this->vertices_threshold)
        {
            // iset internal_triangles;
            // internal_triangles.assign(mesh.get_top_cells_types(),iset());
            pair<iset_iter,bool> coppia;

            n.clear_v_array(); // we have to reset the vertices lists as it contains the internal range of vertices
          
            for(Node_V::child_iterator it=n.begin(); it!=n.end(); ++it)
            {
                if(*it != NULL)
                {
                    // we reinsert the vertices
                    Node_V &son = **it;
                    for(RunIteratorPair itPair = son.make_v_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
                    {
                        RunIterator const& v_id = itPair.first;
                        //if(!mesh.is_vertex_removed(*v_id))
                        n.add_vertex(*v_id);
                    }



                    // for(RunIteratorPair itPair = son.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
                    // {
                    //     RunIterator const& t_id = itPair.first;
                    //     if(!mesh.is_triangle_removed(*t_id))
                    //     {
                    //         coppia = internal_triangles.insert(*t_id);
                    //         if(coppia.second == true)
                    //         {
                    //             n.add_triangle(*t_id);
                    //         }
                    //     }
                    // }
                    
                }
            }

            
        }
    }

    return vertex_counter;

}

// return the level of the ancestor node
Box PRT_Tree::find_lowest_common_ancestor(Node_V& n, Box& domain, int level, itype t, Node_V*  &ancestor, int& anc_level){
    if(n.is_leaf()){
        ancestor = &n;
        anc_level = level;
        return domain;
    }
    else{
        ancestor = &n;
        anc_level = level;
        Box smallest_block = domain;
        for (int i = 0; i < this->subdivision.son_number(); i++)
        {

            if(n.get_son(i)->completely_indexes_triangle_vertices(this->mesh.get_triangle(t)))
            {
                Box son_dom = this->subdivision.compute_domain(domain,level,i);
                int son_level = level +1;
                smallest_block = find_lowest_common_ancestor(*n.get_son(i), son_dom, son_level, t, ancestor,anc_level);
            }

        }
        return smallest_block;
    }

}

void PRT_Tree::add_vertex(Node_V& n, Box& domain, int level, itype v)
{
    if (n.is_leaf())
    {
        n.add_vertex(v);
        if(is_full(n))
            this->split(n,domain,level);
    }
    else
    {
        for (int i = 0; i < this->subdivision.son_number(); i++)
        {
            Box son_dom = this->subdivision.compute_domain(domain,level,i);
            int son_level = level +1;
            if (son_dom.contains(this->mesh.get_vertex(v),this->mesh.get_domain().get_max()))
            {
                this->add_vertex(*n.get_son(i),son_dom,son_level,v);
                break;
            }
        }
    }
}

void PRT_Tree::add_triangle(Node_V& n, Box& domain, int level, itype t)
{
    if (!Geometry_Wrapper::triangle_in_box_build(t,domain,this->mesh)) return;

    if(n.is_leaf())
        n.add_triangle(t);
    else
    {
        for(int i=0;i<this->subdivision.son_number();i++)
        {
            Box son_dom = this->subdivision.compute_domain(domain,level,i);
            int son_level = level +1;
            this->add_triangle(*n.get_son(i),son_dom,son_level,t);
        }
    }
}

void PRT_Tree::split(Node_V& n, Box& domain, int level)
{
    n.init_sons(this->subdivision.son_number());
    //initilize the son nodes
    for(int i=0;i<this->subdivision.son_number();i++)
    {
        Node_V* s = new Node_V();
        n.set_son(s,i);
        s->set_parent(&n);
    }

    //re-insert the vertices
    for(RunIterator runIt = n.v_array_begin_iterator(), runEnd = n.v_array_end_iterator(); runIt != runEnd; ++runIt)
        this->add_vertex(n,domain,level,*runIt);

    //re-insert the triangles
    for(RunIterator runIt = n.t_array_begin_iterator(), runEnd = n.t_array_end_iterator(); runIt != runEnd; ++runIt)
        this->add_triangle(n,domain,level,*runIt);

    //delete the arrays of the node
    n.clear_v_array();
    n.clear_t_array();
}

void PRT_Tree::build_tree(Soup &soup)
{
    this->mesh.set_domain(soup.get_domain());
    int counter = 1;

    for(int i=1; i<=soup.get_triangles_num(); i++)
    {
        Explicit_Triangle &etri = soup.get_triangle(i);
        ivect indexed_tri;
        for(int v=0; v<etri.vertices_num(); v++)
        {
            if (mesh.get_domain().contains(etri.get_vertex(v),mesh.get_domain().get_max()))
            {
                bool first_time = false;
                this->add_vertex_from_soup(this->root,mesh.get_domain(),0,etri.get_vertex(v),counter,indexed_tri,first_time);
                if(first_time)
                    counter++;
            }
        }
        Triangle t = Triangle(indexed_tri);
        mesh.add_triangle(t);
    }

    for(int i=1;i<=this->mesh.get_triangles_num();i++)
    {
        this->add_triangle(this->root,this->mesh.get_domain(),0,i);
    }
}

void PRT_Tree::add_vertex_from_soup(Node_V &n, Box &domain, int level, Vertex &v, itype vertex_index, ivect &indexed_tri, bool &first_time)
{
    if (n.is_leaf())
    {
        int ind = this->is_already_inserted(n,v);
        if(ind == -1)
        {
            ind = vertex_index;
            mesh.add_vertex(v);
            first_time = true;
            n.add_vertex(ind);
            /// in case of split we call the "usual" addVertex because we have all we need into the mesh variable
            if (is_full(n))
            {
                this->split(n,domain,level);
            }
        }
        /// we have to add the index of the vertex to the indexed representation of the index
        indexed_tri.push_back(ind);
    }
    else
    {
        for (int i = 0; i < this->subdivision.son_number(); i++)
        {
            Box son_dom = this->subdivision.compute_domain(domain,level,i);
            int son_level = level +1;
            if (son_dom.contains(v,this->mesh.get_domain().get_max()))
            {
                this->add_vertex_from_soup(*n.get_son(i),son_dom,son_level,v,vertex_index,indexed_tri,first_time);
            }
        }
    }
}

itype PRT_Tree::is_already_inserted(Node_V& n, Vertex &v)
{
    itype ind = -1;
    for(RunIteratorPair itPair = n.make_v_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        Vertex &vert = mesh.get_vertex(*itPair.first);
        if(v == vert)
        {
            ind = *itPair.first;
            break;
        }
    }
    return ind;
}

void PRT_Tree::build_tree_from_cloud(vertex_multifield &multifield)
{
    for(int i=1;i<=this->mesh.get_vertices_num(); i++)
    {
        this->add_vertex_from_cloud(this->root,this->mesh.get_domain(),0,mesh.get_vertex(i),i,multifield);
    }
}

void PRT_Tree::add_vertex_from_cloud(Node_V &n, Box &domain, int level, Vertex &v, itype vertex_index, vertex_multifield &multifield)
{
    if (n.is_leaf())
    {
        int ind = this->is_already_inserted(n,v);
        if(ind != -1) /// already indexed
        {
            multifield[ind].insert(v.get_z());
            ///flag the vertex as removed
            mesh.remove_vertex(vertex_index);
        }
        else
        {
            /// we have to add the index of the vertex to the indexed representation of the index
            n.add_vertex(vertex_index);
            /// and then initialize the corresponding entry in multifield
            dset mf;
            mf.insert(v.get_z());
            multifield.insert(make_pair(vertex_index,mf));

            /// in case of split we call the "usual" addVertex because we have all we need into the mesh variable
            if (is_full(n))
            {
                this->split(n,domain,level);
            }
        }
    }
    else
    {
        for (int i = 0; i < this->subdivision.son_number(); i++)
        {
            Box son_dom = this->subdivision.compute_domain(domain,level,i);
            int son_level = level +1;
            if (son_dom.contains(v,this->mesh.get_domain().get_max()))
            {
                this->add_vertex_from_cloud(*n.get_son(i),son_dom,son_level,v,vertex_index,multifield);
            }
        }
    }
}

void PRT_Tree::compact_vertices_lists(Node_V &n, Mesh &mesh, ivect &surviving_vertices){
     if (n.is_leaf())
    {
        n.compact_vertices_array(mesh,surviving_vertices);
    }
    else
    {
        for(Node_V::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL)
                this->compact_vertices_lists(**it,mesh,surviving_vertices);
        }
    }
}


void PRT_Tree::update_tree(Node_V &n, ivect &new_v_positions, ivect &new_t_positions, bool all_deleted, itype& index_counter)
{
    if (n.is_leaf())
    {
        n.update_vertex_indices(new_v_positions, index_counter);
        n.clear_t_array();
        // if(new_t_positions.size()!=0) // if not all the top simplices have been removed
        //     n.update_and_compress_triangles_arrays(new_t_positions,all_deleted);
    }
    else
    {
        itype start = index_counter;
        for(Node_V::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL)
                this->update_tree(**it,new_v_positions,new_t_positions,all_deleted,index_counter);
        }
        itype end = index_counter;
        n.clear_v_array();
        n.set_v_range(start,end);
    }
}

void PRT_Tree::update_vertex_index(Node_V &n, ivect &new_v_positions, itype& index_counter)
{
    if (n.is_leaf())
    {
        n.update_vertex_indices(new_v_positions, index_counter);
        n.clear_t_array();
        // if(new_t_positions.size()!=0) // if not all the top simplices have been removed
        //     n.update_and_compress_triangles_arrays(new_t_positions,all_deleted);
    }
    else
    {
        itype start = index_counter;
        for(Node_V::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL)
                this->update_vertex_index(**it,new_v_positions,index_counter);
        }
        itype end = index_counter;
        if((end-start)<= this->vertices_threshold){


            n.delete_sons();
        }
        n.clear_v_array();
        if(end-start>0)
            n.set_v_range(start,end);

    }
}



void PRT_Tree::get_leaf_indexing_vertex(Node_V &n, int v_id, Node_V *&res)
{
    if (n.is_leaf())
    {
        res = &n;
    }
    else
    {
        for(Node_V::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL && (*it)->indexes_vertex(v_id)) // NOTA
                this->get_leaf_indexing_vertex(**it,v_id,res);
        }
    }
}




void PRT_Tree::init_leaves_list(Node_V &n)
{
    if(!n.is_leaf())
    {
        for(Node_V::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL)
            {
                if((*it)->is_leaf())
                {
                    this->leaves.push_back(*it);
                }
                else
                {
                    this->init_leaves_list(**it);
                }
            }
        }
    }
    else{
        this->leaves.push_back(&n);
    }


}

void PRT_Tree::init_leaves_list_with_divisions(Node_V &n,  Box &n_dom, int level, Spatial_Subdivision &division, vector<Box>& leaf_divisions){
    if(!n.is_leaf())
    {
        int i = 0;
        for(Node_V::child_iterator it=n.begin(); it!=n.end(); ++it, ++i)
        {
            Box son_dom = division.compute_domain(n_dom,level,i);
            int son_level = level +1;
            if(*it != NULL)
            {
                if((*it)->is_leaf())
                {
                    this->leaves.push_back(*it);
                    leaf_divisions.push_back(son_dom);

                }
                else
                {
                    this->init_leaves_list_with_divisions(**it,son_dom,son_level,division,leaf_divisions);
                }
            }
        }
    }
    else{
        this->leaves.push_back(&n);
        leaf_divisions.push_back(n_dom);

    }


}

void PRT_Tree::add_triangle_new(Node_V& n, Box& domain, int level, itype t)
{

    if (!n.indexes_triangle_vertices(this->mesh.get_triangle(t))&&!Geometry_Wrapper::triangle_in_box_build(t,domain,this->mesh)) return;

    if(n.is_leaf())
        n.add_triangle(t);
    else
    {
        for(int i=0;i<this->subdivision.son_number();i++)
        {
            Box son_dom = this->subdivision.compute_domain(domain,level,i);
            int son_level = level +1;
            this->add_triangle_new(*n.get_son(i),son_dom,son_level,t);
        }
    }
}

#include "terrain_trees/mesh_updater.h"


void Mesh_Updater::clean_vertices_array(Mesh &mesh, ivect &new_v_positions, ivect &surviving_vertices)
{
    int v_counter = 1; // v_counter keeps the new vertex indexing while v is needed to set the new index
    new_v_positions.assign(mesh.get_vertices_num(),-1);

    // vector<Vertex> old_list = mesh.get_vertices_array();
    // mesh.reset_vertices();
    // mesh.reserve_vertices_space(surviving_vertices.size());

    for(ivect_iter it=surviving_vertices.begin(); it!=surviving_vertices.end(); ++it)
    {
        // mesh.add_vertex(old_list[*it-1]);
        new_v_positions[*it-1] = v_counter;
        v_counter++;
    }

    for(int i=1; i<=mesh.get_vertices_num(); i++)
    {
        int j = i -1; // -1 is needed to avoid array lookup error

        if(new_v_positions[j] == i || new_v_positions[j] < 0)
        {
            // mark the last entry visited...
            new_v_positions[j] = -1;
            continue;
        }

        while(new_v_positions[j] != i && new_v_positions[j] > 0)
        {
            mesh.vertices_swap(i,new_v_positions[j]);
            int j_prime = new_v_positions[j] -1; // -1 is needed to avoid array lookup error
            new_v_positions[j] = -1;
            j = j_prime;
        }

        // mark the last entry visited...
        new_v_positions[j] = -1;
    }
    cout<<"Finished reordering"<<endl;

    new_v_positions.assign(mesh.get_vertices_num(),-1);

    mesh.resize_vertex_array(v_counter-1);
    cout<<"Cleaned vertex array"<<endl;

    v_counter=1;
    for(ivect_iter it=surviving_vertices.begin(); it!=surviving_vertices.end(); ++it)
    {
        // mesh.add_vertex(old_list[*it-1]);
        new_v_positions[*it-1] = v_counter;
        v_counter++;
    }    

    // old_list.clear();
    //vector<Vertex>().swap(old_list);
}

 bool Mesh_Updater::update_and_clean_triangles_arrays(Mesh &mesh, ivect &new_v_positions, ivect &new_top_positions,
                                                                                                   itype simplification_counters)
{
    bool all_deleted = false;
 
        //ivect triangle_pos;
        all_deleted = update_and_clean_triangles_array(mesh,new_v_positions,new_top_positions,simplification_counters);
      //  new_top_positions=triangle_pos;
    
    return all_deleted;
}

bool Mesh_Updater::update_and_clean_triangles_array(Mesh &mesh, ivect &new_v_positions, ivect &new_top_positions, int counter)
{
    // step 0: if we have removed all the top d-cells, then, simply clear the corresponding array
    if(counter==mesh.get_triangles_num())
    {
        cout<<"Reset"<<endl;
        mesh.reset_triangles();
        return true;
    }
    else
    {
        cout<<"Clean"<<endl;
        clean_triangles_array(mesh,new_top_positions);
        update_triangles_boundary(mesh,new_v_positions);
    }
    return false;
}

 void Mesh_Updater::clean_triangles_array(Mesh &mesh,  ivect &new_triangle_positions)
{
    int t_counter = 1, t = 0; // t_counter keeps the new top indexing while t is needed to set the new index

    // (1) init the triangles position indices array
    new_triangle_positions.assign(mesh.get_triangles_num(),-1);
    // (2) get the new position indices for the triangles that are not deleted
    cout<<"(2) get the new position indices for the triangles that are not deleted"<<endl;
 
    for(auto it=mesh.get_t_array_begin(); it!=mesh.get_t_array_end(); ++it)  // TODO: Understand why index cannot work, should use iterator
    {
       
        if(!mesh.is_triangle_removed( *it))
        {
            new_triangle_positions[t] = t_counter;
            t_counter++;
        }
        t++;
    }
    // (3) reorder (inline) the triangles arrays
    cout<<"reorder (inline) the triangles arrays"<<endl;
    reorder_triangles_array(mesh,new_triangle_positions);
    // (4) delete the positions from t_counter to the end of the array
    cout<<"delete the positions from t_counter to the end of the array"<<endl;
    cout<<"New size:"<<t_counter-1<<endl;
    cout<<"Old size:"<<mesh.get_triangles_num()<<endl;
    //cout<<"[DEBUG]"<<endl;
    // for(auto it=mesh.get_t_array_begin();it!=mesh.get_t_array_end();it++){

    //     cout<<"Triangle id:"<<it-mesh.get_t_array_begin()<<endl;
    //     cout<<"TV:"<<it->TV(0)<<", "<<it->TV(1)<<", "<<it->TV(2)<<endl;
    // }
    mesh.resize_triangle_array(t_counter-1);
    cout<<"Cleaned triangle array"<<endl;
}



void Mesh_Updater::update_triangles_boundary(Mesh &mesh, ivect &new_v_positions)
{
    int t_id = 1;
    for(auto it=mesh.get_t_array_begin(); it!=mesh.get_t_array_end(); ++it) 
    {
        if(!mesh.is_triangle_removed(*it))
        {
            for(int v=0; v<3; v++)
            {
                if(new_v_positions[abs(it->TV(v))-1]==-1)
                {
                    cout<<"UPDATING A TOP-SIMPLEX WITH A DELETED VERTEX: "<<abs(it->TV(v))-1<<endl;
                   // cout<<t_id<<" "<<*it<<endl;
                    int a; cin>>a;
                }
                it->setTV_keep_border(v, new_v_positions[abs(it->TV(v))-1]);
            }
        }
        else
            cerr<<"[WARNING] update_top_d_cells_array does not removes deleted top d-cells"<<endl;

        t_id++;
    }

     cout<<"Updated triangle boundary"<<endl;
}

 void Mesh_Updater::reorder_vertices_array(Mesh &mesh, ivect &new_v_positions)
{
    // (2) update global vertices array
    //     -1 identifies an already set vertex
    // ** the reordering is done inline the global vertex array ** ///
    for(int i=1;i<=mesh.get_vertices_num();i++)
    {
        int j = i -1; // -1 is needed to avoid array lookup error

        if(new_v_positions[j] == i || new_v_positions[j] < 0)
        {
            // mark the last entry visited...
            new_v_positions[j] = -1;
            continue;
        }

        while(new_v_positions[j] != i)
        {
            mesh.vertices_swap(i,new_v_positions[j]);
            int j_prime = new_v_positions[j] -1; // -1 is needed to avoid array lookup error
            new_v_positions[j] = -1;
            j = j_prime;
        }
        // mark the last entry visited...
        new_v_positions[j] = -1;
    }
}



 void Mesh_Updater::reorder_triangles_array(Mesh &mesh,  ivect new_triangles_positions)
{
    for(int i=1; i<=mesh.get_triangles_num(); i++)
    {
        int j = i -1; // -1 is needed to avoid array lookup error

        if(new_triangles_positions[j] == i || new_triangles_positions[j] < 0)
        {
            // mark the last entry visited...
            new_triangles_positions[j] = -1;
            continue;
        }

        while(new_triangles_positions[j] != i && new_triangles_positions[j] > 0)
        {
            mesh.triangles_swap(i,new_triangles_positions[j]);
            int j_prime = new_triangles_positions[j] -1; // -1 is needed to avoid array lookup error
            new_triangles_positions[j] = -1;
            j = j_prime;
        }

        // mark the last entry visited...
        new_triangles_positions[j] = -1;
    }
    cout<<"Finished reordering"<<endl;
}
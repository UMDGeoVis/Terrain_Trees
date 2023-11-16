#include "sea_ice_processor.h"

simplices_multimap Sea_Ice_Processor::get_processed_triangles(Mesh& mesh){
    vector<pair<itype, itype>> arcs;
    // size_of_simplices.resize(extracted_triangles.size());

    vector<pair<Vertex, Vertex>> edges_of_ridge_lines; 
    // We will filter simplices with peak lower than mode later, so here we can set the original peak to be a value smaller than mode. 
    cout <<"extracted_triangles size: "<<extracted_triangles.size()<<endl;
    int roughness_fid = mesh.get_vertex(1).get_fields_num() - 2;
    cout << roughness_fid<<endl;
    if(this->roughness_filter){
        for(auto simplex = extracted_triangles.begin(); simplex != extracted_triangles.end();){
            if(!has_valid_roughness(*simplex, roughness_limit, roughness_fid, mesh)){
                simplex = extracted_triangles.erase(simplex);
                continue;
            }
            simplex++;
        }
        
    }
    cout << "Filtered triangle list size: "<< extracted_triangles.size() <<endl;
    length_of_simplices.resize(extracted_triangles.size(), 0);
    root_of_simplices.resize(extracted_triangles.size());
    edges_for_each_root.resize(extracted_triangles.size());
    peak_elevation_of_simplices.resize(extracted_triangles.size(), area_mode - 0.01); 
    for(int tid = 0; tid != extracted_triangles.size(); tid++){  
        // No longer use this filter.
        // if(is_lower_than_level(simplex->first, mode, mesh)) {
        //     simplex = extracted_triangles.erase(simplex);
        //     continue;
        // }
        ivect tri = extracted_triangles[tid];
        root_of_simplices[tid] = tid;
        for(int k = 0; k < 3; k++){
            coord_type v_elevation = mesh.get_vertex(tri[k]).get_z();
            peak_elevation_of_simplices[tid] = max(v_elevation, peak_elevation_of_simplices[tid]);
            itype v1 = min(tri[k], tri[(k + 1) % 3]); 
            itype v2 = max(tri[k], tri[(k + 1) % 3]);
            pair<itype, itype> edge = make_pair(v1, v2);
            auto it = root_of_edges.find(edge);
            if(it == root_of_edges.end()){
                root_of_edges[edge] = tid;
            }
            else{
                // each arc is one edge shared by two triangles in the same simplex.
                connect(it->second, tid, extracted_triangles[it->second], tri, mesh);
            }
        }
    }

    cout <<"arc size: "<<arcs.size()<<endl;
    cout <<"Finish connecting"<<endl;
    auto output = update_label();
    cout <<"Finish updating"<<endl;
    return output;
}

bool Sea_Ice_Processor::is_lower_than_level(const ivect& triangle, coord_type mode, Mesh& mesh){
    for(int i = 0; i < 3; i++){
        if(mesh.get_vertex(triangle[i]).get_z() > mode){
            return false;
        }
    }
    return true;
}

bool Sea_Ice_Processor::has_valid_roughness(const ivect& triangle, coord_type threshold, int roughness_fid, Mesh& mesh){
    for(int k = 0; k < 3; k++){
        if(mesh.get_vertex(triangle[k]).get_field(roughness_fid) < threshold){
            return false;
        }
    }
    return true;
}


void Sea_Ice_Processor::connect(itype simplex1, itype simplex2, const ivect& tri1, const ivect& tri2, Mesh& mesh){
    Vertex centroid1 = get_centroid(tri1, mesh);
    Vertex centroid2 = get_centroid(tri2, mesh);
    ridge_paths.push_back(make_pair(centroid1, centroid2));
    double distance = centroid1.distance(centroid2);
    itype root1 = find_root(simplex1);
    itype root2 = find_root(simplex2);
    length_of_simplices[root1] += distance;
    if(root1 == root2) {
        edges_for_each_root[root1].push_back(ridge_paths.size() - 1);
        return;
    }
    
    if(length_of_simplices[root1] < length_of_simplices[root2]){
        root_of_simplices[root1] = root2;
        // size_of_simplices[root2] += size_of_simplices[root1];
        length_of_simplices[root2] += length_of_simplices[root1];
        peak_elevation_of_simplices[root2] = max(peak_elevation_of_simplices[root2], peak_elevation_of_simplices[root1]);
        edges_for_each_root[root2].push_back(ridge_paths.size() - 1);
    }else{
        root_of_simplices[root2] = root1;
        // size_of_simplices[root1] += size_of_simplices[root2];
        length_of_simplices[root1] += length_of_simplices[root2];
        peak_elevation_of_simplices[root1] = max(peak_elevation_of_simplices[root2], peak_elevation_of_simplices[root1]);
        edges_for_each_root[root1].push_back(ridge_paths.size() - 1);
    }

}
itype Sea_Ice_Processor::find_root(itype simplex){
    if(root_of_simplices[simplex] != simplex){
        root_of_simplices[simplex] = find_root(root_of_simplices[simplex]);
    }
    return root_of_simplices[simplex];
}

simplices_multimap Sea_Ice_Processor::update_label(){
    simplices_multimap output;
    // for(auto simplex = extracted_triangles.begin(); simplex != extracted_triangles.end(); ){
    for(int tid = 0; tid < extracted_triangles.size(); tid++){
        itype root = find_root(tid);
        // cout << root << " size: "<<size_of_simplices[root]<<endl;
        // Filter ridges with peak lower than 0.6 meter higher than the level sea ice (to match the results from icesat 2 paper.)
        if( length_of_simplices[root] > length_limit && peak_elevation_of_simplices[root] > 0.6){
            output[extracted_triangles[tid]].push_back(root);
            for(auto eid:edges_for_each_root[tid]){
                filtered_ridge_paths.push_back(ridge_paths[eid]);
            }
        }
    }
    return output;
}
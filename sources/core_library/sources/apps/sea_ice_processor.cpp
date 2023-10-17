#include "sea_ice_processor.h"

simplices_multimap Sea_Ice_Processor::get_processed_triangles(coord_type mode, Mesh& mesh){
    vector<pair<itype, itype>> arcs;
    size_of_simplices.resize(extracted_triangles.size());
    root_of_simplices.resize(extracted_triangles.size());
    int tid = 0;
    cout <<"extracted_triangles size: "<<extracted_triangles.size()<<endl;
    for(auto simplex = extracted_triangles.begin(); simplex != extracted_triangles.end(); ){
        if(is_lower_than_level(simplex->first, mode, mesh)) {
            simplex = extracted_triangles.erase(simplex);
            continue;
        }
        size_of_simplices[tid] = compute_tri_area(simplex->first, mesh);
        root_of_simplices[tid] = tid;
        simplex->second[0] = tid;
        for(int k = 0; k < 3; k++){
            itype v1 = min(simplex->first[k], simplex->first[(k + 1) % 3]);
            itype v2 = max(simplex->first[k], simplex->first[(k + 1) % 3]);
            pair<itype, itype> edge = make_pair(v1, v2);
            auto it = root_of_edges.find(edge);
            if(it == root_of_edges.end()){
                root_of_edges[edge] = tid;
            }
            // else if(it->second != simplex->second[0]){
            else{
                arcs.push_back(make_pair(it->second, tid));
            }
        }
        tid++;
        // for(int i = 0; i < simplex->second.size(); i++){
        //     root_of_simplices[simplex->second[i]] = simplex->second[i];
        //     for(int j = i + 1; j < simplex->second.size(); j++){
        //         arcs.push_back(make_pair(simplex->second[i], simplex->second[j]));
        //     }
        // }
        simplex++;
    }
    cout <<"arc size: "<<arcs.size()<<endl;
    for(auto arc:arcs){
        connect(arc.first, arc.second);
    }
    cout <<"Finish connecting"<<endl;
    update_label();
    cout <<"Finish updating"<<endl;
    return extracted_triangles;
}

bool Sea_Ice_Processor::is_lower_than_level(const ivect& triangle, coord_type mode, Mesh& mesh){
    for(int i = 0; i < 3; i++){
        if(mesh.get_vertex(triangle[i]).get_z() > mode){
            return false;
        }
    }
    return true;
}

void Sea_Ice_Processor::connect(itype simplex1, itype simplex2){
    itype root1 = find_root(simplex1);
    itype root2 = find_root(simplex2);
    if(root1 == root2) return;
    // cout << "root1 "<<root1 << " size: "<<size_of_simplices[root1]<<endl;
    // cout << "root2 "<<root2 << " size: "<<size_of_simplices[root2]<<endl;
    if(size_of_simplices[root1] < size_of_simplices[root2]){
        root_of_simplices[root1] = root2;
        size_of_simplices[root2] += size_of_simplices[root1];
    }else{
        root_of_simplices[root2] = root1;
        size_of_simplices[root1] += size_of_simplices[root2];
    }

}
itype Sea_Ice_Processor::find_root(itype simplex){
    if(root_of_simplices[simplex] != simplex){
        root_of_simplices[simplex] = find_root(root_of_simplices[simplex]);
    }
    return root_of_simplices[simplex];
}

void Sea_Ice_Processor::update_label(){
    for(auto simplex = extracted_triangles.begin(); simplex != extracted_triangles.end(); ){
        itype root = find_root(simplex->second[0]);
        // cout << root << " size: "<<size_of_simplices[root]<<endl;
        if( size_of_simplices[root] <= 50){
            simplex = extracted_triangles.erase(simplex);
        }else{
            simplex->second[0] = root;
            simplex++;
        }
    }
}
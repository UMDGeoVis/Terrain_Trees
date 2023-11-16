/*
Test file for ridge extraction. 
*/

#include "utilities/utility_functions.h"

using namespace utility_functions;
struct Parameters{
    bool use_roughness = false;
    bool consider_Rayleigh_criterion = true;
};

template<class T> void load_tree(T& tree, cli_parameters &cli);
template<class T> void extract_ridges(T& tree, cli_parameters &cli, coord_type& mode_to_correct, double& length_limit, 
                const Parameters& parameters, coord_type roughness_radius, coord_type roughness_limit);
template<class T> void load_terrain(T& tree, cli_parameters &cli);
template<class T> void extract_features(T& tree, cli_parameters &cli, Forman_Gradient_Simplifier& forman_simplifier,
                Forman_Gradient& forman_gradient, string file_name, double& length_limit, const Parameters& parameters, coord_type roughness_limit);
template<class T> void compute_roughness(T& tree, cli_parameters &cli, coord_type& radius, const Point& origin);
Point standarize_input(Mesh& mesh, coord_type& mode_to_correct);
void reverse_mesh_coordinates(Mesh& mesh, const Point& origin);



int main(int argc, char** argv )
{
    cli_parameters cli;
    cli.mesh_path = argv[1];
    cli.division_type = QUAD;
    cli.crit_type = "pr";
    cli.v_per_leaf = atoi(argv[2]);
	cli.app_debug = OUTPUT;
	cli.persistence = atof(argv[3]);
    coord_type mode_to_correct = atof(argv[4]);
    Parameters sea_ice_parameters;
    // bool use_roughness = false;
    // bool consider_Rayleigh_criterion = true;
    coord_type length_limit = atof(argv[5]);
    if(argc >= 7 && strcmp(argv[6],"-r") == 0){
        sea_ice_parameters.use_roughness = true;
    }
    coord_type roughness_limit = atof(argv[7]);
    coord_type radius = 5;
    cerr<<"[OBJECTIVE] this unit-test generates a quadtrees based on the PR-T tree criterion. "
	    <<"then, it saves the index and the mesh in VTK format for visualization purposes and finally "
		<<"it computes the 1-ascending paths after the simplification. Note that the elevations in the outputs are corrected to relative height to the mode for the purpose of sea ice application."<<endl;
		
    cerr<<"[NOTA] all the generated files are saved in the 'data' folder"<<endl;
	
	setParameters(cli);
    PRT_Tree ptree = PRT_Tree(cli.v_per_leaf,cli.division_type);
    cerr<<"[GENERATION] PR-T tree"<<endl;
    
    extract_ridges(ptree,cli, mode_to_correct, length_limit, sea_ice_parameters, radius, roughness_limit);    

    return (EXIT_SUCCESS);
}


template<class T> void load_terrain(T& tree, cli_parameters &cli)
{
    if (!Reader::read_mesh(tree.get_mesh(), cli.mesh_path))
    {
        cout << "[ERROR] Loading mesh file. Execution Stopped." << endl;
        return;
    }
    cerr << "[MEMORY] peak for Indexing the terrain: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
}

template<class T> void extract_ridges(T& tree, cli_parameters &cli, coord_type& mode_to_correct, double& length_limit, 
    const Parameters& parameters, coord_type roughness_radius, coord_type roughness_limit)
{
    stringstream out;
    out << get_path_without_file_extension(cli.mesh_path);
    out << "_" << cli.persistence;
    out << "_l" << length_limit;
    load_terrain(tree, cli);
    Point origin  = standarize_input(tree.get_mesh(), mode_to_correct);

    //Calculate the Forman gradient vector
    Forman_Gradient forman_gradient = Forman_Gradient(tree.get_mesh().get_triangles_num());
    Forman_Gradient_Computation gradient_computation = Forman_Gradient_Computation();

    Timer time = Timer();

    time.start();
    gradient_computation.initial_filtering_IA(tree.get_mesh());
    time.stop();
    time.print_elapsed_time("[TIME] Initial filtering ");

    load_tree(tree,cli);
    if(parameters.use_roughness){
        compute_roughness(tree, cli, roughness_radius, origin);
    }
    reverse_mesh_coordinates(tree.get_mesh(), origin);
    // reverse_mesh_coordinates(tree.get_mesh(), origin);
    /// ---- FORMAN GRADIENT COMPUTATION --- ///
    gradient_computation.reset_filtering(tree.get_mesh(),cli.original_vertex_indices);

    cout<<"[NOTA] Computing the gradient field"<<endl;
    time.start();
    gradient_computation.compute_gradient_vector(forman_gradient,tree.get_root(),tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] computing gradient vector field ");


    /// ---- MORPHOLOGICAL SIMPLIFICATION --- ///        
    {
        Forman_Gradient_Simplifier forman_simplifier;
        forman_simplifier.set_filtration_vec(gradient_computation.get_filtration());

        // cout<<"--- Topological features BEFORE simplification ---"<<endl;
        // extract_features(tree, cli, forman_simplifier, forman_gradient,"before");
        ///
        /// firstly we extract the MIG
        ///
        cout<<"--- Morse Incidence Graph BEFORE simplification ---"<<endl;
        forman_simplifier.get_incidence_graph().init(); /// init again the base of the MIG
        forman_simplifier.extract_incidence_graph(tree.get_root(),tree.get_mesh(),forman_gradient,tree.get_subdivision(),OUTPUT,cli.cache_size); /// we force to keep the MIG structure
    
       // forman_simplifier.get_incidence_graph().print_stats(true);
        forman_simplifier.get_incidence_graph().print_stats(true);
        forman_simplifier.reset_stats();
        forman_simplifier.reset_output_structures(tree.get_mesh());
        cout<<"Simplify the forman gradient vector."<<endl;
        time.start();
        /// otherwise we simplify the gradient computing first a global MIG and then simplifying it and the gradient
        /// default behaviour with alltime!
        forman_simplifier.exec_topological_simplification(tree.get_root(),tree.get_mesh(),forman_gradient,tree.get_subdivision(),
                                                                cli.cache_size,cli.persistence);
        time.stop();
        time.print_elapsed_time("[TIME] simplify the gradient ");    
        forman_simplifier.print_simplification_stats();

        ///
        /// then we compute again and output the simplified mig
        ///
        cout<<"--- Morse Incidence Graph AFTER simplification ---"<<endl;
        forman_simplifier.reset_output_structures(tree.get_mesh());
        forman_simplifier.get_incidence_graph().init(); /// init again the MIG structures
        forman_simplifier.extract_incidence_graph(tree.get_root(),tree.get_mesh(),forman_gradient,tree.get_subdivision(),OUTPUT,cli.cache_size);
        forman_simplifier.get_incidence_graph().print_stats(true);
        forman_simplifier.reset_stats();

        Writer_Morse::write_critical_points_txt(out.str(),"simplified", cli.v_per_leaf, forman_simplifier.get_incidence_graph(),tree.get_mesh(),
                                          cli.original_vertex_indices,cli.original_vertex_fields,cli.rever_to_original);

        forman_simplifier.reset_output_structures(tree.get_mesh());
        forman_simplifier.reset_timer_variables();

        cout<<"--- Topological features AFTER simplification ---"<<endl;
        extract_features(tree, cli, forman_simplifier, forman_gradient, "after", length_limit, parameters, roughness_limit);
    }
}



template<class T> void load_tree(T& tree, cli_parameters &cli)
{
    Timer time;

    stringstream base_info;
    base_info << cli.v_per_leaf << " " << cli.t_per_leaf << " " << cli.crit_type << " ";
    stringstream base;
    base << get_path_without_file_extension(cli.mesh_path);

    stringstream tree_info;
    tree_info << base_info.str();

    if (cli.isTreeFile)
    {
        cout << "tree path: " << cli.tree_path << endl;
        time.start();
        if (!Reader::read_tree(tree, tree.get_root(), cli.tree_path))
        {
            cerr << "[ERROR] Loading .tree file. Regenerate the Terrain tree." << endl;
            return;
        }
        time.stop();
        tree_info << "[TIME] Loading tree from file ";
        time.print_elapsed_time(tree_info.str());
    }
    else
    {
        tree_info << "[TIME] Building ";

        stringstream out;
        out << base.str() << "_" << SpatialDecType2string(cli.division_type) << "_" << cli.crit_type;
        if (cli.crit_type == "pr")
            out << "_v_" << cli.v_per_leaf << "_.tree";
        else if (cli.crit_type == "pm")
            out << "_v_" << cli.v_per_leaf << "_t_" << cli.t_per_leaf << "_.tree";
        else if (cli.crit_type == "pmr")
            out << "_t_" << cli.t_per_leaf << "_.tree";
        cli.tree_path=out.str();

        if (!Reader::read_tree(tree, tree.get_root(), cli.tree_path))
        {
            cerr << "[ERROR] Loading .tree file." << endl;
            cerr << "[GENERATION] tree from triangle mesh" << endl;
            time.start();
            tree.build_tree();
            time.stop();
            time.print_elapsed_time(tree_info.str());
            Writer::write_tree(out.str(), tree.get_root(), tree.get_subdivision());
        }
        else
            cout << "[NOTICE] Found corresponding .tree file. Loaded tree from file successfully"<<endl;
        if(cli.app_debug == OUTPUT)
        {
            stringstream out2;
            out2 << base.str();
            if (cli.crit_type == "pr")
                out2 << "_" << SpatialDecType2string(cli.division_type) << "_" << cli.crit_type << "_v_" << cli.v_per_leaf << "_tree.vtk";
            else if (cli.crit_type == "pm")
                out2 << "_" << SpatialDecType2string(cli.division_type) << "_" << cli.crit_type << "_v_" << cli.v_per_leaf << "_t_" << cli.t_per_leaf << "_tree.vtk";
            else if (cli.crit_type == "pmr")
                out2 << "_" << SpatialDecType2string(cli.division_type) << "_" << cli.crit_type << "_t_" << cli.t_per_leaf << "_tree.vtk";

            Writer::write_tree_VTK(out2.str(),tree.get_root(),tree.get_subdivision(),tree.get_mesh());
        }
    }

    cerr << "[MEMORY] peak for encoding the Terrain tree: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;



    if(cli.reindex)
    {
        cerr<<"[REINDEXING] tree and triangle mesh"<<endl;

        cli.original_vertex_indices.assign(tree.get_mesh().get_vertices_num(),-1);
        if(cli.app_debug == OUTPUT)
            cli.original_triangle_indices.assign(tree.get_mesh().get_triangles_num(),-1);

        time.start();
        Reindexer reindexer = Reindexer();
        reindexer.reindex_tree_and_mesh(tree,(cli.original_vertex_indices.size() == tree.get_mesh().get_vertices_num()),cli.original_vertex_indices,
                                        (cli.original_triangle_indices.size() == tree.get_mesh().get_triangles_num()),cli.original_triangle_indices);
        time.stop();
        time.print_elapsed_time("[TIME] Index and Mesh Reindexing ");
    }


    cerr << "[MEMORY] peak for reindexing the Terrain tree: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;

    if (cli.is_index)
    {
        Statistics stats;
        stats.get_index_statistics(tree,cli.reindex);
    }
}

template<class T> void extract_features(T& tree, cli_parameters &cli,Forman_Gradient_Simplifier& forman_simplifier, Forman_Gradient& forman_gradient, string file_name,
    double& length_limit, const Parameters& parameters, coord_type roughness_limit)
{
    /// --- TOPOLOGY FEATURE EXTRACTION --- ///
    Timer time;
    stringstream out;
    out << get_path_without_file_extension(cli.mesh_path);
    out <<"_"<< file_name<< "_";
    out << cli.persistence;
    out << "_l"<<length_limit;
    if(parameters.use_roughness){
        out<< "_roughness_filtered_r" << roughness_limit;
    }
    if(parameters.consider_Rayleigh_criterion){
        out<<"_Rayleigh_enabled";
    }
    /// ---- ASCENDING 1 MANIFOLD EXTRACTION --- ///
    cout<<"[NOTA] Extract the ascending 1 manifolds."<<endl;
    time.start();
    forman_simplifier.extract_ascending_1cells(tree.get_root(),tree.get_mesh(),forman_gradient,tree.get_subdivision(),tree.get_root(),
                                                cli.app_debug,cli.cache_size);
    time.stop();


    forman_simplifier.print_stats();
    forman_simplifier.reset_stats();
    auto extracted_cells = forman_simplifier.get_extracted_cells(TRIANGLE);

    Sea_Ice_Processor processor(extracted_cells, length_limit, /*area_mode=*/0); // area_mode is zero if standardrze input is applied.
    processor.enable_peak_elevation_filter(true);
    processor.enable_roughness_filter(parameters.use_roughness, roughness_limit);
    auto updated_cells = processor.get_processed_triangles(tree.get_mesh());
    auto ridge_paths_edges = processor.get_ridge_paths_edges();

    Writer_Morse::write_asc1cells_line_VTK(out.str(),"asc1cells", cli.v_per_leaf, updated_cells 
                                        , tree.get_mesh(),ridge_paths_edges);
    Writer_Morse::write_asc1cells_VTK(out.str(),"asc1cells", cli.v_per_leaf, updated_cells 
                                        , tree.get_mesh(), cli.original_triangle_indices,
                                        cli.original_vertex_indices,cli.original_vertex_fields,cli.rever_to_original);
    Writer_Morse::write_asc1cells_PLY(out.str(),"asc1cells", cli.v_per_leaf,
                                        updated_cells , tree.get_mesh(), cli.original_triangle_indices,
                                        cli.original_vertex_indices,cli.original_vertex_fields,cli.rever_to_original);
    Writer_Morse::write_asc1cells_paths_WKT_CSV(out.str(),"asc1cells", cli.v_per_leaf, ridge_paths_edges, tree.get_mesh());
    // Writer_Morse::write_asc1cells_WKT_CSV(out.str(),"asc1cells", cli.v_per_leaf,
    //                                   updated_cells, tree.get_mesh(), cli.original_triangle_indices,
    //                                   cli.original_vertex_indices,cli.original_vertex_fields,cli.rever_to_original);
    Writer_Morse::write_asc1cells_vertices_CSV(out.str(),"asc1cells", cli.v_per_leaf, updated_cells, tree.get_mesh(), cli.original_triangle_indices,
                                        cli.original_vertex_indices,cli.original_vertex_fields,cli.rever_to_original);
    forman_simplifier.reset_output_structures(tree.get_mesh());



}

Point standarize_input(Mesh& mesh, coord_type& mode_to_correct)
{
    if(mesh.get_vertices_num() == 0) return Point(0, 0);
    Point origin = mesh.get_vertex(1);

    for(int i = 1; i <= mesh.get_vertices_num(); i++){
        Vertex old = mesh.get_vertex(i);
        mesh.get_vertex(i).set_c(0, old.get_x() - origin.get_x());
        mesh.get_vertex(i).set_c(1, old.get_y() - origin.get_y());
        mesh.get_vertex(i).set_c(2, old.get_z() - mode_to_correct);
    }
    Box old_domain = mesh.get_domain();
    cout<<"Original data domain:"<< old_domain.get_min() <<" -- "<<old_domain.get_max()<<endl;
    Point new_min = old_domain.get_min() - origin;
    Point new_max = old_domain.get_max() - origin;
    old_domain = Box(new_min, new_max);
    mesh.set_domain(old_domain);
    return origin;
}

void reverse_mesh_coordinates(Mesh& mesh, const Point& origin)
{

    for(int i = 1; i <= mesh.get_vertices_num(); i++){
        Vertex old = mesh.get_vertex(i);
        mesh.get_vertex(i).set_c(0, old.get_x() + origin.get_x());
        mesh.get_vertex(i).set_c(1, old.get_y() + origin.get_y());
        // mesh.get_vertex(i).set_c(2, old.get_z());
    }
    Box old_domain = mesh.get_domain();
    Point new_min = old_domain.get_min() + origin;
    Point new_max = old_domain.get_max() + origin;
    old_domain = Box(new_min, new_max);
    mesh.set_domain(old_domain);
}

template<class T> void compute_roughness(T& tree, cli_parameters &cli, coord_type& radius, const Point& origin)
{
    stringstream out;
    out << get_path_without_file_extension(cli.mesh_path);
    
    cout<<"[NOTA] compute roughness within a range"<<endl;
    Timer time;
    Roughness_circle roughness = Roughness_circle();
    time.start();
    roughness.compute(tree.get_root(),tree.get_mesh(), radius, tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] roughness computation: ");
     cerr << "[MEMORY] peak for computing Roughness: " <<
        to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
    roughness.print_roughness_stats(tree.get_mesh(),tree.get_mesh().get_vertex(1).get_fields_num() - 1);
}
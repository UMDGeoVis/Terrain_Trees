#include "utilities/utility_functions.h"

using namespace utility_functions;


template<class T> void load_terrain(T& tree, cli_parameters &cli);
template<class T> void load_tree(T& tree, cli_parameters &cli);
template<class T> void morse_features_extraction(T& tree, cli_parameters &cli);

int main(int argc, char** argv)
{
    cli_parameters cli;
	cli.mesh_path = "../data/eggs64z4.tri";
    cli.division_type = QUAD;
    cli.crit_type = "pr";
    cli.v_per_leaf = 10;
	cli.app_debug = OUTPUT;
		
    cerr<<"[OBJECTIVE] this unit-test generates a quadtrees based on the PR-T tree criterion. "
	    <<"then, it saves the index and the mesh in VTK format for visualization purposes and finally "
		<<"it computes first the Morse gradient vector and finally it extracts all the morphological features, outputting all the extracted in separate VTK files."<<endl;		
    cerr<<"[NOTA] the input mesh must be pre-filtered, i.e. in a vertices neighborhood the z-coordinate must be unique"<<endl;
    cerr<<"[NOTA2] all the generated files are saved in the 'data' folder"<<endl;
	
	setParameters(cli);
    PRT_Tree ptree = PRT_Tree(cli.v_per_leaf,cli.division_type);
    cerr<<"[GENERATION] PR-T tree"<<endl;
   // load_tree(ptree,cli);
    morse_features_extraction(ptree,cli);    

    return (EXIT_SUCCESS);
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
            cerr << "[ERROR] Loading .tree file. Execution Stopped." << endl;
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
            cerr << "[NOTICE] Cannot find existed .tree file." << endl;
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
            Writer::write_mesh_VTK(base.str(),tree.get_mesh());
        }
    }

    cerr << "[MEMORY] peak for encoding the Terrain tree: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;



    if(cli.reindex)
    {
        cerr<<"[REINDEXING] tree and triangle mesh"<<endl;


        //        if((cli.query_type == MORSE_ANALYSIS || cli.query_type == LOCAL_MORSE_SIMPLIFICATION || cli.query_type == GLOBAL_MORSE_SIMPLIFICATION)
        //                && cli.app_debug == OUTPUT)
        cli.original_vertex_indices.assign(tree.get_mesh().get_vertices_num(),-1);
        if(cli.query_type == MORSE_ANALYSIS && cli.app_debug == OUTPUT)
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

template<class T> void morse_features_extraction(T& tree, cli_parameters &cli)
{
    stringstream out;
    out << get_path_without_file_extension(cli.mesh_path);
    load_terrain(tree,cli);
    //CALCOLO IL FORMAN GRADIENT VECTOR
    Forman_Gradient forman_gradient = Forman_Gradient(tree.get_mesh().get_triangles_num());
    Forman_Gradient_Computation gradient_computation = Forman_Gradient_Computation();
    Forman_Gradient_Features_Extractor features_extractor;    

    Timer time = Timer();

    ///--- INITIAL FILTERING --- ///
    time.start();
    gradient_computation.initial_filtering_IA(tree.get_mesh());
    time.stop();
    time.print_elapsed_time("[TIME] Initial filtering ");

    load_tree(tree,cli);

    //WARNING AFTER THIS THE FILTRATION ARRAY IN forman_gradient AND THE POSITION INDICES OF THE VERTICES OF THE TIN WILL NOT BE ALIGNED..
    // YOU HAVE TO USE original_vertex_indices FOR FETCHING THE CORRECT FILTRATION VALUE OF A VERTEX
    gradient_computation.reset_filtering(tree.get_mesh(),cli.original_vertex_indices, true);


    /// ---- FORMAN GRADIENT COMPUTATION --- ///
    cout<<"[NOTA] Compute the gradient field"<<endl;
    time.start();
    gradient_computation.compute_gradient_vector(forman_gradient,tree.get_root(),tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] computing gradient vector field ");

    /// ---- EXTRACT MORPHOLOGICAL FEATURES ---- ///    
    {
        /// ---- DESCENDING 2 MANIFOLD EXTRACTION --- ///
        cout<<"[NOTA] Extract the descending 2 manifolds."<<endl;
        if(cli.app_debug == OUTPUT)
            features_extractor.init_segmentation_vector(tree.get_mesh());
        time.start();
        features_extractor.extract_descending_2cells(tree.get_root(),tree.get_mesh(),forman_gradient,tree.get_subdivision(),tree.get_root(),
                                                     cli.app_debug,cli.cache_size);
        time.stop();

        if(cli.app_debug == OUTPUT) //get statistics
        {
            features_extractor.print_stats();
            features_extractor.reset_stats();
            Writer_Morse::write_desc2cells_VTK(out.str(),"desc2cells", cli.v_per_leaf,
                                         features_extractor.get_segmentation_vector(),tree.get_mesh(),cli.original_triangle_indices,
                                         cli.original_vertex_indices,cli.original_vertex_fields,cli.rever_to_original);
            features_extractor.reset_output_structures(tree.get_mesh());
        }
        
        /// ---- DESCENDING 1 MANIFOLD EXTRACTION --- ///
        cout<<"[NOTA] Extract the descending 1 manifolds."<<endl;
        time.start();
        features_extractor.extract_descending_1cells(tree.get_root(),tree.get_mesh(),forman_gradient,tree.get_subdivision(),tree.get_root(),
                                                     cli.app_debug,cli.cache_size);
        time.stop();

        if(cli.app_debug == OUTPUT)
        {
            features_extractor.print_stats();
            features_extractor.reset_stats();
            Writer_Morse::write_desc1cells_VTK(out.str(),"desc1cells", cli.v_per_leaf,
                                         features_extractor.get_extracted_cells(EDGE), tree.get_mesh(), cli.original_vertex_indices,
                                         cli.original_vertex_fields,cli.rever_to_original);
            features_extractor.reset_output_structures(tree.get_mesh());
        }
        
        /// ---- ASCENDING 2 MANIFOLD EXTRACTION --- ///
        cout<<"[NOTA] Extract the ascending 2 manifolds."<<endl;
        if(cli.app_debug == OUTPUT)
            features_extractor.init_ascending_segmentation_vector(tree.get_mesh());
        time.start();
        features_extractor.extract_ascending_2cells(tree.get_root(),tree.get_mesh(),forman_gradient,tree.get_subdivision(),tree.get_root(),
                                                    cli.app_debug,cli.cache_size);
        time.stop();

        if(cli.app_debug == OUTPUT)
        {
            features_extractor.print_stats();
            features_extractor.reset_stats();
            Writer_Morse::write_asc2cells_VTK(out.str(),"asc2cells", cli.v_per_leaf,
                                        features_extractor.get_ascending_segmentation(), tree.get_mesh(), cli.original_vertex_indices,
                                        cli.original_vertex_fields,cli.rever_to_original);
            features_extractor.reset_output_structures(tree.get_mesh());
        }
        
        /// ---- ASCENDING 1 MANIFOLD EXTRACTION --- ///
        cout<<"[NOTA] Extract the ascending 1 manifolds."<<endl;
        time.start();
        features_extractor.extract_ascending_1cells(tree.get_root(),tree.get_mesh(),forman_gradient,tree.get_subdivision(),tree.get_root(),
                                                    cli.app_debug,cli.cache_size);
        time.stop();

        if(cli.app_debug == OUTPUT)
        {
            features_extractor.print_stats();
            features_extractor.reset_stats();
            Writer_Morse::write_asc1cells_VTK(out.str(),"asc1cells", cli.v_per_leaf,
                                         features_extractor.get_extracted_cells(TRIANGLE), tree.get_mesh(), cli.original_triangle_indices,
                                         cli.original_vertex_indices,cli.original_vertex_fields,cli.rever_to_original);
            features_extractor.reset_output_structures(tree.get_mesh());
        }
        
        /// ---- MORSE INCIDENCE GRAPH COMPUTATION --- ///
        cout<<"[NOTA] Extract the Morse Incidence Graph."<<endl;
        time.start();
        features_extractor.extract_incidence_graph(tree.get_root(),tree.get_mesh(),forman_gradient,tree.get_subdivision(),cli.app_debug,cli.cache_size);
        time.stop();

        if(cli.app_debug == OUTPUT)
        {
            features_extractor.print_stats();
            features_extractor.reset_stats();
            Writer_Morse::write_incidence_graph_VTK(out.str(),"mig", cli.v_per_leaf, features_extractor.get_incidence_graph(),tree.get_mesh(),
                                              cli.original_vertex_indices,cli.original_vertex_fields,cli.rever_to_original);
            features_extractor.reset_output_structures(tree.get_mesh());
        }       
    }
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
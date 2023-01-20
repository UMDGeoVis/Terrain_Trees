/*
Test file for topological simplification module. Currently only the global version is working.
*/

#include "utilities/utility_functions.h"

using namespace utility_functions;

template<class T> void load_tree(T& tree, cli_parameters &cli);
template<class T> void morse_simplification(T& tree, cli_parameters &cli);
template<class T> void load_terrain(T& tree, cli_parameters &cli);


int main(int argc, char** argv )
{
    cli_parameters cli;
    cli.mesh_path = argv[2];
    cli.division_type = QUAD;
    cli.crit_type = "pr";
    cli.v_per_leaf = 200;
	cli.app_debug = OUTPUT;
	cli.persistence = atof(argv[3]);
    if(strcmp(argv[1],"local")==0)
       cli.query_type = LOCAL_MORSE_SIMPLIFICATION;//LOCAL_ or GLOBAL_MORSE_SIMPLIFICATION
	else if(strcmp(argv[1],"global")==0)
        cli.query_type=GLOBAL_MORSE_SIMPLIFICATION;
    else
    {
        cout<<"Please enter the type of morse simplification"<<endl;
        return 1;
    }
    
    cerr<<"[OBJECTIVE] this unit-test generates a quadtrees based on the PR-T tree criterion. "
	    <<"then, it saves the index and the mesh in VTK format for visualization purposes and finally "
		<<"it computes the Morse gradient vector and simplifies it, outputting the initial Morse IG and the simplified one in VTK format."<<endl;
		
    cerr<<"[NOTA] all the generated files are saved in the 'data' folder"<<endl;
	
	setParameters(cli);
    PRT_Tree ptree = PRT_Tree(cli.v_per_leaf,cli.division_type);
    cerr<<"[GENERATION] PR-T tree"<<endl;

    morse_simplification(ptree,cli);    

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

template<class T> void morse_simplification(T& tree, cli_parameters &cli)
{
    stringstream out;
    out << get_path_without_file_extension(cli.mesh_path);

    load_terrain(tree,cli);


    //CALCOLO IL FORMAN GRADIENT VECTOR
    Forman_Gradient forman_gradient = Forman_Gradient(tree.get_mesh().get_triangles_num());
    Forman_Gradient_Computation gradient_computation = Forman_Gradient_Computation();

    Timer time = Timer();

    time.start();
    gradient_computation.initial_filtering_IA(tree.get_mesh());
    time.stop();
    time.print_elapsed_time("[TIME] Initial filtering ");

    load_tree(tree,cli);

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

        ///
        /// firstly we extract the MIG
        ///
        cout<<"--- Morse Incidence Graph BEFORE simplification ---"<<endl;
        forman_simplifier.get_incidence_graph().init(); /// init again the base of the MIG
        forman_simplifier.extract_incidence_graph(tree.get_root(),tree.get_mesh(),forman_gradient,tree.get_subdivision(),OUTPUT,cli.cache_size); /// we force to keep the MIG structure
        

        Writer_Morse::write_incidence_graph_VTK(out.str(),"mig", cli.v_per_leaf, forman_simplifier.get_incidence_graph(),tree.get_mesh(),
                                          cli.original_vertex_indices,cli.original_vertex_fields,cli.rever_to_original); /// and we save it

       // features_extractor.get_incidence_graph().print_stats(true);
        forman_simplifier.get_incidence_graph().print_stats(true);
        forman_simplifier.reset_stats();
        forman_simplifier.reset_output_structures(tree.get_mesh());
        ///
        /// then we execute effectively the topological simplification
        ///
        /// we have chosen a fully local simplification. i.e. the MIG is computed locally
        if (cli.query_type==LOCAL_MORSE_SIMPLIFICATION)
        {
            cout<<"[LOCALLY] Simplify the forman gradient vector."<<endl;
            time.start();
            forman_simplifier.exec_local_topological_simplification(tree.get_root(),tree.get_mesh(),forman_gradient,tree.get_subdivision(),
                                                                    cli.app_debug,cli.cache_size,cli.persistence);
            time.stop();
            time.print_elapsed_time("[TIME] simplify the gradient ");
        }
        else if(cli.query_type==GLOBAL_MORSE_SIMPLIFICATION){

             cout<<"[GLOBALLY] Simplify the forman gradient vector."<<endl;
                      time.start();
            /// otherwise we simplify the gradient computing first a global MIG and then simplifying it and the gradient
            /// default behaviour with alltime!
            forman_simplifier.exec_global_topological_simplification(tree.get_root(),tree.get_mesh(),forman_gradient,tree.get_subdivision(),
                                                                     cli.app_debug,cli.cache_size,cli.persistence);
            time.stop();
            time.print_elapsed_time("[TIME] simplify the gradient ");
        }
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

        if(cli.query_type==LOCAL_MORSE_SIMPLIFICATION)
        Writer_Morse::write_incidence_graph_VTK(out.str(),"simplified_mig_local", cli.v_per_leaf, forman_simplifier.get_incidence_graph(),tree.get_mesh(),
                                          cli.original_vertex_indices,cli.original_vertex_fields,cli.rever_to_original);
        else if (cli.query_type==GLOBAL_MORSE_SIMPLIFICATION)
         Writer_Morse::write_incidence_graph_VTK(out.str(),"simplified_mig_global", cli.v_per_leaf, forman_simplifier.get_incidence_graph(),tree.get_mesh(),
                                          cli.original_vertex_indices,cli.original_vertex_fields,cli.rever_to_original);

        forman_simplifier.reset_output_structures(tree.get_mesh());
        forman_simplifier.reset_timer_variables();
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
            cerr << "[ERROR] Loading .tree file." << endl;
            cerr << "[GENERATION] tree from triangle mesh" << endl;
            time.start();
            tree.build_tree();
            time.stop();
            time.print_elapsed_time(tree_info.str());
         //   Writer::write_tree(out.str(), tree.get_root(), tree.get_subdivision());
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
            // Writer::write_mesh_VTK(base.str(),tree.get_mesh());
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

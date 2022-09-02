/*
            This file is part of the Terrain Trees library.

            Author(s): Riccardo Fellegara (riccardo.fellegara@gmail.com)

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


#include "utilities/utility_functions.h"
#include "core_library/sources/queries/topological_queries.h" ///TEMPORARY COMMAND FOR JOURNAL PAPER experiments

using namespace utility_functions;

template<class T> void load_terrain(T& tree, cli_parameters &cli);
template<class T> void load_tree(T& tree, cli_parameters &cli);
template<class T> void load_tree_lite(T& tree, cli_parameters &cli);
template<class T> void morse_features_extraction_and_simplification(T& tree, cli_parameters &cli);
template<class T> void multi_morse_terrain_analysis(T& tree, cli_parameters &cli);
void multi_morse_terrain_analysis(PMRT_Tree &tree, cli_parameters &cli);

/// CUSTOM COMMAND --- TEMPORARY COMMAND FOR JOURNAL PAPER experiments
template<class T> void custom_execution(T& tree, cli_parameters &cli);
void custom_execution(PMRT_Tree& tree, cli_parameters &cli);

template<class T> void SF_test(T& tree, cli_parameters &cli);
void SF_test(PMRT_Tree& tree, cli_parameters &cli);
int main(int argc, char** argv)
{
    if(argc == 1)
    {
        print_help();
        return 0;
    }

    cli_parameters cli = cli_parameters();

    if (read_arguments(argc, argv,cli) == -1)
    {
        print_usage();
        return (EXIT_FAILURE);
    }

    if (cli.isTreeFile)
    {
        //nel caso di lettura da file recupero le info dal nome del file
        setParameters(cli);
    }

    //controllo che tutto sia inizializzato correttamente
    if (!cli.is_getInput)
    {
        if (!checkParameters(cli))
            return (EXIT_FAILURE);
    }

    if (cli.crit_type == "pr")
    {
        PRT_Tree tree = PRT_Tree(cli.v_per_leaf,cli.division_type);
        //        load_tree(tree,cli);
        if(cli.reindex)
        {
             morse_features_extraction_and_simplification(tree,cli);
            //  multi_morse_terrain_analysis(tree,cli);
            //  custom_execution(tree,cli);
            SF_test(tree,cli);
        }
    }
    else if (cli.crit_type == "pm")
    {
        PMT_Tree tree = PMT_Tree(cli.v_per_leaf, cli.t_per_leaf,cli.division_type);
        //        load_tree(tree,cli);
        if(cli.reindex)
        {
               morse_features_extraction_and_simplification(tree,cli);
            // multi_morse_terrain_analysis(tree,cli);
            //   custom_execution(tree,cli);
            SF_test(tree,cli);
        }
    }
    else if (cli.crit_type == "pmr")
    {
        PMRT_Tree tree = PMRT_Tree(cli.t_per_leaf,cli.division_type);
        //        load_tree(tree,cli);
        if(cli.reindex)
        {
            //   multi_morse_terrain_analysis(tree,cli);

            //    custom_execution(tree,cli);
            SF_test(tree,cli);
        }
    }

    return (EXIT_SUCCESS);
}

template<class T> void load_tree(T& tree, cli_parameters &cli)
{
    Timer time;
    Soup soup;

    //[1] read the input
    if(get_file_extension(cli.mesh_path) == "soup")
    {
        /// triangle mesh as soup
        if (!Reader::read_soup(soup, cli.mesh_path))
        {
            cerr << "[ERROR] Loading soup file. Execution Stopped." << endl;
            return;
        }
    }
    else if(cli.query_type == FILTER)
    {
        ///read a point cloud
        if(!Reader::read_vertices(tree.get_mesh(),cli.mesh_path))
        {
            cerr << "[ERROR] Loading point cloud file. Execution Stopped." << endl;
            return;
        }
    }
    /// tri or off format triangle mesh
    else if (!Reader::read_mesh(tree.get_mesh(), cli.mesh_path))
    {
        cout << "[ERROR] Loading mesh file. Execution Stopped." << endl;
        return;
    }

    cerr << "[MEMORY] peak for Indexing the terrain: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;

    stringstream base_info;
    base_info << cli.v_per_leaf << " " << cli.t_per_leaf << " " << cli.crit_type << " ";
    stringstream base;
    base << get_path_without_file_extension(cli.mesh_path);

    stringstream tree_info;
    tree_info << base_info.str();

    /// if we want to filter a point clouds we force the construction of the tree
    if (cli.isTreeFile && cli.query_type != FILTER)
    {
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
//        time.start();
        if(get_file_extension(cli.mesh_path) == "soup")
        {
            cerr<<"[GENERATION] tree from soup of triangles"<<endl;
            time.start();
            tree.build_tree(soup);
            soup.clear();
            time.stop();
            Writer::write_mesh(base.str(),"from_soup",tree.get_mesh(),false);
        }
        else if(cli.query_type == FILTER)
        {
            /// here we generate a PR tree on the points
            /// we filter them and finally we output two files containing
            /// the clouds in a SpatialHadoop can understand
            /// and a multifield version of the cloud for future usage
            cerr<<"[GENERATION] tree from points cloud"<<endl;
            vertex_multifield multifield;
            time.start();
            tree.build_tree_from_cloud(multifield);
            time.stop();
            time.print_elapsed_time(tree_info.str());
            cerr<<"[STATS] vertices: "<<tree.get_mesh().get_vertices_num()<<" vs. multifield_vertices: "<<multifield.size()<<endl;
            cerr<<"[OUTPUT] writing points cloud files"<<endl;
            Writer::write_filtered_points_cloud(base.str(),tree.get_mesh());
            Writer::write_multifield_points_cloud(base.str(),multifield,tree.get_mesh());
            return; /// after generated these files we have done..
        }
        else
        {
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
                cerr<<"[GENERATION] tree from triangle mesh"<<endl;
                time.start();
                tree.build_tree();
                time.stop();
                Writer::write_tree(out.str(), tree.get_root(), tree.get_subdivision());
            }
            else
                cout << "[NOTICE] Found corresponding .tree file. Loaded tree from file successfully"<<endl;
        }

        time.print_elapsed_time(tree_info.str());

        //        stringstream out;
        //        out << base.str() << "_" << SpatialDecType2string(cli.division_type) << "_" << cli.crit_type;
        //        if (cli.crit_type == "pr")
        //            out << "_v_" << cli.v_per_leaf << "_.tree";
        //        else if (cli.crit_type == "pm")
        //            out << "_v_" << cli.v_per_leaf << "_t_" << cli.t_per_leaf << "_.tree";
        //        else if (cli.crit_type == "pmr")
        //            out << "_t_" << cli.t_per_leaf << "_.tree";
        //        Writer::write_tree(out.str(), tree.get_root(), tree.get_subdivision());

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

    if(cli.query_type == MORSE_ANALYSIS)
    {
        cli.rever_to_original = Reader::read_noised_field_value(tree.get_mesh(),base.str(),cli.original_vertex_fields);
    }

    if(cli.reindex)
    {
        cerr<<"[REINDEXING] tree and triangle mesh"<<endl;


        if((cli.query_type == MORSE_ANALYSIS || cli.query_type == LOCAL_MORSE_SIMPLIFICATION || cli.query_type == GLOBAL_MORSE_SIMPLIFICATION)
                && cli.app_debug == OUTPUT)
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

template<class T> void morse_features_extraction_and_simplification(T& tree, cli_parameters &cli)
{
    /// TO-DO: define a better granularity
    if(cli.query_type != MORSE_ANALYSIS && cli.query_type != LOCAL_MORSE_SIMPLIFICATION && cli.query_type != GLOBAL_MORSE_SIMPLIFICATION)
        return;

    stringstream out;
    out << get_path_without_file_extension(cli.mesh_path);
    Timer time = Timer();

    
    load_terrain(tree,cli);
    //CALCOLO IL FORMAN GRADIENT VECTOR
    Forman_Gradient forman_gradient = Forman_Gradient(tree.get_mesh().get_triangles_num());
    //    Forman_Gradient_Computation gradient_computation = Forman_Gradient_Computation(cli.original_vertex_indices);
    Forman_Gradient_Features_Extractor features_extractor;

    Forman_Gradient_Computation gradient_computation = Forman_Gradient_Computation();

    time.start();
    gradient_computation.initial_filtering_IA(tree.get_mesh());
    time.stop();
    time.print_elapsed_time("[TIME] Initial filtering ");

    load_tree_lite(tree,cli);
    gradient_computation.reset_filtering(tree.get_mesh(),cli.original_vertex_indices, true);

    /// ---- FORMAN GRADIENT COMPUTATION --- ///
    cout<<"[NOTA] Compute the gradient field"<<endl;
    time.start();
    gradient_computation.compute_gradient_vector(forman_gradient,tree.get_root(),tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] computing gradient vector field ");
    cerr << "[MEMORY] peak for computing the Forman Gradient: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;

//    ivect edge = {11376, 11377};
//    ET et = {22741, 22728};
//    cout<<"IS THE TARGET EDGE CRITICAL? "<<forman_gradient.is_edge_critical(edge,et,tree.get_mesh())<<endl;

    /// ---- EXTRACT MORPHOLOGICAL FEATURES ---- ///
    if(cli.query_type == MORSE_ANALYSIS)
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
        else //get timings
        {
            time.print_elapsed_time("[TIME] extract descending 2-cells ");
            if(cli.app_debug == TIME_VERBOSE)
            {
                features_extractor.print_feature_extraction_time();
                features_extractor.reset_timer_variables();
            }
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
        else //get timings
        {
            time.print_elapsed_time("[TIME] extract descending 1-cells ");
            if(cli.app_debug == TIME_VERBOSE)
            {
                features_extractor.print_feature_extraction_time();
                features_extractor.reset_timer_variables();
            }
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
        else //get timings
        {
            time.print_elapsed_time("[TIME] extract ascending 2-cells ");
            if(cli.app_debug == TIME_VERBOSE)
            {
                features_extractor.print_feature_extraction_time();
                features_extractor.reset_timer_variables();
            }
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
        else //get timings
        {
            time.print_elapsed_time("[TIME] extract ascending 1-cells ");
            if(cli.app_debug == TIME_VERBOSE)
            {
                features_extractor.print_feature_extraction_time();
                features_extractor.reset_timer_variables();
            }
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
        else
        {
            time.print_elapsed_time("[TIME] extract incidence graph ");
            if(cli.app_debug == TIME_VERBOSE)
            {
                features_extractor.print_feature_extraction_time();
                features_extractor.reset_timer_variables();
            }
        }
    }

    /// ---- MORPHOLOGICAL SIMPLIFICATION --- ///
    if(cli.query_type == LOCAL_MORSE_SIMPLIFICATION || cli.query_type == GLOBAL_MORSE_SIMPLIFICATION)
    {
        Forman_Gradient_Simplifier forman_simplifier;

        ///
        /// firstly we extract the MIG
        ///
        cout<<"--- Morse Incidence Graph BEFORE simplification ---"<<endl;
        forman_simplifier.get_incidence_graph().init(); /// init again the base of the MIG
        forman_simplifier.extract_incidence_graph(tree.get_root(),tree.get_mesh(),forman_gradient,tree.get_subdivision(),OUTPUT,cli.cache_size);
        /// we force to keep the MIG structure

        if(cli.app_debug == TIME_VERBOSE)
            Writer_Morse::write_incidence_graph_VTK(out.str(),"mig", cli.v_per_leaf, forman_simplifier.get_incidence_graph(),tree.get_mesh(),
                                                    cli.original_vertex_indices,cli.original_vertex_fields,cli.rever_to_original); /// and we save it

        forman_simplifier.reset_stats();
        forman_simplifier.reset_output_structures(tree.get_mesh());
        ///
        /// then we execute effectively the topological simplification
        ///
        if(cli.query_type == LOCAL_MORSE_SIMPLIFICATION) /// we have chosen a fully local simplification. i.e. the MIG is computed locally
        {
            cli.persistence = 0.8;
            cout<<"[LOCALLY] Simplify the forman gradient vector."<<endl;
            time.start();
            forman_simplifier.exec_local_topological_simplification(tree.get_root(),tree.get_mesh(),forman_gradient,tree.get_subdivision(),
                                                                    cli.app_debug,cli.cache_size,cli.persistence);
            time.stop();
            time.print_elapsed_time("[TIME] simplify the gradient ");
        }
        else if(cli.query_type == GLOBAL_MORSE_SIMPLIFICATION)
        {
            cout<<"[GLOBALLY] Simplify the forman gradient vector."<<endl;
            /// otherwise we simplify the gradient computing first a global MIG and then simplifying it and the gradient
            /// default behaviour with alltime!
            forman_simplifier.exec_global_topological_simplification(tree.get_root(),tree.get_mesh(),forman_gradient,tree.get_subdivision(),
                                                                     cli.app_debug,cli.cache_size,cli.persistence);
        }

        forman_simplifier.print_simplification_stats();

        if(cli.app_debug == TIME_VERBOSE)
        {
            forman_simplifier.print_feature_extraction_time();
            forman_simplifier.reset_timer_variables();
            cerr<<"--- --- --- --- ---"<<endl;
            forman_simplifier.print_stats();
            forman_simplifier.reset_stats();
        }

        ///
        /// then we compute again and output the simplified mig
        ///
        cout<<"--- Morse Incidence Graph AFTER simplification ---"<<endl;
        forman_simplifier.reset_output_structures(tree.get_mesh());
        forman_simplifier.get_incidence_graph().init(); /// init again the MIG structures
        forman_simplifier.extract_incidence_graph(tree.get_root(),tree.get_mesh(),forman_gradient,tree.get_subdivision(),OUTPUT,cli.cache_size);
        forman_simplifier.reset_stats();

        if(cli.app_debug == TIME_VERBOSE)
            Writer_Morse::write_incidence_graph_VTK(out.str(),"simplified_mig", cli.v_per_leaf, forman_simplifier.get_incidence_graph(),tree.get_mesh(),
                                                    cli.original_vertex_indices,cli.original_vertex_fields,cli.rever_to_original);

        forman_simplifier.reset_output_structures(tree.get_mesh());
        forman_simplifier.reset_timer_variables();
    }
}

template<class T> void multi_morse_terrain_analysis (T& tree, cli_parameters &cli)
{
    if(cli.query_type != MULTIVARIATE_MORSE_ANALYSIS)
        return;

    Timer time = Timer();

    if(tree.get_mesh().get_vertex(1).get_fields_num() == 1)
        /// only the z-coordinate.. we compute at least the curvature to have a second value for each vertex
    {
        C_Curvature ccurvature = C_Curvature(MEAN,true);
        time.start();
        ccurvature.compute(tree.get_root(),tree.get_mesh(),tree.get_subdivision());
        time.stop();
        time.print_elapsed_time("[TIME] mean-cCurvature computation: ");
    }

    Forman_Gradient forman_gradient = Forman_Gradient(tree.get_mesh().get_triangles_num());
    Forman_Gradient_Computation gradient_computation = Forman_Gradient_Computation();

    cout<<"[NOTA] Compute the gradient field"<<endl;
    time.start();
    gradient_computation.compute_gradient_vector(forman_gradient,tree.get_root(),tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] computing gradient vector field ");
    cerr << "[MEMORY] peak for computing the MultiVariate Forman Gradient: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;

    /// this extraction must be executed at the end..
    /// as it invalidates the critical d.s.
    cout<<"[NOTA] Extracting Critical Clusters"<<endl;
    Forman_Gradient_Features_Extractor features_extractor;
    features_extractor.extract_critical_clusters(tree.get_root(),tree.get_mesh(),forman_gradient,gradient_computation.get_critical_simplices(),
                                                 tree.get_subdivision(), cli.app_debug,cli.cache_size,get_path_without_file_extension(cli.mesh_path));
    cerr << "[MEMORY] peak for extracting the Critical Clusters: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;

    if(cli.app_debug == OUTPUT)
    {
        cout<<"[OUTPUT] saving the terrain with the extra fields"<<endl;
        Writer::write_mesh(get_path_without_file_extension(cli.mesh_path),"extra_fields",tree.get_mesh(),true);
    }
}

void multi_morse_terrain_analysis(PMRT_Tree &tree, cli_parameters &cli)
{
    if(cli.query_type != MULTIVARIATE_MORSE_ANALYSIS)
        return;

    Timer time = Timer();

    if(tree.get_mesh().get_vertex(1).get_fields_num() == 1)
        /// only the z-coordinate.. we compute at least the curvature to have a second value for each vertex
    {
        C_Curvature ccurvature = C_Curvature(MEAN,true);
        time.start();
        ccurvature.compute(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_subdivision());
        time.stop();
        time.print_elapsed_time("[TIME] mean-cCurvature computation: ");
    }

    Forman_Gradient forman_gradient = Forman_Gradient(tree.get_mesh().get_triangles_num());
    Forman_Gradient_Computation gradient_computation = Forman_Gradient_Computation();

    cout<<"[NOTA] Compute the gradient field"<<endl;
    time.start();
    gradient_computation.compute_gradient_vector(forman_gradient,tree.get_root(),tree.get_mesh().get_domain(),tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] computing gradient vector field ");
    cerr << "[MEMORY] peak for computing the MultiVariate Forman Gradient: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;

    /// this extraction must be executed at the end..
    /// as it invalidates the critical d.s.
    cout<<"[NOTA] Extracting Critical Clusters"<<endl;
    Forman_Gradient_Features_Extractor features_extractor;
    features_extractor.extract_critical_clusters(tree.get_root(),tree.get_mesh(),forman_gradient,gradient_computation.get_critical_simplices(),
                                                 tree.get_subdivision(), cli.app_debug,cli.cache_size,get_path_without_file_extension(cli.mesh_path));
    cerr << "[MEMORY] peak for extracting the Critical Clusters: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;

    if(cli.app_debug == OUTPUT)
    {
        cout<<"[OUTPUT] saving the terrain with the extra fields"<<endl;
        Writer::write_mesh(get_path_without_file_extension(cli.mesh_path),"extra_fields",tree.get_mesh(),true);
    }
}

template<class T> void custom_execution(T& tree, cli_parameters &cli)
{
    if(cli.query_type != CUSTOM)
        return;

    Timer time = Timer();

    Topological_Queries tq;
    tq.batched_VT(tree.get_root(),tree.get_mesh().get_domain(),tree.get_mesh(),tree.get_subdivision(),cli.reindex);

    Forman_Gradient forman_gradient = Forman_Gradient(tree.get_mesh().get_triangles_num());
    Forman_Gradient_Computation gradient_computation = Forman_Gradient_Computation();

    cout<<"[NOTA] Computing the gradient field"<<endl;
    time.start();
    gradient_computation.initial_filtering(tree.get_mesh());
    time.stop();
    time.print_elapsed_time("[TIME] Initial filtering");
    time.start();
    gradient_computation.compute_gradient_vector(forman_gradient,tree.get_root(),tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] computing gradient vector field ");
    cerr << "[MEMORY] peak for computing the Forman Gradient: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;

    stringstream out;
    out << get_path_without_file_extension(cli.mesh_path);

    Forman_Gradient_Features_Extractor features_extractor;
    /// ---- MORSE INCIDENCE GRAPH COMPUTATION --- ///
    cout<<"[NOTA] Extracting the Critical Net."<<endl;
    time.start();
    features_extractor.extract_incidence_graph(tree.get_root(),tree.get_mesh(),forman_gradient,tree.get_subdivision(),cli.app_debug,cli.cache_size);
    time.stop();
    time.print_elapsed_time("[TIME] computing critical net ");
    cerr << "[MEMORY] peak for computing the critical net: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;

    if(cli.app_debug == OUTPUT)
    {
        features_extractor.print_stats();
        features_extractor.reset_stats();
        Writer_Morse::write_incidence_graph_VTK(out.str(),"mig", cli.v_per_leaf, features_extractor.get_incidence_graph(),tree.get_mesh(),
                                                cli.original_vertex_indices,cli.original_vertex_fields,cli.rever_to_original);
        features_extractor.reset_output_structures(tree.get_mesh());
    }
}

void custom_execution(PMRT_Tree& tree, cli_parameters &cli)
{
    if(cli.query_type != CUSTOM)
        return;

    Timer time = Timer();

    Topological_Queries tq;
    tq.batched_VT(tree.get_root(),tree.get_mesh().get_domain(),tree.get_mesh(),tree.get_subdivision(),cli.reindex);

    Forman_Gradient forman_gradient = Forman_Gradient(tree.get_mesh().get_triangles_num());
    Forman_Gradient_Computation gradient_computation = Forman_Gradient_Computation();

    cout<<"[NOTA] Computing the gradient field"<<endl;
    time.start();
    gradient_computation.initial_filtering(tree.get_mesh());
    time.stop();
    time.print_elapsed_time("[TIME] Initial filtering");
    time.start();
    gradient_computation.compute_gradient_vector(forman_gradient,tree.get_root(),tree.get_mesh().get_domain(),tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] computing gradient vector field ");
    cerr << "[MEMORY] peak for computing the Forman Gradient: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;

    stringstream out;
    out << get_path_without_file_extension(cli.mesh_path);

    Forman_Gradient_Features_Extractor features_extractor;
    /// ---- MORSE INCIDENCE GRAPH COMPUTATION --- ///
    cout<<"[NOTA] Extracting the Critical Net."<<endl;
    time.start();
    features_extractor.extract_incidence_graph(tree.get_root(),tree.get_mesh().get_domain(),tree.get_mesh(),forman_gradient,tree.get_subdivision(),cli.app_debug,cli.cache_size);
    time.stop();
    time.print_elapsed_time("[TIME] computing critical net ");
    cerr << "[MEMORY] peak for computing the critical net: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;

    if(cli.app_debug == OUTPUT)
    {
        features_extractor.print_stats();
        features_extractor.reset_stats();
        Writer_Morse::write_incidence_graph_VTK(out.str(),"mig", cli.v_per_leaf, features_extractor.get_incidence_graph(),tree.get_mesh(),
                                                cli.original_vertex_indices,cli.original_vertex_fields,cli.rever_to_original);
        features_extractor.reset_output_structures(tree.get_mesh());
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

template<class T> void load_tree_lite(T& tree, cli_parameters &cli)
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

template<class T> void SF_test(T& tree, cli_parameters &cli){
    Timer time = Timer();
    cli.app_debug=OUTPUT;
    //    Topological_Queries tq;
    //    tq.batched_VT(tree.get_root(),tree.get_mesh().get_domain(),tree.get_mesh(),tree.get_subdivision(),cli.reindex);

    load_terrain(tree,cli);

    Forman_Gradient forman_gradient = Forman_Gradient(tree.get_mesh().get_triangles_num());
    Forman_Gradient_Computation gradient_computation = Forman_Gradient_Computation();

    time.start();
    gradient_computation.initial_filtering_IA(tree.get_mesh());
    time.stop();
    time.print_elapsed_time("[TIME] Initial filtering ");

    load_tree_lite(tree,cli);

    //WARNING AFTER THIS THE FILTRATION ARRAY IN forman_gradient AND THE POSITION INDICES OF THE VERTICES OF THE TIN WILL NOT BE ALIGNED..
    // YOU HAVE TO USE original_vertex_indices FOR FETCHING THE CORRECT FILTRATION VALUE OF A VERTEX
    gradient_computation.reset_filtering(tree.get_mesh(),cli.original_vertex_indices, true);

    cout<<"[NOTA] Computing the gradient field"<<endl;
    time.start();
    gradient_computation.compute_gradient_vector(forman_gradient,tree.get_root(),tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] computing gradient vector field ");
    cerr << "[MEMORY] peak for computing the Forman Gradient: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;

    stringstream out;
    out << get_path_without_file_extension(cli.mesh_path);

    Forman_Gradient_Features_Extractor features_extractor;
    /// ---- MORSE INCIDENCE GRAPH COMPUTATION --- ///
    cout<<"[NOTA] Extracting the Critical Net."<<endl;
    time.start();
    features_extractor.extract_incidence_graph(tree.get_root(),tree.get_mesh(),forman_gradient,tree.get_subdivision(),cli.app_debug,cli.cache_size);
    time.stop();
    time.print_elapsed_time("[TIME] computing critical net ");
    cerr << "[MEMORY] peak for computing the critical net: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;

    // FOR DEBUG ONLY
    features_extractor.get_incidence_graph().print_stats(true);

    if(cli.app_debug == TIME_VERBOSE)
    {
        features_extractor.print_feature_extraction_time();
        features_extractor.reset_timer_variables();
    }

    if(cli.app_debug == OUTPUT)
    {
        features_extractor.print_stats();
        features_extractor.reset_stats();
        Writer_Morse::write_incidence_graph_VTK(out.str(),"mig", cli.v_per_leaf, features_extractor.get_incidence_graph(),tree.get_mesh(),
                                                cli.original_vertex_indices,cli.original_vertex_fields,cli.rever_to_original);
        features_extractor.reset_output_structures(tree.get_mesh());
    }
}


void SF_test(PMRT_Tree& tree, cli_parameters &cli){

    Timer time = Timer();
    cli.app_debug=OUTPUT;
    //    Topological_Queries tq;
    //    tq.batched_VT(tree.get_root(),tree.get_mesh().get_domain(),tree.get_mesh(),tree.get_subdivision(),cli.reindex);

    load_terrain(tree,cli);

    Forman_Gradient forman_gradient = Forman_Gradient(tree.get_mesh().get_triangles_num());
    Forman_Gradient_Computation gradient_computation = Forman_Gradient_Computation();

    time.start();
    gradient_computation.initial_filtering_IA(tree.get_mesh());
    time.stop();
    time.print_elapsed_time("[TIME] Initial filtering");

    load_tree_lite(tree,cli);

    //WARNING AFTER THIS THE FILTRATION ARRAY IN forman_gradient AND THE POSITION INDICES OF THE VERTICES OF THE TIN WILL NOT BE ALIGNED..
    // YOU HAVE TO USE original_vertex_indices FOR FETCHING THE CORRECT FILTRATION VALUE OF A VERTEX
    gradient_computation.reset_filtering(tree.get_mesh(),cli.original_vertex_indices, true);

    time.start();
    gradient_computation.compute_gradient_vector(forman_gradient,tree.get_root(),tree.get_mesh().get_domain(),tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] computing gradient vector field ");
    cerr << "[MEMORY] peak for computing the Forman Gradient: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;

    stringstream out;
    out << get_path_without_file_extension(cli.mesh_path);

    Forman_Gradient_Features_Extractor features_extractor;
    /// ---- MORSE INCIDENCE GRAPH COMPUTATION --- ///
    cout<<"[NOTA] Extracting the Critical Net."<<endl;
    time.start();
    features_extractor.extract_incidence_graph(tree.get_root(),tree.get_mesh().get_domain(),tree.get_mesh(),forman_gradient,tree.get_subdivision(),cli.app_debug,cli.cache_size);
    time.stop();
    time.print_elapsed_time("[TIME] computing critical net ");
    cerr << "[MEMORY] peak for computing the critical net: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;

    // FOR DEBUG ONLY
    features_extractor.get_incidence_graph().print_stats(true);

    if(cli.app_debug == TIME_VERBOSE)
    {
        features_extractor.print_feature_extraction_time();
        features_extractor.reset_timer_variables();
    }

    if(cli.app_debug == OUTPUT)
    {
        features_extractor.print_stats();
        features_extractor.reset_stats();
        Writer_Morse::write_incidence_graph_VTK(out.str(),"mig", cli.v_per_leaf, features_extractor.get_incidence_graph(),tree.get_mesh(),
                                                cli.original_vertex_indices,cli.original_vertex_fields,cli.rever_to_original);
        features_extractor.reset_output_structures(tree.get_mesh());
    }

}

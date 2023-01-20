#include "utilities/utility_functions.h"
using namespace utility_functions;

template<class T> void load_tree(T& tree, cli_parameters &cli);
template<class T> void compute_terrain_features(T& tree);
void compute_terrain_features(PMRT_Tree& tree);

int main(int , char** )
{
    cerr<<"[OBJECTIVE] this unit-test generates a quadtrees based on the P_Tree criterion. "
        <<"it computes some terrain features, such as the triangle and edge slopes and the critical points."<<endl;
    cerr<<"[TEST1] INPUT MESH eggs64z4.tri"<<endl;
    cli_parameters cli;
    cli.mesh_path = "../data/devil_0.tri";
	cli.division_type = QUAD;
    cli.crit_type = "pr";
    cli.v_per_leaf = 50;
    PRT_Tree ptree = PRT_Tree(cli.v_per_leaf,cli.division_type);
    cerr<<"[GENERATION] PR-T tree"<<endl;
    load_tree(ptree,cli);
    compute_terrain_features(ptree);
    cli.t_per_leaf = 100;
    PMRT_Tree rtree = PMRT_Tree(cli.t_per_leaf,cli.division_type);
    cerr<<"[GENERATION] PMR-T tree"<<endl;
    load_tree(rtree,cli);
    compute_terrain_features(rtree);

    return (EXIT_SUCCESS);
}

template<class T> void load_tree(T& tree, cli_parameters &cli)
{
    Timer time;
    if (!Reader::read_mesh(tree.get_mesh(), cli.mesh_path))
    {
        cout << "[ERROR] Loading mesh file. Execution Stopped." << endl;
        return;
    }

    stringstream base_info;
    base_info << cli.v_per_leaf << " " << cli.t_per_leaf << " " << cli.crit_type << " ";
    stringstream base;
    base << get_path_without_file_extension(cli.mesh_path);

    stringstream tree_info;
    tree_info << base_info.str() << "[TIME] Building ";
    time.start();
    tree.build_tree();
    time.stop();
    time.print_elapsed_time(tree_info.str());    
    
    time.start();
    Reindexer reindexer = Reindexer();
    reindexer.reindex_tree_and_mesh(tree,false,cli.original_vertex_indices,
                                    false,cli.original_triangle_indices);
    time.stop();
    time.print_elapsed_time("[TIME] Index and Mesh Reindexing ");    
}

template<class T> void compute_terrain_features(T& tree)
{
    Timer time;

    {
        cout<<"[TEST] Extract triangle slopes."<<endl;
        Slope_Extractor se;
        time.start();
        se.compute_triangles_slopes(tree.get_root(),tree.get_mesh(),tree.get_subdivision());
        time.stop();
        time.print_elapsed_time("[TIME] triangle-slopes computation: ");
        se.print_slopes_stats(tree.get_mesh().get_triangles_num());
        se.reset_stats();
    }

    // {
    //     cout<<"[TEST] New way to extract triangle slopes."<<endl;
    //     Slope_Extractor se;
    //     time.start();
    //     se.compute_triangles_slopes_new(tree.get_root(),tree.get_mesh(),tree.get_subdivision());
    //     time.stop();
    //     time.print_elapsed_time("[TIME] triangle-slopes computation: ");
    //     se.print_slopes_stats(tree.get_mesh().get_triangles_num());
    //     se.reset_stats();
    // }
    {
        cout<<"[TEST] Extract edge slopes."<<endl;
        Slope_Extractor se;
        time.start();
        se.compute_edges_slopes(tree.get_root(),tree.get_mesh(),tree.get_subdivision());
        time.stop();
        time.print_elapsed_time("[TIME] edge-slopes computation: ");
        se.print_slopes_stats();
        se.reset_stats();
    }
    
    {
        cout<<"[TEST] Extract triangle aspects."<<endl;
        Aspect aspect;
        time.start();
        aspect.compute_triangles_aspects(tree.get_root(),tree.get_mesh(),tree.get_subdivision());
        time.stop();
        time.print_elapsed_time("[TIME] triangle-aspects computation: ");
        aspect.print_aspects_stats(tree.get_mesh().get_triangles_num());
        aspect.reset_stats();
    }
    {
        cout<<"[TEST] Extract critical points."<<endl;
        Critical_Points_Extractor cpe = Critical_Points_Extractor();
        time.start();
        cpe.compute_critical_points(tree.get_root(),tree.get_mesh(),tree.get_subdivision());
        time.stop();
        time.print_elapsed_time("[TIME] critical points extraction: ");
        cpe.print_stats();
    }
}

void compute_terrain_features(PMRT_Tree& tree)
{
    Timer time;

    {
        cout<<"[TEST] Extract triangle slopes."<<endl;
        Slope_Extractor se;
        time.start();
        se.compute_triangles_slopes(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_subdivision());
        time.stop();
        time.print_elapsed_time("[TIME] triangle-slopes computation: ");
        se.print_slopes_stats(tree.get_mesh().get_triangles_num());
        se.reset_stats();
    }

    // {
    //     cout<<"[TEST] New way to Extract triangle slopes."<<endl;
    //     Slope_Extractor se;
    //     time.start();
    //     se.compute_triangles_slopes_new(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_subdivision());
    //     time.stop();
    //     time.print_elapsed_time("[TIME] triangle-slopes computation: ");
    //     se.print_slopes_stats(tree.get_mesh().get_triangles_num());
    //     se.reset_stats();
    // }
    
    
    {
        cout<<"[TEST] Extract edge slopes."<<endl;
        Slope_Extractor se;
        time.start();
        se.compute_edges_slopes(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_subdivision());
        time.stop();
        time.print_elapsed_time("[TIME] edge-slopes computation: ");
        se.print_slopes_stats();
        se.reset_stats();
    }
    
    {
        cout<<"[TEST] Extract triangle aspects."<<endl;
        Aspect aspect;
        time.start();
        aspect.compute_triangles_aspects(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_subdivision());
        time.stop();
        time.print_elapsed_time("[TIME] triangle-aspects computation: ");
        aspect.print_aspects_stats(tree.get_mesh().get_triangles_num());
        aspect.reset_stats();
    }
    {
        cout<<"[TEST] Extract critical points."<<endl;
        Critical_Points_Extractor cpe = Critical_Points_Extractor();
        time.start();
        cpe.compute_critical_points(tree.get_root(),tree.get_mesh().get_domain(),tree.get_mesh(),tree.get_subdivision());
        time.stop();
        time.print_elapsed_time("[TIME] critical points extraction: ");
        cpe.print_stats();
    }
}

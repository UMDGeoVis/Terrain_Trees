#include "utilities/utility_functions.h"

using namespace utility_functions;

void load_tree(PRT_Tree& tree, cli_parameters &cli);
void calc_stats(PRT_Tree& tree, cli_parameters &cli, Mesh &mesh, string& csv_path);



int main(int argc, char** argv)
{
	cli_parameters cli;
	cli.mesh_path = argv[1];
	cli.debug_mode = false;
	cerr<<"[OBJECTIVE] this unit-test reads a PR-quadtree on the input TIN dataset and "
	    <<"a corresponding csv file which stores the contracted edge number for each leaf node."<<endl;
		
    cerr<<"[NOTA] all the generated files are saved in the 'data' folder"<<endl;
	
	cli.division_type = QUAD;
    cli.crit_type = "pr";
    cli.v_per_leaf = atoi(argv[2]);
    PRT_Tree ptree = PRT_Tree(cli.v_per_leaf,cli.division_type);
    load_tree(ptree,cli);


    if(strcmp(argv[3], "-d") == 0){
        string csv_path = argv[4];
        stringstream output_name;
        output_name << get_path_without_file_extension(cli.mesh_path);
        output_name << "_edge_distribution.vtk";
        Writer::write_edge_cost_distribution(ptree.get_root(), ptree.get_subdivision(), output_name.str(), ptree.get_mesh(), csv_path);
    }
    else if(strcmp(argv[3], "-s") == 0){
        cout<<"[NOTA]Checking Delaunay property"<<endl;
        ptree.init_leaves_list(ptree.get_root()); 
        Contraction_Simplifier simplifier;
        simplifier.check_delaunay(ptree, ptree.get_mesh());
        cout<<"[NOTA]Compute compactness"<<endl;
        simplifier.compute_compactness(ptree, ptree.get_mesh());
    }

    return (EXIT_SUCCESS);
}

void load_tree(PRT_Tree& tree, cli_parameters &cli)
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
        out << "_v_" << cli.v_per_leaf << "_.tree";
 
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

    }

    cerr << "[MEMORY] peak for encoding the Terrain tree: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;


    stringstream out2;
    out2 << base.str();
    out2 << "_" << SpatialDecType2string(cli.division_type) << "_" << cli.crit_type << "_v_" << cli.v_per_leaf << "_tree.vtk";
    


    time.start();
    Reindexer reindexer = Reindexer();
    reindexer.reindex_tree_and_mesh(tree,false,cli.original_vertex_indices,
                                    false,cli.original_triangle_indices);
    time.stop();
    time.print_elapsed_time("[TIME] Index and Mesh Reindexing ");
         cerr << "[MEMORY] peak for Index and Mesh Reindexing: " <<
        to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
}  


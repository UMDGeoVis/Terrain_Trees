#include "utilities/utility_functions.h"

using namespace utility_functions;

void load_tree(PRT_Tree& tree, cli_parameters &cli);
void interpolate_vertices(PRT_Tree& tree, cli_parameters& cli, vector<Vertex>& input_vertices, vector<coord_type>& distances);

int main(int argc, char** argv)
{
	cli_parameters cli;
	cli.mesh_path = argv[1]; // simplified mesh
	cli.debug_mode = false;
	cerr<<"[OBJECTIVE] this unit-test reads the simplified mesh and the original mesh. For each vertex in the original mesh, it calculates its vertical distance to the  simplified mesh."
	    <<"The output is the maximum of all vertice distances. "<<endl;
		
	
	cli.division_type = QUAD;
    cli.crit_type = "pr";
    cli.v_per_leaf = atoi(argv[2]);
    PRT_Tree ptree = PRT_Tree(cli.v_per_leaf,cli.division_type);
    load_tree(ptree,cli);

    ptree.init_leaves_list(ptree.get_root()); 

    string orig_mesh_path = argv[3];
    stringstream output_name;
    Mesh orig_mesh;
    if (!Reader::read_vertices(orig_mesh, orig_mesh_path))
    {
        cout << "[ERROR] Loading original mesh file. Execution Stopped." << endl;
        return -1;
    }
    vector<coord_type> distances(orig_mesh.get_vertices_num(), 0);
    vector<Vertex> vertices = orig_mesh.get_vertices_array();
    
    Distance_Calculator distance_calculator;
    distance_calculator.vertical_distance(ptree, vertices);
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

    stringstream out;
    out << base.str() << "_" << SpatialDecType2string(cli.division_type) << "_" << cli.crit_type;
    out << "_v_" << cli.v_per_leaf << "_.tree";

    cli.tree_path=out.str();


    cerr << "[GENERATION] tree from triangle mesh" << endl;
    time.start();
    tree.build_tree();
    time.stop();
    time.print_elapsed_time(tree_info.str());
    // Writer::write_tree(out.str(), tree.get_root(), tree.get_subdivision());

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





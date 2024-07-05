#include "utilities/utility_functions.h"

using namespace utility_functions;

template<class T> void load_tree(T& tree, cli_parameters &cli, coord_type& mode_to_correct, Point& origin);
template<class T> void compute_interpolated_elevations(T& tree, cli_parameters &cli, const Point& origin, string& point_file_path);
Point standarize_input_mesh(Mesh& mesh, coord_type& mode_to_correct);
// void reverse_mesh_coordinates(Mesh& mesh, const Point& origin);
void output_triangle_area(Mesh& mesh, string path);
void standarize_input_points(vector<Point>& query_points, const Point& reference_origin);
void reverse_input(vector<Point>& query_points,  const Point& reference_origin);

int main(int argc , char** argv)
{
	cli_parameters cli;
    cli.mesh_path = argv[1];
    
    cerr<<"[OBJECTIVE] this unit-test generates a quadtrees based on the PR-T tree criterion. "
	    <<"then, it computes the roughness of all the vertices in the mesh. "<<endl;
		
    cerr<<"[NOTA] all the generated files are saved in the 'data' folder"<<endl;
    
    cli.division_type = QUAD;
    cli.crit_type = "pr";
    cli.v_per_leaf = atoi(argv[2]);    
    coord_type mode_to_correct = atof(argv[3]);
    string point_file_path = argv[4];

    PRT_Tree ptree = PRT_Tree(cli.v_per_leaf,cli.division_type);
    cerr<<"[GENERATION] PR-T tree"<<endl;
    Point origin; 
    load_tree(ptree,cli, mode_to_correct, origin);
    compute_interpolated_elevations(ptree,cli, origin, point_file_path);
    // output_triangle_area(ptree.get_mesh(), get_path_without_file_extension(cli.mesh_path));
    return (EXIT_SUCCESS);
}

template<class T> void load_tree(T& tree, cli_parameters &cli, coord_type& mode_to_correct, Point& origin)
{
    Timer time;
    if (!Reader::read_mesh(tree.get_mesh(), cli.mesh_path))
    {
        cout << "[ERROR] Loading mesh file. Execution Stopped." << endl;
        return;
    }

    origin = standarize_input_mesh(tree.get_mesh(), mode_to_correct);

    stringstream base_info;
    base_info << cli.v_per_leaf << " " << cli.t_per_leaf << " " << cli.crit_type << " ";
    stringstream base;
    base << get_path_without_file_extension(cli.mesh_path);

    stringstream tree_info;
    tree_info << base_info.str() << "[TIME] Building ";
    stringstream out;
    out << base.str() << "_" << SpatialDecType2string(cli.division_type) << "_" << cli.crit_type;
    if (cli.crit_type == "pr")
        out << "_v_" << cli.v_per_leaf << "_.tree";
    else if (cli.crit_type == "pm")
        out << "_v_" << cli.v_per_leaf << "_t_" << cli.t_per_leaf << "_.tree";
    else if (cli.crit_type == "pmr")
        out << "_t_" << cli.t_per_leaf << "_.tree";
    cli.tree_path=out.str();
    time.start();
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

    time.stop();
    time.print_elapsed_time(tree_info.str());

    stringstream out2;
    out2 << base.str();
    out2 << "_" << SpatialDecType2string(cli.division_type) << "_" << cli.crit_type << "_v_" << cli.v_per_leaf << "_tree.vtk";
    
    // Writer::write_tree_VTK(out2.str(),tree.get_root(),tree.get_subdivision(),tree.get_mesh());
    // Writer::write_mesh_VTK(base.str(),tree.get_mesh());        
    time.start();
    Reindexer reindexer = Reindexer();
    reindexer.reindex_tree_and_mesh(tree,false,cli.original_vertex_indices,
                                    false,cli.original_triangle_indices);
    time.stop();
    time.print_elapsed_time("[TIME] Index and Mesh Reindexing ");
         cerr << "[MEMORY] peak for Index and Mesh Reindexing: " <<
        to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
}  


Point standarize_input_mesh(Mesh& mesh, coord_type& mode_to_correct)
{
    if(mesh.get_vertices_num() == 0) return Point(0, 0);
    Point origin = mesh.get_vertex(1);

    cout<< mesh.get_vertex(1)<<endl;
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

void standarize_input_points(vector<Point>& query_points, const Point& reference_origin)
{
    for(int i = 0; i < query_points.size(); i++){
        Point old = query_points[i];
        old.set_c(0, old.get_x() - reference_origin.get_x());
        old.set_c(1, old.get_y() - reference_origin.get_y());
        query_points[i] = old;
    }

}

template<class T> void compute_interpolated_elevations(T& tree, cli_parameters &cli, const Point& origin, string& point_file_path)
{
    stringstream out;
    out << get_path_without_file_extension(cli.mesh_path);
    out << "_oib_icesat2_";
    vector<Point> query_points;
    Reader::read_queries(query_points, point_file_path);
    std::cout <<"Read input point list: " << query_points.size() << " query points in total"<<endl;
    standarize_input_points(query_points, origin);
    vector<coord_type> elevations(query_points.size());
    vector<bool> intersect(query_points.size(), false);
    Spatial_Queries sq;
    #pragma omp paralllel for
    for(unsigned int i = 0; i < query_points.size(); i++){
        intersect[i] = sq.exec_point_interpolation(tree, query_points[i], tree.get_mesh(), tree.get_subdivision(), elevations[i]);
    }
    reverse_input(query_points, origin);
    Writer::write_interpolation_results(out.str(), elevations, query_points, intersect);
}

void reverse_input(vector<Point>& query_points, const Point& origin)
{
    for(int i = 0; i < query_points.size(); i++){
        Point old = query_points[i];
        old.set_c(0, old.get_x() + origin.get_x());
        old.set_c(1, old.get_y() + origin.get_y());
        query_points[i] = old;
    }
}

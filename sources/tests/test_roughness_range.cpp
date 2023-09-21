#include "utilities/utility_functions.h"

using namespace utility_functions;

template<class T> void load_tree(T& tree, cli_parameters &cli, coord_type& radius);
template<class T> void compute_roughness(T& tree, cli_parameters &cli, coord_type& radius);
void standarize_input(Mesh& mesh, coord_type& radius);
void output_triangle_area(Mesh& mesh, string path);
coord_type compute_area(Triangle& t, Mesh& mesh){
    coord_type area;
    
    Vertex &vi = mesh.get_vertex(t.TV(0));
    Vertex &vj = mesh.get_vertex(t.TV(1));
    Vertex &vk = mesh.get_vertex(t.TV(2));

    dvect ki={vi.get_x()-vk.get_x(),vi.get_y()-vk.get_y()};
    dvect ij={vj.get_x()-vi.get_x(),vj.get_y()-vi.get_y()};
        // compute the vectors of the edges after rotate by 90 degrees
    dvect ki_vert= { vk.get_y()-vi.get_y() , vi.get_x()-vk.get_x()};
    dvect ij_vert= {vi.get_y()-vj.get_y() , vj.get_x()-vi.get_x()};
    // compute the area of the triangle 
    area=abs(0.5*(ki[0]*ij[1]-ki[1]*ij[0]));
    return area;
}
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
    coord_type radius = atof(argv[3]);
    coord_type mode_to_correct = atof(argv[4]);
    PRT_Tree ptree = PRT_Tree(cli.v_per_leaf,cli.division_type);
    cerr<<"[GENERATION] PR-T tree"<<endl;
    load_tree(ptree,cli, mode_to_correct);
    compute_roughness(ptree,cli, radius);
    // output_triangle_area(ptree.get_mesh(), get_path_without_file_extension(cli.mesh_path));

    return (EXIT_SUCCESS);
}

template<class T> void load_tree(T& tree, cli_parameters &cli, coord_type& mode_to_correct)
{
    Timer time;
    if (!Reader::read_mesh(tree.get_mesh(), cli.mesh_path))
    {
        cout << "[ERROR] Loading mesh file. Execution Stopped." << endl;
        return;
    }

    standarize_input(tree.get_mesh(), mode_to_correct);

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

template<class T> void compute_roughness(T& tree, cli_parameters &cli, coord_type& radius)
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
    roughness.print_roughness_stats(tree.get_mesh(),tree.get_mesh().get_vertex(1).get_fields_num()-1);

    // Writer::write_roughness_txt(out.str(),tree.get_mesh(),tree.get_mesh().get_vertex(1).get_fields_num() - 2, radius);
    Writer::write_elevation_txt(out.str(),tree.get_mesh());

    Writer::write_mesh_roughness_VTK(out.str(),tree.get_mesh(),tree.get_mesh().get_vertex(1).get_fields_num() - 2, radius, true);
    Writer::write_mesh_roughness_VTK(out.str(),tree.get_mesh(),tree.get_mesh().get_vertex(1).get_fields_num() - 2, radius);
   }

void standarize_input(Mesh& mesh, coord_type& mode_to_correct)
{
    if(mesh.get_vertices_num() == 0) return;
    Point origin = mesh.get_vertex(1);
    // cout<< mesh.get_vertex(0)<<endl;

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
}

void output_triangle_area(Mesh& mesh, string path)
{
    
    dvect areas(mesh.get_triangles_num(), 0);
    #pragma omp parallel for
    for(int i = 1; i <= mesh.get_triangles_num(); i++){
        Triangle t = mesh.get_triangle(i);
        areas[i - 1] = compute_area(t, mesh);
    }
    Writer::write_tri_area_VTK(path, mesh, areas);
}
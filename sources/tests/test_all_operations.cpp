/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "utilities/utility_functions.h"
#include "gradient/Gradient.h"
using namespace utility_functions;
template<class T> void load_tree(T& tree, cli_parameters &cli);


template<class T> void compute_roughness(T& tree, cli_parameters &cli);
void compute_roughness(PMRT_Tree& tree, cli_parameters &cli);

template<class T> void compute_curvature(T& tree, cli_parameters &cli);
void compute_curvature(PMRT_Tree& tree, cli_parameters &cli);

template<class T> void compute_multifield(T& tree, cli_parameters &cli, string mode);
void compute_multifield(PMRT_Tree& tree, cli_parameters &cli,string mode);

template<class T> void compute_terrain_features(T& tree);
void compute_terrain_features(PMRT_Tree& tree);

template<class T> void relation_extraction(T& tree, cli_parameters &cli);
void relation_extraction(PMRT_Tree& tree, cli_parameters &cli);

template<class T> void field_output(T& tree, cli_parameters &cli);
void field_output(PMRT_Tree& tree, cli_parameters &cli);


int main(int argc, char** argv )
{
    cli_parameters cli;
    cli.mesh_path = argv[1];

    cerr<<"[OBJECTIVE] This unit is designed to put all the operations(except for the morse related operations) included in the experiments in one executable file.";

    cerr<<"[NOTA] all the generated files are saved in the 'data' folder"<<endl;

    cli.division_type = QUAD;

    cli.v_per_leaf = atoi(argv[2]);
    cli.t_per_leaf= atoi(argv[3]);
    //string mode="rRGB";
    string operation = argv[5];
    string mode="All";
    // int factor=atoi(argv[3]);
    string tree_type=argv[4];


    if(tree_type=="PR")
    {
    cli.crit_type = "pr";
    PRT_Tree ptree = PRT_Tree(cli.v_per_leaf,cli.division_type);
    cerr<<"[GENERATION] PR-T tree"<<endl;
    load_tree(ptree,cli);
    if(operation=="multifield"){
    compute_roughness(ptree,cli);
    compute_multifield(ptree,cli,mode);}
    else if (operation=="roughness"){
        compute_roughness(ptree,cli);
    }
    else if(operation=="slope"){
          compute_terrain_features(ptree);
    }
    else if (operation=="curvature"){
            compute_curvature(ptree,cli);
    }
    else if(operation=="relation"){
        relation_extraction(ptree,cli);
    }

    }
    else if(tree_type=="PMR"){
      cli.crit_type = "pmr";
    cli.v_per_leaf = -1;

    PMRT_Tree rtree = PMRT_Tree(cli.t_per_leaf,cli.division_type);
    cerr<<"[GENERATION] PMR-T tree"<<endl;
    load_tree(rtree,cli);
    if(operation=="multifield"){
    compute_roughness(rtree,cli);
    compute_multifield(rtree,cli,mode);}
    else if (operation=="roughness"){
        compute_roughness(rtree,cli);
    }
    else if(operation=="slope"){
        compute_terrain_features(rtree);
    }
    else if (operation=="curvature"){
        compute_curvature(rtree,cli);
    }
        else if(operation=="relation"){
        relation_extraction(rtree,cli);
    }
    }
    else if(tree_type=="PM"){
            cli.crit_type = "pm";
    PMT_Tree pttree = PMT_Tree(cli.v_per_leaf,cli.t_per_leaf,cli.division_type);
    cerr<<"[GENERATION] PM_Tree"<<endl;
    load_tree(pttree,cli);
    if(operation=="multifield")
    {compute_roughness(pttree,cli);
    compute_multifield(pttree,cli,mode);}
    else if (operation=="roughness"){
        compute_roughness(pttree,cli);
    }
    else if(operation=="slope"){
          compute_terrain_features(pttree);
    }
    else if (operation=="curvature"){
            compute_curvature(pttree,cli);
    }
        else if(operation=="relation"){
        relation_extraction(pttree,cli);
    }
    }




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
        cerr << "[ERROR] Loading .tree file. Execution Stopped." << endl;
            tree.build_tree();
              time.start();
             Writer::write_tree(out.str(), tree.get_root(), tree.get_subdivision());
  time.stop();
        time.print_elapsed_time(tree_info.str());
    }






    stringstream out2;
    out2 << base.str();
    out2 << "_" << SpatialDecType2string(cli.division_type) << "_" << cli.crit_type << "_v_" << cli.v_per_leaf << "_tree.vtk";

    Writer::write_tree_VTK(out2.str(),tree.get_root(),tree.get_subdivision(),tree.get_mesh());
    Writer::write_mesh_VTK(base.str(),tree.get_mesh());

    time.start();
    Reindexer reindexer = Reindexer();
    reindexer.reindex_tree_and_mesh(tree,false,cli.original_vertex_indices,
                                    false,cli.original_triangle_indices);
    time.stop();
    time.print_elapsed_time("[TIME] Index and Mesh Reindexing ");

          cerr << "[MEMORY] peak for indexing: " <<
        to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
}




template<class T> void compute_roughness(T& tree, cli_parameters &cli)
{
    stringstream out;
    out << get_path_without_file_extension(cli.mesh_path);
       Timer time;

    cout<<"[NOTA]Border checking"<<endl;
    Border_Checker border_checker=Border_Checker();
           time.start();
    border_checker.compute_borders(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_subdivision());
           time.stop();
           time.print_elapsed_time("[TIME] Border Checking time:");

    cout<<"[NOTA] compute roughness"<<endl;

    Roughness roughness = Roughness(tree.get_mesh());
    time.start();
    roughness.compute(tree.get_root(),tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] roughness computation: ");
     cerr << "[MEMORY] peak for computing Roughness: " <<
        to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
        //    roughness.store_result(tree.get_mesh());
    roughness.print_roughness_stats(tree.get_mesh(),tree.get_mesh().get_vertex(1).get_fields_num()-1);
 
 Writer::write_mesh_roughness_VTK(out.str(),tree.get_mesh(),tree.get_mesh().get_vertex(1).get_fields_num()-1);
   }

void compute_roughness(PMRT_Tree& tree, cli_parameters &cli)
{
      stringstream out;
    out << get_path_without_file_extension(cli.mesh_path);

    Timer time;
     cout<<"[NOTA]Border checking"<<endl;
    Border_Checker border_checker=Border_Checker();
    time.start();
    border_checker.compute_borders(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] Border Checking time:");
    cout<<"[NOTA] compute roughness"<<endl;
    Roughness roughness = Roughness(tree.get_mesh());
    time.start();
    roughness.compute(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] roughness computation: ");
         cerr << "[MEMORY] peak for computing Roughness: " <<
        to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
            //    roughness.store_result(tree.get_mesh());
    roughness.print_roughness_stats(tree.get_mesh(),tree.get_mesh().get_vertex(1).get_fields_num()-1);
    }


template<class T> void compute_curvature(T& tree, cli_parameters &cli)
{
    stringstream out;
    out << get_path_without_file_extension(cli.mesh_path);

     cout<<"[NOTA] compute concentrated-curvature"<<endl;
    Timer time;
    Concentrated_Curvature curvature = Concentrated_Curvature(CONCENTRATED);
    time.start();
    curvature.compute_borders(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] border computation: ");
    time.start();
    curvature.compute(tree.get_root(),tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] concentrated-curvature computation: ");
    curvature.print_curvature_stats(tree.get_mesh(),tree.get_mesh().get_vertex(1).get_fields_num()-1);
    Writer::write_mesh_curvature_VTK(out.str(),tree.get_mesh(),"concentrated",tree.get_mesh().get_vertex(1).get_fields_num()-1);

}

void compute_curvature(PMRT_Tree& tree, cli_parameters &cli)
{
    stringstream out;
    out << get_path_without_file_extension(cli.mesh_path);


    cout<<"[NOTA] compute concentrated-curvature"<<endl;
    Timer time;
    Concentrated_Curvature curvature = Concentrated_Curvature(CONCENTRATED);
    time.start();
    curvature.compute_borders(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] border computation: ");
    time.start();
    curvature.compute(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] concentrated-curvature computation: ");
    curvature.print_curvature_stats(tree.get_mesh(),tree.get_mesh().get_vertex(1).get_fields_num()-1);
    Writer::write_mesh_curvature_VTK(out.str(),tree.get_mesh(),"concentrated",tree.get_mesh().get_vertex(1).get_fields_num()-1);

}



template<class T> void compute_terrain_features(T& tree)
{
    Timer time;


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
        cout<<"[TEST] Extract critical points."<<endl;
        Critical_Points_Extractor cpe = Critical_Points_Extractor();
        time.start();
        cpe.compute_critical_points(tree.get_root(),tree.get_mesh(),tree.get_subdivision());
        time.stop();
        time.print_elapsed_time("[TIME] critical points extraction: ");
        cpe.print_stats();
    }

    {
        cout<<"[TEST] New way to extract triangle slopes."<<endl;
        Slope_Extractor se;
        time.start();
        se.compute_triangles_slopes(tree.get_root(),tree.get_mesh(),tree.get_subdivision());
        time.stop();
        time.print_elapsed_time("[TIME] triangle-slopes computation: ");
        se.print_slopes_stats(tree.get_mesh().get_triangles_num());
        se.reset_stats();
    }
}

void compute_terrain_features(PMRT_Tree& tree)
{
    Timer time;


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
        cout<<"[TEST] Extract critical points."<<endl;
        Critical_Points_Extractor cpe = Critical_Points_Extractor();
        time.start();
        cpe.compute_critical_points(tree.get_root(),tree.get_mesh().get_domain(),tree.get_mesh(),tree.get_subdivision());
        time.stop();
        time.print_elapsed_time("[TIME] critical points extraction: ");
        cpe.print_stats();
    }

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

}




template<class T> void compute_multifield(T& tree, cli_parameters &cli,string mode)
{
    stringstream out;
    out << get_path_without_file_extension(cli.mesh_path);
    cout<<"[NOTA] compute multifield-based value"<<endl;
    Timer time;

 //   Built-in parameters
    int field_pos=0;  // The starting position for field values, 0 means the first field value is z value. 
    int factor=100;

    Gradient gradient(mode);
 //   
    time.start();
    gradient.compute_field_stats(tree.get_mesh(),factor);
    gradient.multi_field(tree.get_root(),tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] multi-field computation: ");

    cerr << "[MEMORY] peak for computing Multi field measure: " <<
    to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;

    gradient.print_multifield_stats(tree.get_mesh(),tree.get_mesh().get_vertex(1).get_fields_num()-1);

    Writer::write_mesh_multifield_VTK(out.str(),tree.get_mesh(),tree.get_mesh().get_vertex(1).get_fields_num()-1,mode);
  // Writer::write_field_csv(out.str(),tree.get_mesh());

}

void compute_multifield(PMRT_Tree& tree, cli_parameters &cli,string mode)
{
    stringstream out;
    out << get_path_without_file_extension(cli.mesh_path);
    Timer time;
    time.start();
     //   Built-in parameters
    int field_pos=0;
    int factor=100;


    Gradient gradient(mode);

    gradient.compute_field_stats(tree.get_mesh(),factor);
    cout<<"[NOTA] compute multifield-based value"<<endl;
    gradient.multi_field(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_subdivision());

    time.stop();
    time.print_elapsed_time("[TIME] multi-field computation: ");

    cerr << "[MEMORY] peak for computing Multi field measure: " <<
    to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
    gradient.print_multifield_stats(tree.get_mesh(),tree.get_mesh().get_vertex(1).get_fields_num()-1);
    Writer::write_mesh_multifield_VTK(out.str(),tree.get_mesh(),tree.get_mesh().get_vertex(1).get_fields_num()-1,mode);
    //Writer::write_field_csv(out.str(),tree.get_mesh());
    }


template<class T> void relation_extraction(T& tree, cli_parameters &cli)
{
    stringstream out;
    out << get_path_without_file_extension(cli.mesh_path);
    Timer time;
    Gradient gradient;
    cout<<"[NOTA] extract VT relation"<<endl;
    time.start();
    gradient.VT_relation(tree.get_root(),tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] VT relation extraction time:");

    cout<<"[NOTA] extract VV relation"<<endl;
    time.start();
    gradient.VV_relation(tree.get_root(),tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] VV relation extraction time:");

   }

void relation_extraction(PMRT_Tree& tree, cli_parameters &cli)
{
    stringstream out;
    out << get_path_without_file_extension(cli.mesh_path);

    Timer time;
    Gradient gradient;
    cout<<"[NOTA] extract VT relation"<<endl;
    time.start();
    gradient.VT_relation(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] VT relation extraction time:");

    cout<<"[NOTA] extract VV relation"<<endl;
        time.start();
    gradient.VV_relation(tree.get_root(),tree.get_mesh().get_domain(),0,tree.get_mesh(),tree.get_subdivision());
    time.stop();
    time.print_elapsed_time("[TIME] VV relation extraction time:");
 
    }

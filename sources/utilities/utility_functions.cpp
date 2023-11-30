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

#include "utility_functions.h"

namespace utility_functions
{

string SpatialDecType2string(SpatialDecType sdt)
{
    if(sdt == KD)
        return "kd";
    else if(sdt == QUAD)
        return "quad";
    else
        return DEFAULT;
}

void setParameters(cli_parameters &variables)
{
    vector<string> tokens;
    tokenize(variables.tree_path, tokens, "_");

    for (unsigned int i = 0; i < tokens.size(); i++)
    {
        if (tokens.at(i) == "pm" || tokens.at(i) == "pr" || tokens.at(i) == "pmr")
            variables.crit_type = tokens.at(i);
        if (tokens.at(i) == "quad")
            variables.division_type = QUAD;
        if (tokens.at(i) == "kd")
            variables.division_type = KD;
        if (tokens.at(i) == "v")
            variables.v_per_leaf = atoi(tokens.at(i + 1).c_str());
        if (tokens.at(i) == "t")
            variables.t_per_leaf = atoi(tokens.at(i + 1).c_str());
    }
}

bool checkParameters(cli_parameters &variables)
{
    if (variables.crit_type == DEFAULT || variables.division_type == UNSET)
    {
        cout << "[ERROR] initializing criterion or division type. Execution Stopped." << endl;
        print_usage();
        return false;
    }

    if (variables.crit_type == "pm" && (variables.v_per_leaf == DEFAULT_X_PER_LEAF || variables.t_per_leaf == DEFAULT_X_PER_LEAF))
    {
        cout << "[ERROR] initializing vertices_per_leaf or triangles_per_leaf. Execution Stopped." << endl;
        print_usage();
        return false;
    }
    if (variables.crit_type == "pr" && variables.v_per_leaf == DEFAULT_X_PER_LEAF)
    {
        cout << "[ERROR] initializing vertices_per_leaf. Execution Stopped." << endl;
        print_usage();
        return false;
    }
    if ((variables.crit_type == "pmr") && variables.t_per_leaf == DEFAULT_X_PER_LEAF)
    {
        cout << "[ERROR] initializing triangles_per_leaf. Execution Stopped." << endl;
        print_usage();
        return false;
    }
    return true;
}

int read_arguments(int argc, char** argv, cli_parameters &variables)
{
    string trash;
    variables.exe_name = argv[0];
    variables.exe_name = strip_path(variables.exe_name);

    for(int i=1; i<(argc); i++)
    {
        char* tag = argv[i];
        if(strcmp(tag, "-i") == 0)
        {
            variables.mesh_path = argv[i+1];
            i++;
        }
        else if(strcmp(tag, "-f") == 0)
        {
            variables.tree_path = argv[i+1];
            variables.isTreeFile = true;
            i++;
        }
        else if(strcmp(tag, "-d") == 0)
        {
            if(strcmp(argv[i+1],"quad") == 0)
                variables.division_type = QUAD;
            else if(strcmp(argv[i+1],"kd") == 0)
                variables.division_type = KD;
            i++;
        }
        else if(strcmp(tag, "-c") == 0)
        {
            variables.crit_type = argv[i+1];
            i++;
        }
        else if(strcmp(tag, "-v") == 0)
        {
            variables.v_per_leaf = atoi(argv[i+1]);
            if (variables.v_per_leaf < 1) {
                cerr << "[ERROR] the limit of vertices per leaf must be greater than 0" << endl;
                return -1;
            }
            i++;
        }
        else if(strcmp(tag, "-t") == 0)
        {
            variables.t_per_leaf = atoi(argv[i+1]);
            if (variables.t_per_leaf < 1) {
                cerr << "[ERROR] the limit of triangles per leaf must be greater than 0" << endl;
                return -1;
            }
            i++;
        }
        else if(strcmp(tag, "-s") == 0)
        {
            variables.is_index = true;
        }
        else if(strcmp(tag, "-noR") == 0)
        {
            variables.reindex = false;
        }
        else if(strcmp(tag, "-q") == 0)
        {
            trash = argv[i+1];
            vector<string> tok;
            tokenize(trash,tok,"-");
            if(tok.size()==1)
            {               
                /// TERRAIN TOPOLOGICAL FEATURES EXTRACTION
                if(tok[0]=="batch")
                    variables.query_type = BATCH;
                /// CURVATURE EXTRACTION
                else if(tok[0]=="concurv")
                    variables.query_type = CONCENTRATED_CURVATURE;
                else if(tok[0]=="mccurv")
                    variables.query_type = MEAN_CCURVATURE;
                else if(tok[0]=="gccurv")
                    variables.query_type = GAUSS_CCURVATURE;
                else if(tok[0]=="morse")
                    variables.query_type = MORSE_ANALYSIS;
                // else if(tok[0]=="simpl")
                //     variables.query_type = LOCAL_MORSE_SIMPLIFICATION;
                else if(tok[0]=="gsimpl")
                     variables.query_type = GLOBAL_MORSE_SIMPLIFICATION;
                else if(tok[0]=="multiv")
                    variables.query_type = MULTIVARIATE_MORSE_ANALYSIS;
                else if(tok[0]=="filter")
                    variables.query_type = FILTER;
                else if(tok[0]=="crit")
                    variables.query_type = CRITICAL_POINTS;
                else if(tok[0]=="slopes")
                    variables.query_type = SLOPES;
                else if(tok[0]=="rough")
                    variables.query_type = ROUGHNESS;
                else if(tok[0]=="multifield")
                    variables.query_type = MULTIFIELD;
                /// CUSTOM COMMAND --- TEMPORARY COMMAND FOR JOURNAL PAPER
                else if(tok[0]=="custom")
                    variables.query_type = CUSTOM;
            }
            else if(tok.size()<2)
                cerr<<"[ERROR] [-q argument] when reading arguments"<<endl;
            else
            {
                if(tok[0] == "wvt")
                    variables.query_type = WINDVT;
                else if(tok[0] == "wtt")
                    variables.query_type = WINDTT;
                else if(tok[0] == "point")
                    variables.query_type = POINT;
                else if(tok[0] == "box")
                    variables.query_type = BOX;
                else if(tok[0]=="simpl")
                    variables.query_type = LOCAL_MORSE_SIMPLIFICATION;
                else if(tok[0]=="gsimpl")
                    variables.query_type = GLOBAL_MORSE_SIMPLIFICATION;

                if(variables.query_type == WINDVT || variables.query_type == WINDTT ||
                        variables.query_type == POINT || variables.query_type == BOX)
                    variables.query_path = tok[1];
                else
                    variables.persistence = atof(tok[1].c_str());
            }
            i++;
        }        
        else if(strcmp(tag, "-output") == 0)
        {
            variables.app_debug = OUTPUT;
        }
        else if(strcmp(tag, "-vtime") == 0)
        {
            variables.app_debug = TIME_VERBOSE;
        }
    }

    return 0;
}

void print_usage()
{
    cerr<<"[ERROR] Wrong Usage. Run ./terrain_trees for detailed instructions."<<endl;
}

void print_help()
{
    //annoying stuff to get the dimension of the output shell (!!! not sure it works on Mac,
    //everything based on the command tput. If it doesn't work the dimension si setted to 80 by default)
    FILE* fp;
    char path[1035];

    int cols;
    fp = popen("tput cols", "r");
    if(fp != NULL){
        fgets(path, sizeof(path)-1, fp);
        cols = atoi(path);
    }
    else{
        cols = 80;
    }

    printf(BOLD "\n  NAME:\n\n" RESET);
    printf("\tTerrain Trees library - Multivariate Morse Analysis\n\n" RESET);

    printf(BOLD "  USAGE: \n\n" RESET);
    printf(BOLD "    ./terrain_trees {<-v [kv] -t [kt] -c [crit] -d [div] | -f [tree_file]>\n" RESET);
    printf(BOLD "                       -q [op-file | app] -s -noR} | {-g [query-ratio-quantity-type]}\n" RESET);
    printf(BOLD "                       -i [mesh_file]\n" RESET);

    printf(BOLD "    -v [kv]\n" RESET);
    print_paragraph("kv is the vertices threshold per leaf. This parameter is needed by PR-T and PM-T trees.", cols);
    printf(BOLD "    -t [kt]\n" RESET);
    print_paragraph("kt is the triangles threshold per leaf. This parameter is needed by PMR-T and PM-T trees.", cols);
    printf(BOLD "    -c [crit]\n" RESET);
    print_paragraph("crit is the criterion type of the index. This can be PR-T tree (pr), PMR-T tree (pmr) and  PM-T tree (pm).", cols);
    printf(BOLD "    -d [div]\n" RESET);
    print_paragraph("div is the division type of the index. This can be quadtree (quad) or kD-tree (kd).", cols);

    print_paragraph("NOTA: these arguments must be used in conjunction in order to create an index. "
                    "This operation generate as output a file containing the triangles index.", cols);

    printf(BOLD "    -f [tree_file]\n" RESET);
    print_paragraph("reads an spatial index from an input file", cols);
    print_paragraph("tree_file contains a terrain tree. This file has a fixed syntax of the name "
                    "that allows to recover the informations needed to initialize the index (i.e., kv, kt, division and criterion types)", cols);

    print_paragraph("NOTA: you can use -f argument [OR] {-v / -t / -c / -d} accordingly to the chosen criterion.", cols);

    printf(BOLD "    -q [app]\n" RESET);
    print_paragraph("'app' can be: batch - concurv - mccurv - gccurv - slopes - crit - morse - simpl - gsimpl - multiv - filter "
                    "'batch' extracts VT and TT relations on the whole mesh, "                    
                    "'concurv' extracts the Concentrated Curvature, "
                    "'mccurv' extracts the Mean CCurvature, "
                    "'gccurv' extracts the Gaussian CCurvature, "
                    "'slopes' extracts the slope values for the triangles and edges of the terrain, "
                    "'rough' extracts the roughness values for the vertices of the terrain, "
                    "'multifield' extracts the multifield measure for the vertices of the terrain (Note that the input data should have more than one field), "
                    "'crit' extracts the critical points of the terrain, "                    
                    "'morse' computes the Forman gradient vector, extracts 2/1 descending/ascending manifolds and the Morse indicence graph, "
                    "'simpl' computes the Forman gradient vector and simplify it using a local Morse incidence graph, [Developing]"
                    "'gsimpl' computes the Forman gradient vector and simplify it using a global Morse incidence graph, "
                    "'multiv' computes the Multivariate Morse analysis by generating the Multivariate gradient and extracting the Critical Clusters, "
                    "'filter' reads a points cloud generates a PR tree, extracts the multifield of each 2D point and finally outputs both"
                    " the multifield file and a points cloud file compatible with SpatialHadoop.", cols);
    print_paragraph("NOTA: simpl and gsimpl can have an optional persistence parameter, with range from x to y, following the format "
                    "'(g)simpl'-'persistence'.If not set, it is used the default persistence value 0.65.", cols);
    printf(BOLD "    -g [query-ratio-quantity-type]\n" RESET);
    print_paragraph("generates a given number of input data for a specific query", cols);
    print_paragraph("query can be: point - box. "
                    "'ratio' is a number between 0 and 1, and and represents the percentage of the maximum side of the domain to pick. "
                    "'quantity' is a positive number that indicate the number of inputs to generate. "
                    "type can be: near - rand. 'near' stands for a randomly generated point that is near the mesh, "
                    "while 'rand' stands for a randomly generated point that is inside the domain.", cols);
    print_paragraph("If 'query' is equal to point 'ratio' must be equal to 0, otherwise 'ratio' must be greater than 0.", cols);

    printf(BOLD "    -s\n" RESET);
    print_paragraph("computes the statistics of a tree.", cols);
    printf(BOLD "    -noR\n" RESET);
    print_paragraph("disable the procedures that exploits the spatial coherence of the index and the mesh. "
                    "Notice that only spatial queries can be executed on this type of index.", cols);
    printf(BOLD "    - i [mesh_file]\n" RESET);
    print_paragraph("reads the mesh_file containing the triangle mesh. The mesh can be in .tri, .off or .soup formats. "
                    "The .soup format contains a triangulated terrain in which the vertices coordinates are represented within each triangle in their star.", cols);

    printf(BOLD "    -output\n" RESET);
    print_paragraph("if the application has this feature, save to file the executed analysis.", cols);
    printf(BOLD "    -vtime\n" RESET);
    print_paragraph("if the application has this feature, outputs detailed execution timing.", cols);

    printf(BOLD "  EXAMPLE[1]: \n" RESET);
    printf("          ./terrain_trees -v 20 -c pr -d quad -s -i mesh.off\n");
    print_paragraph("First it reads the mesh [mesh.off]. Then, builds a PR-T tree with kv=20 and quadtree subdivision, "
                    "on which it is exploited the spatial coherence of the mesh and index. "
                    "Finally, it computes the index statistics (-s).", cols);

    printf(BOLD "  EXAMPLE[2]: \n" RESET);
    printf("          ./terrain_trees -f tree_file -q wvt-boxfile -i mesh.tri\n");
    print_paragraph("First it reads the mesh [mesh.tri]. Then, it reads the spatial index from .tree file (gathering the tree "
                    "parameters direcly from the file name) and exploits the spatial coherence of the mesh and index. "
                    "Finally, it executes the windowed VT queries, using the boxes into 'boxfile'.", cols);

    printf(BOLD "  EXAMPLE[3]: \n" RESET);
    printf("          ./terrain_trees -v 20 -c pr -d quad -q crit -i mesh.soup -output\n");
    print_paragraph("First it reads the soup [mesh.soup].  Then, it builds a PR-T tree index with kv=20 and with quadtree subdivision. "
                    "As last generation step, it exploits the spatial coherence of the mesh and index. "
                    "Then, it extracts the indexed mesh representation of the soup and it saves that in a .off file. "
                    "Finally, it computes the critical points, outputting them in .vtk files for visualization purposes (-output parameter).", cols);


    printf(BOLD "  EXAMPLE[4]: \n" RESET);
    printf("          ./terrain_trees -v 20 -c pr -d quad -q morse -i mesh.soup -output\n");
    print_paragraph("First it reads the soup [mesh.soup].  Then, it builds a PR-T tree index with kv=20 and with quadtree subdivision. "
                    "As last generation step, it exploits the spatial coherence of the mesh and index. "
                    "Then, it extracts the indexed mesh representation of the soup and it saves that in a .off file. "
                    "Finally, it computes the Forman gradient and extracts the morphological features, "
                    "outputting them in .vtk files for visualization purposes (-output parameter).", cols);

    printf(BOLD "  IMPLEMENTATION:\n" RESET);
    printf("          Author: Riccardo Fellegara, Yunting Song\n");
    printf("          Group: UMD GeoVisualization Group\n");
    printf("          Man-page Last Update: September 2021\n\n");

    printf(BOLD "  DESCRIPTION: \n" RESET);
    print_paragraph("We propose a new in-core family of spatial indexes for the analysis "
                    "of Triangulated Irregular Networks (TINs). We call such indexes "
                    "Terrain trees. A Terrain tree combines a minimal encoding of the "
                    " connectivity of the underlying triangle mesh with a hierarchical spatial "
                    "index implicitly encoding the other topological relations among "
                    "the mesh elements. Topological relations are extracted locally within "
                    "each leaf block of the hierarchy at runtime, based on specific application needs. "
                    "We introduce a new tool for the multivariate analysis of the surface by combining different scalar "
                    "fields (i.e., the elevation and the curvature values). By computing a "
                    "combinatorial discrete vector field we are able to study the mutual "
                    "relationships between the different fields providing new insight on "
                    "the terrain morphology. Moreover, we have developed other state-of-the-art estimators, "
                    "such as slope estimation, curvature computation, "
                    "and the extraction of the terrain critical points on the Terrain trees.", cols);
}

void print_paragraph(string stringa, int cols){
    if((int)stringa.size() < cols-20){
        printf("          %s\n\n",stringa.c_str());
    }
    else{
        float dim = (float)(stringa.size()/((float)cols-20));
        int dim_int = dim;
        if(dim > dim_int) dim_int++;
        for(int i=0; i<dim_int; i++){
            if(stringa.at((cols-20)*i) == ' ')
                printf("         %s\n", stringa.substr( (cols-20)*i, cols-20).c_str());
            else printf("          %s\n", stringa.substr( (cols-20)*i, cols-20).c_str());
        }
        printf("\n");
    }
}

}

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

#ifndef CLI_PARAMETERS
#define CLI_PARAMETERS

#include <string>
#include "basic_types/basic_wrappers.h"
using namespace std;

//this define the operation type that can executed during the forman gradient computation and during the feature extraction
// TIME_VERBOSE - gets the extra verbose timings of the intermidiate operations without saving nothing into the global structure
// OUTPUT - saves globally the results to allow a visual debugging
// NOTHING - is the default behaviour
typedef enum { TIME_VERBOSE, OUTPUT, NOTHING } OpType;

#define DEFAULT_X_PER_LEAF -1
#define DEFAULT "null"
#define CACHESIZE 200
enum QueryType { POINT, BOX, WINDVT, WINDTT, BATCH,
                 CONCENTRATED_CURVATURE, MEAN_CCURVATURE, GAUSS_CCURVATURE,
                 MORSE_ANALYSIS, LOCAL_MORSE_SIMPLIFICATION, GLOBAL_MORSE_SIMPLIFICATION, MULTIVARIATE_MORSE_ANALYSIS,
                 FILTER,
                 SLOPES/*, ESLOPE, TSLOPE*/, CRITICAL_POINTS,
                 CUSTOM,
                 NULL_QUERY };
enum SpatialDecType {KD = 2, QUAD = 4, UNSET = -1};
#define BOLD  "\033[1m\033[33m" //for dark background shell
//#define BOLD "\033[1m\033[31m"  //for white background shell
#define RESET   "\033[0m"

/// GLOBAL VARIABLES
struct cli_parameters
{
    string mesh_path, query_path, exe_name, tree_path;
    SpatialDecType division_type;
    string crit_type;
    bool is_index, is_getInput, isTreeFile, reindex, is_multified;
    int v_per_leaf;
    int t_per_leaf;
    bool debug_mode;
    bool debug_prints_mode;

    //for simplification
    bool QEM_based;
    bool contract_all_edges;
    double maximum_limit; 
    int num_of_threads;
    bool evaluation_mode;

    int num_input_entries;
    coord_type ratio;
    QueryType query_type;
    string input_gen_type;

    //global variables (used only during a printing phase) which store the association newSimplIndex -> oldSimplIndex
    ivect original_vertex_indices;
    ivect original_triangle_indices;
    dvect original_vertex_fields;
    bool rever_to_original;
    OpType app_debug;
    int cache_size;
    coord_type persistence;
    string weights_file;

    // for gradient-aware simplification validation
    ivect new_vertex_indices;

    cli_parameters()
    {
        division_type = UNSET;
        crit_type = DEFAULT;
        v_per_leaf = DEFAULT_X_PER_LEAF;
        t_per_leaf = DEFAULT_X_PER_LEAF;
        query_type = NULL_QUERY;
        is_index = false;
        is_getInput = false;
        isTreeFile = false;
        reindex = true; /// the default behavior is to exploit the spatial coherence

        num_input_entries = 0;
        input_gen_type = DEFAULT;

        weights_file = DEFAULT;

        app_debug = NOTHING;
        cache_size = CACHESIZE;
        rever_to_original = false;
        persistence = 0.65;

        is_multified = false;

        debug_mode = true;
        debug_prints_mode = false;
        QEM_based = false;
        contract_all_edges = false;
        evaluation_mode = true;

        num_of_threads = 1;
    }

     inline bool exec_weighted_simplification() { return weights_file != DEFAULT; }
};
///

#endif // CLI_PARAMETERS


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

#ifndef MAIN_UTILITY_FUNCTIONS_H
#define MAIN_UTILITY_FUNCTIONS_H

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <string>
#include <sstream>


#include "geometry/geometry_wrapper.h"
#include "io/reader.h"
#include "io/writer_morse.h"
#include "queries/spatial_queries.h"
#include "queries/topological_queries.h"
#include "statistics/statistics.h"
#include "terrain_trees/spatial_subdivision.h"
#include "terrain_trees/prt_tree.h"
#include "terrain_trees/pmt_tree.h"
#include "terrain_trees/pmrt_tree.h"
#include "terrain_trees/reindexer.h"
#include "utilities/input_generator.h"
#include "utilities/string_management.h"
#include "utilities/timer.h"
#include "curvature/concentrated_curvature.h"
#include "curvature/c_curvature.h"
#include "terrain_features/slope_extractor.h"
#include "terrain_features/critical_points_extractor.h"
#include "edge_contraction/contraction_simplifier.h"
#include "edge_contraction/gradient_aware_contraction.h"
#include "edge_contraction/smooth.h"
#include "edge_contraction/distance_calculator.h"

#include "morse/forman_gradient.h"
#include "morse/forman_gradient_computation.h"
#include "morse/forman_gradient_features_extractor.h"
#include "morse/forman_gradient_simplifier.h"
#include "morse/boundary_matrix.h"
#include "utilities/cli_parameters.h"
#include "utilities/usage.h"

using namespace std;
using namespace string_management;

namespace utility_functions
{
    extern string SpatialDecType2string(SpatialDecType sdt);
    extern int read_arguments(int argc, char **argv, cli_parameters &variables);
    extern void print_usage();
    extern void setParameters(cli_parameters &variables);
    extern bool checkParameters(cli_parameters &variables);
    extern void print_help();
    extern void print_paragraph(string stringa, int cols);
} // namespace utility_functions

#endif // MAIN_UTILITY_FUNCTIONS_H

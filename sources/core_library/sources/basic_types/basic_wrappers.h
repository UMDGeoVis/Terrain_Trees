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

#ifndef BASIC_STRUCTURE
#define BASIC_STRUCTURE

#include <vector>
#include <set>
#include <queue>
#include <map>
#include <unordered_map>
#include <omp.h>

#ifdef LONG_TYPES
typedef size_t utype;
typedef ssize_t itype;
#else
typedef unsigned utype;
typedef int itype;
#endif

typedef double coord_type;

typedef std::queue<itype> iqueue;


typedef std::vector<itype> ivect;
typedef ivect::iterator ivect_iter;
typedef ivect::const_iterator ivect_citer;

typedef std::vector<utype> uvect;
typedef uvect::iterator uvect_iter;
typedef uvect::const_iterator uvect_citer;

typedef std::vector<coord_type> dvect;
typedef dvect::iterator dvect_iter;
typedef dvect::const_iterator dvect_citer;

typedef std::set<itype> iset;
typedef iset::iterator iset_iter;
typedef iset::const_iterator iset_citer;

typedef std::set<ivect> ivect_set;

typedef std::set<utype> uset;
typedef uset::iterator uset_iter;
typedef uset::const_iterator uset_citer;

typedef std::set<coord_type> dset;
typedef dset::iterator dset_iter;
typedef dset::const_iterator dset_citer;

typedef std::map<itype,dset> vertex_multifield;

typedef ivect VT;
typedef iset leaf_ids;
typedef iset VV;
typedef ivect VV_vec;

typedef std::vector<VT> leaf_VT;
typedef std::vector<leaf_ids> lists_leafs;
typedef std::vector<VV> leaf_VV;
typedef std::vector<VV_vec> leaf_VV_vec;

typedef ivect leaf_VTstar;

typedef std::pair<itype,itype> ET;
typedef std::map<ivect,ET> leaf_ET;
typedef std::map<ivect,itype> leaf_ETstar;

typedef std::map<itype,ivect> wTT;
typedef std::map<itype,VT> wVT;

typedef dvect FG;//store the gradient of one function.
typedef std::vector<dvect>vect_FG;



#endif // BASIC_STRUCTURE


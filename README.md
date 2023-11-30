# Terrain Trees library #

Terrain trees are a new in-core family of spatial indexes for the representation 
and analysis of Triangulated Irregular Networks (TINs).
Terrain trees combine a minimal encoding of the connectivity of the
underlying triangle mesh with a hierarchical spatial index, implicitly
representing the other topological relations among vertices, edges
and vertices. Topological relations are extracted locally within each
leaf block of the hierarchal index at runtime, based on specific application 
needs. We have developed a tool based on Terrain trees for
terrain analysis, which includes state-of-the-art estimators for slope
and curvature, and for the extraction of critical points. 
By working on TINs generated from big LiDAR (Light, Detection and Ranging) 
data sets, we demonstrate the effectiveness and scalability of the 
Terrain trees against a state-of-the-art compact data structures.

### Reference Paper ###

Riccardo Fellegara, Yunting Song, Federico Iuricich, and Leila De Floriani.
*Terrain trees: a framework for representing, analyzing and visualizing triangulated terrains*. GeoInformatica 27.3 (2023): 525-564. [doi](https://link.springer.com/article/10.1007/s10707-022-00472-3)

Riccardo Fellegara, Federico Iuricich, and Leila De Floriani. 
*Efficient representation and analysis of triangulated terrains*.
In Proceedings of SIGSPATIAL’17, Los Angeles Area, CA, USA, November 7–10, 2017, 4 pages.
[doi](https://dl.acm.org/citation.cfm?id=3139958.3140050)

### Features ###

+ Three spatial indexes based on
    * point threshold (PR-T tree)
    * triangle threshold (PMR-T tree)
    * point and triangle thresholds (PM-T tree)
+ Two spatial decomposition
    * quadtree
    * kD-tree
+ Spatial queries
    * point location
    * box query
    * incremental nearest neighbor ([reference paper](http://link.springer.com/chapter/10.1007%2F3-540-60159-7_6))
+ Terrain Features
    * Triangle/Edges slope computation
    * Critical Points extraction
    * Roughness computation ([reference paper](https://doi.org/10.1029/2000JE001429))
+ Curvature computation ([reference1](http://dl.acm.org/citation.cfm?id=1463498) and [reference2](http://www.umiacs.umd.edu/~deflo/papers/2010grapp/2010grapp.pdf))
    * Concentrated curvature
    * Mean and Gaussian CCurvature 
+ Soup to indexed mesh conversion
+ Points cloud indexing
    * multifield extraction ([reference paper](https://onlinelibrary.wiley.com/doi/10.1111/j.1467-8659.2011.01959.x))
+ Topology Data Analytics ([reference1](https://doi.org/10.1111/cgf.12123) and [reference2](https://doi.org/10.1145/2666310.2666412))
    * Gradient computation
    * Critical points extraction
    * Data segmentation / Critical net extraction
    * Topological simplification (PR-T tree only)

### How to compile ###

The library requires the [boost library](http://www.boost.org/) (for dynamic_bitset class), [Eigen3 library](http://eigen.tuxfamily.org/index.php?title=Main_Page) (for computing  eigenvalues) and [cmake](https://cmake.org/) installed in your system.

Once in the root of the repository type from the command line
```
cmake CMakeList.txt
```
and once configured
```
make
```
This command generates a portable library file, located into *lib* folder, as well as an executable into the `bin` folder.

The compilation process has been test on linux and mac systems.

### Execute unit tests ###

To compile the unit-tests run from command line the following command
```
make tests
```
Once compiled, it is possible to test the main functionalities running a script located into the `bin` folder.

From command line executing 
```
sh run_tests.sh
```
checks the functionalities of the main implemented features.
The output files are saved into the `data` folder (where the input datasets are located).

It is possible to clean the output files running, from command line (from `data` folder), the following command
```
sh clean_up.sh
```

### Use the main library ###

In the `bin` folder there is the main executable file named `terrain_trees` that contains the whole library. For a complete list of the command line options refer the [wiki](https://github.com/FellegaraR/Terrain_Trees/wiki/Command-line-parameters) page.

### Supported Input Files ###

The Terrain trees framework supports three input formats for the triangulated irregular network:
+ off
+ tri
+ soup

For a detailed description of the input formats refer the [wiki](https://github.com/FellegaraR/Terrain_Trees/wiki/Supported-Input-Formats) page.

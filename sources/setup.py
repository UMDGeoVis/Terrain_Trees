from setuptools import setup

from Cython.Build import cythonize

from distutils.core import setup, Extension

import glob

cpp_files = glob.glob('**/*.cpp', recursive=True)

cpp_files.remove('targets/main.cpp')

tests_files = glob.glob('tests/*.cpp')
for tests_name in tests_files:
    cpp_files.remove(tests_name)

if 'py_terrain_trees.cpp' in cpp_files:
    cpp_files.remove('py_terrain_trees.cpp')

sources_print = ["py_terrain_trees.pyx"] + cpp_files

ext = [Extension("PythonMain",
    sources= sources_print,
    language="c++",
    include_dirs=['core_library/sources/', 'core_library/sources/terrain_trees',
                    'core_library/sources/utilities','core_library/sources/basic_types',
                    'core_library/sources/curvature','core_library/sources/geometry',
                    # 'core_library/sources/gradient',
                    'core_library/sources/io',
                    'core_library/sources/queries','core_library/sources/roughness',
                    'core_library/sources/statistics','core_library/sources/terrain_features','/usr/include/eigen3'])]
setup(ext_modules=cythonize(ext))

# use command 
# python setup.py build_ext --inplace
# to compile and create wrapper

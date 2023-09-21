set(INCLUDE_DIRS
		basic_types
		curvature
		geometry
		io		
		queries
		statistics
		terrain_trees
		utilities
		morse
		roughness
		gradient
		apps
	)
	
INCLUDE_DIRECTORIES( ${INCLUDE_DIRS} )

file(GLOB BASIC_TYPES_SRC
		basic_types/*.h
		basic_types/*.cpp
   )

file(GLOB CURVATURE_SRC
		curvature/*.h
        curvature/*.cpp
	)

file(GLOB ROUGHNESS_SRC
		roughness/*.h
        roughness/*.cpp
	)
file(GLOB GRADIENT_SRC
		gradient/*.h
        gradient/*.cpp
	)

file(GLOB GEOMETRY_SRC
		geometry/*.h
		geometry/*.cpp
	)
	
file(GLOB IO_SRC
		io/*.h
		io/*.cpp
	)

	
file(GLOB QUERIES_SRC		
		queries/*.h
		queries/*.cpp
	)
    
file(GLOB STATISTICS_SRC
		statistics/*.h
		statistics/*.cpp
	)

file(GLOB TERRAIN_FEATURES_SRC
        terrain_features/*.h
        terrain_features/*.cpp
    )
	
file(GLOB TRIANGLE_TREES_SRC
		terrain_trees/*.h
		terrain_trees/*.cpp
	)
    

file(GLOB UTILITIES_SRC 
		utilities/*.h
        utilities/*.cpp
	)
	
file(GLOB MORSE_SRC 
		morse/*.h
        morse/*.cpp
	)

file(GLOB GEOG_APPS_SRC
		apps/*.h
		apps/*.cpp
)

set(TERRAIN_TREES_CORE_SRC
        ${BASIC_TYPES_SRC}
		${GEOMETRY_SRC}
		${IO_SRC}
		${STATISTICS_SRC}
		${TRIANGLE_TREES_SRC}
		${UTILITIES_SRC}
	)

set(TERRAIN_TREES_APPS_SRC
                ${CURVATURE_SRC}
                ${ROUGHNESS_SRC}
				${GRADIENT_SRC}
                ${QUERIES_SRC}
                ${TERRAIN_FEATURES_SRC}
				${MORSE_SRC}
				${GEOG_APPS_SRC}
        )

add_library(_terrain_trees_core ${TERRAIN_TREES_CORE_SRC})
add_library(_terrain_trees_apps ${TERRAIN_TREES_APPS_SRC})

target_link_libraries(_terrain_trees_apps _terrain_trees_core)




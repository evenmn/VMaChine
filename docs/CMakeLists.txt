find_package(Doxygen REQUIRED)
 
# Find all the public headers
# get_target_property(VMACHINE_PUBLIC_HEADER_DIR SourceFiles INTERFACE_INCLUDE_DIRECTORIES)
# file(GLOB_RECURSE VMACHINE_PUBLIC_HEADERS VMACHINE_EIGEN_PUBLIC_HEADERS ${VMACHINE_PUBLIC_HEADER_DIR}/*.h)
 
#This will be the main output of our command
set(DOXYGEN_INDEX_FILE ${CMAKE_CURRENT_SOURCE_DIR}/html/index.html)
 
# add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
#                    DEPENDS ${VMACHINE_PUBLIC_HEADERS}
#                    COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile
#                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#                    MAIN_DEPENDENCY Doxyfile
#                    COMMENT "Generating docs")
 
add_custom_target(Doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})

cmake_minimum_required (VERSION 3.10)
project(SPADES-Docs NONE)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Setup paths
set(SPADES_PROJECT_DIR ${CMAKE_SOURCE_DIR}/..)
set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/sphinx)
set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/sphinx)
set(DOXYGEN_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/doxygen)
set(DOXYGEN_INPUT_DIR ${SPADES_PROJECT_DIR}/Source)
set(DOXYGEN_OUTPUT_DIR ${SPHINX_SOURCE}/doxygen)
set(DOXYGEN_HTML_OUTPUT_DIR ${DOXYGEN_OUTPUT_DIR}/html)
file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR})

# Doxygen target
find_package(Doxygen REQUIRED)

file(GLOB_RECURSE SPADES_PUBLIC_SOURCE ${DOXYGEN_INPUT_DIR}/*.H ${DOXYGEN_INPUT_DIR}/*.cpp)
set(DOXYGEN_INDEX_FILE ${DOXYGEN_HTML_OUTPUT_DIR}/index.html)
set(DOXYGEN_TAG_FILE ${DOXYGEN_HTML_OUTPUT_DIR}/tagfile.xml)
set(DOXYGEN_AWESOME_STYLE_FILE ${SPADES_PROJECT_DIR}/Submodules/doxygen-awesome-css/doxygen-awesome.css)
set(DOXYFILE_IN ${DOXYGEN_SOURCE_DIR}/Doxyfile.in)
set(DOXYFILE_OUT ${DOXYGEN_OUTPUT_DIR}/Doxyfile)
configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
                   DEPENDS ${SPADES_PUBLIC_SOURCE}
                   DEPENDS ${DOXYGEN_SOURCE_DIR}/main_page.md
                   COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
                   MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
                   COMMENT "Generating documentation with Doxygen")

add_custom_target(doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})

# Doxysphinx target
find_package(Doxysphinx REQUIRED)

add_custom_target(doxysphinx ALL
  COMMAND ${DOXYSPHINX_EXECUTABLE} build
  ${SPHINX_SOURCE} ${SPHINX_BUILD}/html ${DOXYGEN_HTML_OUTPUT_DIR}
  COMMENT "Converting doxygen documentation to sphinx with doxysphinx")
add_dependencies(doxysphinx doxygen)


# Sphinx target
find_package(Sphinx REQUIRED)

add_custom_target(sphinx ALL
  COMMAND ${SPHINX_EXECUTABLE} -M html
  ${SPHINX_SOURCE} ${SPHINX_BUILD}
  -c ${SPHINX_SOURCE} -W --keep-going -n
  COMMENT "Generating documentation with Sphinx")
add_dependencies(sphinx doxysphinx)

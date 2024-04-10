target_sources(${spades_lib_name}
  PRIVATE
    CellSortedParticleContainer.H
    CellSortedParticleContainer.cpp
    Constants.H
    Factory.H
    IC.H
    IC.cpp
    SPADES.H
    SPADES.cpp
    ParticleData.H
    ParticleOps.H
    Utilities.H
    Utilities.cpp
)

target_sources(${spades_exe_name} PRIVATE main.cpp)

target_include_directories(${spades_lib_name} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>)

include(AMReXBuildInfo)
generate_buildinfo(${spades_lib_name} ${CMAKE_SOURCE_DIR})

configure_file("${CMAKE_SOURCE_DIR}/CMake/SpadesVersion.H.in"
  "${CMAKE_CURRENT_BINARY_DIR}/SpadesVersion.H" @ONLY)

target_link_libraries_system(${spades_lib_name} PUBLIC AMReX::amrex)
target_link_libraries(${spades_exe_name} PRIVATE ${spades_lib_name})

install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/SpadesVersion.H
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

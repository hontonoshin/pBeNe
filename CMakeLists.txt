cmake_minimum_required(VERSION 3.16)
project(pBN_multichannel VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#──────────────────────────────────────────────────────────────────────────────
#  CC=/usr/bin/gcc CXX=/usr/bin/g++ cmake ..
#──────────────────────────────────────────────────────────────────────────────
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

#──────────────────────────────────────────────────────────────────────────────
# Sources & headers
#──────────────────────────────────────────────────────────────────────────────
set(SOURCES
    src/ActionInitialization.cc
    src/B11ProtonInelasticXS.cc
    src/ChannelCrossSectionData.cc
    src/ChannelFinalStateModel.cc
    src/DetectorConstruction.cc
    src/EventAction.cc
    src/main.cc
    src/MultiChannelSampler.cc
    src/PhysicsListWithTALYS.cc
    src/PrimaryGeneratorAction.cc
    src/RunAction.cc
    src/SteppingAction.cc
    src/TALYSFinalStateModel.cc
    src/TALYSNeutronSampler.cc
)

add_executable(pBN ${SOURCES})

target_include_directories(pBN PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${Geant4_INCLUDE_DIRS}
)

target_link_libraries(pBN PRIVATE ${Geant4_LIBRARIES})

#──────────────────────────────────────────────────────────────────────────────
# Copy macro files and data symlink to build directory
#──────────────────────────────────────────────────────────────────────────────
foreach(_mac run.mac vis.mac)
    if(EXISTS ${CMAKE_SOURCE_DIR}/${_mac})
        configure_file(${CMAKE_SOURCE_DIR}/${_mac}
                       ${CMAKE_CURRENT_BINARY_DIR}/${_mac} COPYONLY)
    endif()
endforeach()

# Create a symlink to data/ so the binary finds it at data/
if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/data)
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E create_symlink
            ${CMAKE_SOURCE_DIR}/data
            ${CMAKE_CURRENT_BINARY_DIR}/data
    )
endif()

#──────────────────────────────────────────────────────────────────────────────
# Install
#──────────────────────────────────────────────────────────────────────────────
install(TARGETS pBN DESTINATION bin)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/data DESTINATION share/pBN)

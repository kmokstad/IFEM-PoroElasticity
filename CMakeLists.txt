project(PoroElasticity)

cmake_minimum_required(VERSION 2.6)

# Add local modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${PROJECT_SOURCE_DIR}/../../cmake/Modules)

if(NOT IFEM_CONFIGURED)
  find_package(IFEM REQUIRED)
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IFEM_CXX_FLAGS}")

# Generate regtest script with correct paths
configure_file(${PROJECT_SOURCE_DIR}/../Common/scripts/regtest.sh.in
               regtest.sh)

include_directories(${IFEM_INCLUDES} ../Common)

set(PoroElasticity_SRCS main_PoroElasticity.C
                        PoroElasticity.C
                        PoroMaterial.C)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

add_executable(PoroElasticity ${PoroElasticity_SRCS})
target_link_libraries(PoroElasticity ${IFEM_LIBRARIES})

# Installation
install(TARGETS PoroElasticity DESTINATION bin COMPONENT bin)

enable_testing()

# Regression tests
if(IFEM_USE_PARALLEL_PETSC)
  # Add parallel tests here
else()
  file(GLOB POROEL_TESTFILES RELATIVE ${PROJECT_SOURCE_DIR}/Test
                                      ${PROJECT_SOURCE_DIR}/Test/*.reg)
  foreach(TESTFILE ${POROEL_TESTFILES})
    ifem_add_test(${TESTFILE} PoroElasticity)
  endforeach()
endif()
list(APPEND TEST_APPS PoroElasticity)
if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
  set(UNIT_TEST_NUMBER ${UNIT_TEST_NUMBER} PARENT_SCOPE)
else()
  add_check_target()
endif()

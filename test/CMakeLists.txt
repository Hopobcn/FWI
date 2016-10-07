add_executable(fwi-tests
    fwi_tests.c
    fwi_common_tests.c
    fwi_propagator_tests.c
    fwi_kernel_tests.c
)

target_include_directories(fwi-tests PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(fwi-tests
    fwi-core
    Unity
    m
)

add_test(fwi_test_1 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/fwi-tests)

# (option1: use 'make ctest') LAUNCH CTEST with build dependences
add_custom_target(ctest
    COMMAND ${CMAKE_CTEST_COMMAND}
    DEPENDS fwi-tests
)

# (option2: use 'make utest') LAUNCH Unity binary with build deps
add_custom_target(utest
    COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/fwi-tests -v
    DEPENDS fwi-tests
)

add_executable(fwi-tests
    fwi-tests.c
)

target_link_libraries(fwi-tests
    fwi-core
    Unity
    m
)

add_test(fwi_test_1 fwi-tests)

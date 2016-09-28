add_library(Unity STATIC 
    Unity/src/unity.c
    Unity/extras/fixture/src/unity_fixture.c
)

target_include_directories(Unity PUBLIC 
    Unity/src
    Unity/extras/fixture/src
    ${PROJECT_SOURCE_DIR}/include/test
)

target_compile_definitions(Unity PUBLIC
    "-DUNITY_INCLUDE_CONFIG_H"
)

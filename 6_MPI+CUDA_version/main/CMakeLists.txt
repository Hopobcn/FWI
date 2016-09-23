add_executable(fwi
    fwi_main.c
)

target_link_libraries(fwi
    fwi-core
    m
)

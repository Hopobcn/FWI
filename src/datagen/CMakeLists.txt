add_executable(fwi-data-generator
    fwi_generatemodel.cpp
)

target_link_libraries(fwi-data-generator
    fwi-core
    m
)

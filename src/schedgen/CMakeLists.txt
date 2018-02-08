add_executable(fwi-sched-generator
    fwi_schedule.c
)

target_include_directories(fwi-sched-generator PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(fwi-sched-generator
    fwi-core
    m
)

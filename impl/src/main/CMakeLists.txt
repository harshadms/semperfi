# Copyright (C) 2012-2020  (see AUTHORS file for a list of contributors)
#
# GNSS-SDR is a software-defined Global Navigation Satellite Systems receiver
#
# This file is part of GNSS-SDR.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

if(USE_CMAKE_TARGET_SOURCES)
    add_executable(semperfi)
    target_sources(semperfi PRIVATE main.cc)
else()
    add_executable(semperfi main.cc)
endif()


if(${FILESYSTEM_FOUND})
    target_compile_definitions(semperfi PRIVATE -DHAS_STD_FILESYSTEM=1)
    if(${find_experimental})
        target_compile_definitions(semperfi PRIVATE -DHAS_STD_FILESYSTEM_EXPERIMENTAL=1)
    endif()
    target_link_libraries(semperfi PRIVATE std::filesystem)
else()
    target_link_libraries(semperfi PRIVATE Boost::filesystem Boost::system)
endif()

target_link_libraries(semperfi
    PRIVATE
        algorithms_libs
        core_receiver
        Boost::headers
        Boost::thread
        Gflags::gflags
        Glog::glog
        Threads::Threads
)

target_compile_definitions(semperfi PRIVATE -DGNSS_SDR_VERSION="${VERSION}")

if(ENABLE_CUDA)
    if(NOT CMAKE_VERSION VERSION_GREATER 3.11)
        target_link_libraries(semperfi
            PRIVATE
                ${CUDA_LIBRARIES}
        )
        target_include_directories(semperfi
            PRIVATE
                ${CUDA_INCLUDE_DIRS}
        )
    endif()
    target_compile_definitions(semperfi PRIVATE -DCUDA_GPU_ACCEL=1)
endif()

if(ENABLE_GPERFTOOLS)
    if(GPERFTOOLS_FOUND)
        target_link_libraries(semperfi
            PRIVATE
                Gperftools::profiler
                Gperftools::tcmalloc
        )
    endif()
endif()

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        target_link_libraries(semperfi PRIVATE "-lc++")
    endif()
endif()

if(ENABLE_STRIP)
    set_target_properties(semperfi PROPERTIES LINK_FLAGS "-s")
endif()

if(ENABLE_CLANG_TIDY)
    if(CLANG_TIDY_EXE)
        set_target_properties(semperfi
            PROPERTIES
                CXX_CLANG_TIDY "${DO_CLANG_TIDY}"
        )
    endif()
endif()

add_custom_command(TARGET semperfi
    POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:semperfi>
        ${LOCAL_INSTALL_BASE_DIR}/install/$<TARGET_FILE_NAME:semperfi>
)

install(TARGETS semperfi
    RUNTIME DESTINATION bin
    COMPONENT "semperfi"
)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/conf
    DESTINATION share/semperfi
    FILES_MATCHING PATTERN "*.conf"
)

install(FILES ${CMAKE_SOURCE_DIR}/conf/gnss-sdr.conf
    DESTINATION share/gnss-sdr/conf
    RENAME default.conf
)

if(NOT VOLKGNSSSDR_FOUND)
    install(PROGRAMS ${CMAKE_BINARY_DIR}/volk_gnsssdr_module/build/apps/volk_gnsssdr_profile
        DESTINATION bin
        COMPONENT "volk_gnsssdr"
    )

    install(PROGRAMS ${CMAKE_BINARY_DIR}/volk_gnsssdr_module/build/apps/volk_gnsssdr-config-info
        DESTINATION bin
        COMPONENT "volk_gnsssdr"
    )
endif()

find_program(GZIP
    gzip
    /bin
    /usr/bin
    /usr/local/bin
    /opt/local/bin
    /sbin
)

if(NOT GZIP_NOTFOUND)
    execute_process(COMMAND gzip -9 -c ${CMAKE_SOURCE_DIR}/docs/manpage/semperf-manpage
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        OUTPUT_FILE "${CMAKE_BINARY_DIR}/gnss-sdr.1.gz"
    )

    install(FILES ${CMAKE_BINARY_DIR}/gnss-sdr.1.gz
        DESTINATION share/man/man1
    )

    execute_process(COMMAND gzip -9 -c ${CMAKE_SOURCE_DIR}/docs/changelog.md
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        OUTPUT_FILE "${CMAKE_BINARY_DIR}/changelog.gz"
    )

    install(FILES ${CMAKE_BINARY_DIR}/changelog.gz
        DESTINATION share/doc/gnss-sdr
    )

    if(NOT VOLKGNSSSDR_FOUND)
        execute_process(COMMAND gzip -9 -c ${CMAKE_SOURCE_DIR}/src/algorithms/libs/volk_gnsssdr_module/volk_gnsssdr/cmake/Packaging/volk_gnsssdr_profile-manpage
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            OUTPUT_FILE "${CMAKE_BINARY_DIR}/volk_gnsssdr_profile.1.gz"
        )

        execute_process(COMMAND gzip -9 -c ${CMAKE_SOURCE_DIR}/src/algorithms/libs/volk_gnsssdr_module/volk_gnsssdr/cmake/Packaging/volk_gnsssdr-config-info-manpage
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            OUTPUT_FILE "${CMAKE_BINARY_DIR}/volk_gnsssdr-config-info.1.gz"
        )

        install(FILES ${CMAKE_BINARY_DIR}/volk_gnsssdr_profile.1.gz
            DESTINATION share/man/man1
        )

        install(FILES ${CMAKE_BINARY_DIR}/volk_gnsssdr-config-info.1.gz
            DESTINATION share/man/man1
        )
    endif()
endif()

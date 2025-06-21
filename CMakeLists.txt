cmake_minimum_required(VERSION 3.15)

project(SoLoud
    VERSION 1.0.3
    DESCRIPTION "SoLoud audio library"
    LANGUAGES C CXX)

# set c++ standard
set(CMAKE_CXX_STANDARD 17) # I'm not sure how to make 17 the default but still allow 11, since the project __is__ compatible with C++11.
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON) # will this cause problems?

option(BUILD_SHARED_LIBS "Build using shared libraries" OFF)
option(SOLOUD_TESTS "Build test executables" OFF)
option(SOLOUD_ASAN "Enable AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
option(SOLOUD_TSAN "Enable ThreadSanitizer and UndefinedBehaviorSanitizer" OFF)

# check for conflicting sanitizer options
if(SOLOUD_ASAN AND SOLOUD_TSAN)
    message(FATAL_ERROR "SOLOUD_ASAN and SOLOUD_TSAN cannot be enabled simultaneously")
endif()

# set default build type if none specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# find optional dependencies
find_package(PkgConfig QUIET)

# check for sdl3
set(SOLOUD_WITH_SDL3 OFF)
if(PkgConfig_FOUND)
    pkg_check_modules(SDL3 QUIET sdl3>=3.2.0)
    if(SDL3_FOUND)
        set(SOLOUD_WITH_SDL3 ON)
        message(STATUS "Found SDL3 ${SDL3_VERSION}")
    endif()
endif()

# check for libmpg123
set(SOLOUD_WITH_LIBMPG123 OFF)
if(PkgConfig_FOUND)
    pkg_check_modules(MPG123 QUIET libmpg123>=1.30.0)
    if(MPG123_FOUND)
        set(SOLOUD_WITH_LIBMPG123 ON)
        message(STATUS "Found libmpg123 ${MPG123_VERSION}")
    endif()
endif()

# collect source files
file(GLOB_RECURSE SOLOUD_CORE_SOURCES "src/core/*.cpp")
file(GLOB_RECURSE SOLOUD_AUDIOSOURCE_SOURCES "src/audiosource/*.cpp" "src/audiosource/*.c")
file(GLOB_RECURSE SOLOUD_FILTER_SOURCES "src/filter/*.cpp")

# backend sources (always include these)
set(SOLOUD_BACKEND_SOURCES
    src/backend/miniaudio/soloud_miniaudio.cpp
    src/backend/nosound/soloud_nosound.cpp
    src/backend/null/soloud_null.cpp
)

# add sdl3 backend if available
if(SOLOUD_WITH_SDL3)
    list(APPEND SOLOUD_BACKEND_SOURCES src/backend/sdl3/soloud_sdl3.cpp)
endif()

# vendored sources
file(GLOB SOLOUD_VENDORED_SOURCES "vendored/*.c" "vendored/*.cpp")

# combine all sources
set(SOLOUD_SOURCES
    ${SOLOUD_CORE_SOURCES}
    ${SOLOUD_AUDIOSOURCE_SOURCES}
    ${SOLOUD_FILTER_SOURCES}
    ${SOLOUD_BACKEND_SOURCES}
    ${SOLOUD_VENDORED_SOURCES}
)

# create the library target
if(BUILD_SHARED_LIBS)
    add_library(soloud SHARED ${SOLOUD_SOURCES})
    # for shared libraries on msvc, also need to _STATIC_CPPLIB
    if(MSVC)
        target_compile_definitions(soloud PRIVATE _STATIC_CPPLIB)
    endif()
else()
    add_library(soloud STATIC ${SOLOUD_SOURCES})
endif()

# include directories
target_include_directories(soloud
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/vendored
)

# compile definitions
target_compile_definitions(soloud
    PUBLIC
        $<$<CONFIG:Debug>:_DEBUG>
    PRIVATE
        $<$<BOOL:${SOLOUD_WITH_SDL3}>:WITH_SDL3>
        $<$<BOOL:${SOLOUD_WITH_LIBMPG123}>:WITH_LIBMPG123>
)

# msvc specific settings
if(MSVC)
    target_compile_definitions(soloud PRIVATE
        _CRT_SECURE_NO_WARNINGS
        _HAS_EXCEPTIONS=0
        WINDOWS_VERSION
    )
    # remove default cmake exceptions setting
    string(REGEX REPLACE "/EH[a-z]+" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
else()
    if(WIN32)
        target_compile_definitions(soloud PRIVATE
            WINDOWS_VERSION
        )
    endif()
    # disable exceptions for gcc/clang
    target_compile_options(soloud PRIVATE -fno-exceptions)
endif()

# threading
find_package(Threads REQUIRED)
target_link_libraries(soloud PRIVATE Threads::Threads)

# sanitizers
if(SOLOUD_ASAN)
    if(MSVC)
        target_compile_options(soloud PRIVATE /fsanitize=address)
        target_link_options(soloud PRIVATE /fsanitize=address)
    else()
        target_compile_options(soloud PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
        target_link_options(soloud PRIVATE -fsanitize=address,undefined)
    endif()
elseif(SOLOUD_TSAN)
    if(MSVC)
        message(WARNING "ThreadSanitizer is not supported on MSVC")
    else()
        target_compile_options(soloud PRIVATE -fsanitize=thread,undefined)
        target_link_options(soloud PRIVATE -fsanitize=thread,undefined)
    endif()
endif()

# link libraries
if(SOLOUD_WITH_SDL3)
    target_include_directories(soloud PRIVATE ${SDL3_INCLUDE_DIRS})
    target_link_libraries(soloud PRIVATE ${SDL3_LIBRARIES})
endif()

if(SOLOUD_WITH_LIBMPG123)
    target_include_directories(soloud PRIVATE ${MPG123_INCLUDE_DIRS})
    target_link_libraries(soloud PRIVATE ${MPG123_LIBRARIES})
endif()

# platform specific libraries (don't actually link to them, just add the includes)
if(WIN32)
    # i don't think we need anything specifically?
elseif(APPLE)
    # macos core audio for miniaudio
    find_library(AUDIOTOOLBOX_FRAMEWORK AudioToolbox)
    if(AUDIOTOOLBOX_FRAMEWORK)
        target_include_directories(soloud PRIVATE ${AUDIOTOOLBOX_FRAMEWORK_INCLUDE_DIRS})
    endif()
elseif(UNIX)
    # -ldl (for miniaudio to find required libraries at runtime)
    target_link_libraries(soloud PRIVATE ${CMAKE_DL_LIBS})

    # alsa/pulse for miniaudio
    if(PkgConfig_FOUND)
        pkg_check_modules(ALSA QUIET alsa)
        if(ALSA_FOUND)
            target_include_directories(soloud PRIVATE ${ALSA_INCLUDE_DIRS})
        endif()
        pkg_check_modules(PULSEAUDIO QUIET libpulse)
        if(PULSEAUDIO_FOUND)
            target_include_directories(soloud PRIVATE ${PULSEAUDIO_INCLUDE_DIRS})
        endif()
    endif()
endif()

# c api shared library (CURRENTLY OUTDATED, SHOULD BE RE-GENERATED!)
if(BUILD_SHARED_LIBS)
    add_library(soloud_c SHARED src/c_api/soloud_c.cpp)
    target_link_libraries(soloud_c PRIVATE soloud)
    target_include_directories(soloud_c
        PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:include>
    )

    # export definitions for windows dll
    if(WIN32)
        set_target_properties(soloud_c PROPERTIES
            WINDOWS_EXPORT_ALL_SYMBOLS ON
        )
    endif()
endif()

# test executables
if(SOLOUD_TESTS)
    file(GLOB TEST_DIRS "src/tests/*")
    foreach(TEST_DIR ${TEST_DIRS})
        if(IS_DIRECTORY ${TEST_DIR})
            get_filename_component(TEST_NAME ${TEST_DIR} NAME)
            file(GLOB_RECURSE TEST_SOURCES "${TEST_DIR}/*.cpp" "${TEST_DIR}/*.c")

            if(TEST_SOURCES)
                add_executable(${TEST_NAME} ${TEST_SOURCES})
                target_link_libraries(${TEST_NAME} PRIVATE soloud)
                target_include_directories(${TEST_NAME} PRIVATE
                    ${CMAKE_CURRENT_SOURCE_DIR}/include
                    ${CMAKE_CURRENT_SOURCE_DIR}/src
                )

                # apply same compiler settings as main library
                if(MSVC)
                    target_compile_definitions(${TEST_NAME} PRIVATE
                        _CRT_SECURE_NO_WARNINGS
                        _HAS_EXCEPTIONS=0
                    )
                    string(REGEX REPLACE "/EH[a-z]+" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
                else()
                    target_compile_options(${TEST_NAME} PRIVATE -fno-exceptions)
                endif()

                # apply sanitizers if enabled
                if(SOLOUD_ASAN)
                    if(MSVC)
                        target_compile_options(${TEST_NAME} PRIVATE /fsanitize=address)
                        target_link_options(${TEST_NAME} PRIVATE /fsanitize=address)
                    else()
                        target_compile_options(${TEST_NAME} PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
                        target_link_options(${TEST_NAME} PRIVATE -fsanitize=address,undefined)
                    endif()
                elseif(SOLOUD_TSAN)
                    if(MSVC)
                        message(WARNING "ThreadSanitizer is not supported on MSVC")
                    else()
                        target_compile_options(${TEST_NAME} PRIVATE -fsanitize=thread,undefined)
                        target_link_options(${TEST_NAME} PRIVATE -fsanitize=thread,undefined)
                    endif()
                endif()
            endif()
        endif()
    endforeach()
endif()

# installation
include(GNUInstallDirs)

# .pc generation
set(SOLOUD_PC_REQUIRES "")
set(SOLOUD_PC_REQUIRES_PRIVATE "")
set(SOLOUD_PC_LIBS_PRIVATE ${CMAKE_THREAD_LIBS_INIT})

# add optional dependencies to requires
if(SOLOUD_WITH_SDL3)
    set(SOLOUD_PC_REQUIRES_PRIVATE "${SOLOUD_PC_REQUIRES_PRIVATE} sdl3")
endif()

if(SOLOUD_WITH_LIBMPG123)
    set(SOLOUD_PC_REQUIRES_PRIVATE "${SOLOUD_PC_REQUIRES_PRIVATE} libmpg123")
endif()

# only add libraries that are actually linked
if(UNIX AND NOT APPLE)
    # only -ldl is actually linked on unix systems
    set(SOLOUD_PC_LIBS_PRIVATE "${SOLOUD_PC_LIBS_PRIVATE} -ldl")
endif()

# trim leading spaces
string(STRIP "${SOLOUD_PC_REQUIRES}" SOLOUD_PC_REQUIRES)
string(STRIP "${SOLOUD_PC_REQUIRES_PRIVATE}" SOLOUD_PC_REQUIRES_PRIVATE)
string(STRIP "${SOLOUD_PC_LIBS_PRIVATE}" SOLOUD_PC_LIBS_PRIVATE)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/soloud.pc.in"
    "${CMAKE_CURRENT_BINARY_DIR}/soloud.pc"
    @ONLY
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/soloud.pc"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
)

install(TARGETS soloud
    EXPORT SoLoudTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

if(BUILD_SHARED_LIBS AND TARGET soloud_c)
    install(TARGETS soloud_c
        EXPORT SoLoudTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.h"
)

# export configuration
install(EXPORT SoLoudTargets
    FILE SoLoudTargets.cmake
    NAMESPACE SoLoud::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SoLoud
)

# create config file
include(CMakePackageConfigHelpers)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/SoLoudConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/SoLoudConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SoLoud
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/SoLoudConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/SoLoudConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/SoLoudConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SoLoud
)

# display configuration summary
message(STATUS "SoLoud configuration:")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  Build shared: ${BUILD_SHARED_LIBS}")
message(STATUS "  Build tests: ${SOLOUD_TESTS}")
message(STATUS "  AddressSanitizer: ${SOLOUD_ASAN}")
message(STATUS "  ThreadSanitizer: ${SOLOUD_TSAN}")
message(STATUS "  SDL3 support: ${SOLOUD_WITH_SDL3}")
message(STATUS "  libmpg123 support: ${SOLOUD_WITH_LIBMPG123}")

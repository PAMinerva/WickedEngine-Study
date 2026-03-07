cmake_minimum_required(VERSION 3.12)
project(16-GPUCloth)

file(GLOB HEADER_FILES CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/incl/*.h
)

file(GLOB SRC_FILES CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
)

file(GLOB MAIN_FILES CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/main_*.cpp
)

set(SOURCE_FILES ${HEADER_FILES} ${SRC_FILES} ${MAIN_FILES})

list(FILTER SOURCE_FILES EXCLUDE REGEX ${CMAKE_CURRENT_SOURCE_DIR}/main_.*)
list(APPEND SOURCE_FILES main_${PLATFORM}.cpp)

if (WIN32)

	list (APPEND SOURCE_FILES
            ${CMAKE_CURRENT_SOURCE_DIR}/resource.rc
	)

	add_executable(${PROJECT_NAME} WIN32 ${SOURCE_FILES})

    add_compile_definitions(WIN32 _WINDOWS)

	set_property(TARGET ${PROJECT_NAME} PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
    set_target_properties(${PROJECT_NAME} PROPERTIES WIN32_EXECUTABLE YES)

	set(LIB_DXCOMPILER "dxcompiler.dll")

else ()

	add_executable(${PROJECT_NAME} ${SOURCE_FILES})

	set(LIB_DXCOMPILER "libdxcompiler.so")

endif ()

target_link_libraries(${PROJECT_NAME} PUBLIC WickedEngine)

target_include_directories(${PROJECT_NAME} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/incl
    ${CMAKE_CURRENT_SOURCE_DIR}/shaders
)

if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
endif()

# Pass source shader directory to C++ so the engine reads .hlsl from source tree directly
target_compile_definitions(${PROJECT_NAME} PRIVATE
    SHADER_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/shaders/"
)

set(LIBDXCOMPILER_PATH "${WICKED_ROOT_DIR}/WickedEngine/${LIB_DXCOMPILER}")
set(STARTUP_LUA "${WICKED_ROOT_DIR}/Editor/startup.lua")

add_custom_command(
    TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${LIBDXCOMPILER_PATH} ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${STARTUP_LUA} ${CMAKE_CURRENT_BINARY_DIR}/startup.lua
)

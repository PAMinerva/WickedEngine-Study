cmake_minimum_required(VERSION 3.12)
project(01-HelloTriangle)

file(GLOB SOURCE_FILES CONFIGURE_DEPENDS *.cpp *.h)
list(FILTER SOURCE_FILES EXCLUDE REGEX ${CMAKE_CURRENT_SOURCE_DIR}/main_.*)
list(APPEND SOURCE_FILES main_${PLATFORM}.cpp)

if (WIN32)

	list (APPEND SOURCE_FILES
            ${CMAKE_CURRENT_SOURCE_DIR}/resource.rc
	)

    # Quando si utilizza WIN32 come argomento in add_executable, si indica a CMake che l'eseguibile target è un'applicazione Windows GUI. 
    # Questo ha i seguenti effetti:
    #   - CMake imposterà l'opzione WIN32_EXECUTABLE per l'eseguibile target, che indica al compilatore di creare un'applicazione Windows GUI,
    #     sopprimendo la finestra della console che altrimenti verrebbe visualizzata quando si esegue l'applicazione.
    #   - Il punto di ingresso predefinito è WinMain, invece che main.
	add_executable(${PROJECT_NAME} WIN32 ${SOURCE_FILES})

	target_link_libraries(${PROJECT_NAME} PUBLIC
            WickedEngine_Windows
	)

    add_compile_definitions(WIN32 _WINDOWS)

	set_property(TARGET ${PROJECT_NAME} PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
    set_target_properties(${PROJECT_NAME} PROPERTIES WIN32_EXECUTABLE YES)

	set(LIB_DXCOMPILER "dxcompiler.dll")

else ()

	add_executable(${PROJECT_NAME} ${SOURCE_FILES})

	target_link_libraries(${PROJECT_NAME} PUBLIC WickedEngine_Linux)

	set(LIB_DXCOMPILER "libdxcompiler.so")

endif ()

if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
endif()

set(LIBDXCOMPILER_PATH "${WICKED_ROOT_DIR}/WickedEngine/${LIB_DXCOMPILER}")
set(STARTUP_LUA "${WICKED_ROOT_DIR}/Editor/startup.lua")
message("dxcompiler found at ${LIBDXCOMPILER_PATH}")
message("startup lua found at ${STARTUP_LUA}")

# Copy shader compiler library in build/Samples/Tutorial/01-Hello-Triangle/[Debug, Release]
add_custom_command(
    TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${LIBDXCOMPILER_PATH} ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${STARTUP_LUA} ${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>/startup.lua
)

cmake_minimum_required(VERSION 3.19)

# CMAKE_SOURCE_DIR: directory di CMakeList.txt radice del progetto
# CMAKE_BINARY_DIR: directory dove CMake genera i file di build (in genere build, se esegui cmake . -B build)
if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR
        "In-source builds are not supported!\n"
        "Run `git clean -d -f` to clean up the files CMake has created (stash "
        "your changes first, if you have made any), then run `cmake -B build "
        "<other_options>` followed by `cmake --build build --parallel`"
    )
endif()

project(WickedEngine)

include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported)

# CMAKE_CURRENT_SOURCE_DIR: directory di CMakeList.txt corrente
# (in questo caso coincide con quello radice)
# set memorizza una variabile, in questo caso la directory radice del progetto
set(WICKED_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# option: fornisce una opzione booleana che può essere attivata dall'utente
option(WICKED_DYNAMIC_LIBRARY "Build WickedEngine as a dynamic library" OFF)
option(WICKED_PIC "Build WickedEngine as position-independent code" WICKED_DYNAMIC_LIBRARY)
option(USE_LIBCXX "Link WickedEngine to llvm libc++ library - only available with the Clang compiler" OFF)
option(WICKED_EDITOR "Build WickedEngine editor" ON)
option(WICKED_TESTS "Build WickedEngine tests" ON)
option(WICKED_IMGUI_EXAMPLE "Build WickedEngine imgui example" ON)
option(WICKED_USE_IPO "Enable IPO/LTO in non-debug builds" ${ipo_supported})

if (CMAKE_HOST_WIN32)
    set(symlink_default OFF)
else()
    set(symlink_default ON)
endif()
option(WICKED_USE_SYMLINKS "Prefer symlinking over copying directories" ${symlink_default})

# check that IPO is supported when turned on
if(WICKED_USE_IPO)
    check_ipo_supported()
endif()

if(WICKED_USE_SYMLINKS)
    # check for symlink support (on windows it requires admin or developer mode)
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E create_symlink CMakeLists.txt ${CMAKE_CURRENT_BINARY_DIR}/symlink-test
        ERROR_QUIET
        RESULT_VARIABLE symlink_check_result
    )
    if(symlink_check_result EQUAL 0)
        file(REMOVE ${CMAKE_BINARY_DIR}/symlink-test)
        set(COPY_OR_SYMLINK_DIR_CMD create_symlink)
    else()
        message(FATAL_ERROR "Symlinks are not supported. Either disable them using -DWICKED_USE_SYMLINKS=OFF or enable Windows' developer mode")
    endif()
else()
    if(CMAKE_VERSION VERSION_LESS "3.26.0")
        set(COPY_OR_SYMLINK_DIR_CMD copy_directory)
    else()
        set(COPY_OR_SYMLINK_DIR_CMD copy_directory_if_different)
    endif()
endif()

# option(WICKED_TUTORIAL_SAMPLES "Build tutorial samples using WickedEngine library" ON)

# include il modulo CMakeDependentOption, che consente di definire opzioni che dipendono dal valore di altre variabili o opzioni.
# In questo caso forse non era necessario, dato che opzioni in blocchi if o elseif non necessitano di tale modulo anche se in effetti
# dipendono da variabili usate come consizioni di if ed elseif. Comunque è tutto da verificare meglio.
# include(CMakeDependentOption)

if(UNIX)
    option(WICKED_LINUX_TEMPLATE "Build WickedEngine Linux template" ON)
elseif(WIN32)
    option(WICKED_WINDOWS_TEMPLATE "Build WickedEngine Windows template" ON)
endif()

# Configure CMake global variables
# Standard del linguaggio C++:
# Imposta lo standard minimo a C++17.
# Richiede obbligatoriamente C++17.
# Disabilita le estensioni specifiche del compilatore.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE WICKED_PIC)

# Use solution folders to organize projects
# Crea od usa cartelle virtuali visibili solo nell'ID (solo VS, no VScode) e per posizionare i file generati (da verificare)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)


# Definisce il nome del progetto
project(WickedEngine)

if (WIN32)
    set(PLATFORM "Windows")
    add_compile_definitions(WIN32=1) # aggiunge definizione da usare in fase di compilazione
    # add_compile_definitions(_WIN32=1) this is a given from the compiler
elseif(UNIX)
    set(PLATFORM "SDL2")
    add_compile_definitions(SDL2=1)

    # Common compiler options and warning level for CLANG and GCC:
    add_compile_options(
        -Wall

        -Wno-unused-variable
        -Wno-unused-function
        -Wno-unused-but-set-variable
        -Wno-sign-compare
    )
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    # CLANG specific compile options:
    add_compile_options(
        -fdeclspec
        -fms-extensions

        -Wno-nullability-completeness
        -Wno-unused-private-field

        # Increase Ccache hit rate, note that you also need to set config options
        # in your Ccache config, see https://ccache.dev/manual/latest.html#_precompiled_headers
        "SHELL:-Xclang -fno-pch-timestamp"
    )
    if (USE_LIBCXX)
        add_compile_options(
            $<$<COMPILE_LANGUAGE:CXX>:-stdlib=libc++>
        )
        add_link_options(
            -stdlib=libc++
        )
    endif()
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # GCC specific compile options:
    add_compile_options(
        -Wno-strict-aliasing
    )
endif()


# add_subdirectory aggiunge una sottodirectory al progetto, permettendo di gestirla come parte della build principale.
# In questo caso CMake cerca il file CMakeLists.txt in sottocartella WickedEngine.
# Configura i target definiti in WickedEngine/CMakeLists.txt.
# I file risultanti (eseguibili, librerie, ecc.) saranno generati in build/WickedEngine.
add_subdirectory(WickedEngine)
add_custom_target(Content
   COMMAND ${CMAKE_COMMAND} -E ${COPY_OR_SYMLINK_DIR_CMD} ${WICKED_ROOT_DIR}/Content ${CMAKE_CURRENT_BINARY_DIR}/Content
   COMMENT "$<IF:$<BOOL:${WICKED_USE_SYMLINKS}>,Symlinking,Copying> Content directory"
   VERBATIM
)

if (WICKED_EDITOR)
    add_subdirectory(Editor)
endif()

if (WICKED_TESTS)
    add_subdirectory(Samples/Tests)
endif()

if (WICKED_IMGUI_EXAMPLE)
    add_subdirectory(Samples/Example_ImGui)
    add_subdirectory(Samples/Example_ImGui_Docking)
endif()

if (WICKED_LINUX_TEMPLATE)
    add_subdirectory(Samples/Template_Linux)
endif()

if (WICKED_WINDOWS_TEMPLATE)
    add_subdirectory(Samples/Template_Windows)
endif()

# My samples
add_subdirectory(Samples/Tutorial/01-HelloTriangle)

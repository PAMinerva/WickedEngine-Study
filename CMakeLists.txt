cmake_minimum_required(VERSION 3.12)

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
option(WICKED_TUTORIAL_SAMPLES "Build tutorial samples using WickedEngine library" ON)

# include il modulo CMakeDependentOption, che consente di definire opzioni che dipendono dal valore di altre variabili o opzioni.
# In questo caso forse non era necessario, dato che opzioni in blocchi if o elseif non necessitano di tale modulo anche se in effetti
# dipendono da variabili usate come consizioni di if ed elseif. Comunque è tutto da verificare meglio.
include(CMakeDependentOption)

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
    set(DXC_TARGET "${CMAKE_CURRENT_SOURCE_DIR}/WickedEngine/dxc.exe")
elseif(UNIX)
    set(PLATFORM "SDL2")
    add_compile_definitions(SDL2=1)
    set(DXC_TARGET "dxc")
endif()


if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdeclspec -fms-extensions")
	if (USE_LIBCXX)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
		set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libc++ -lc++abi")
	endif()
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
else()
endif()


# add_subdirectory aggiunge una sottodirectory al progetto, permettendo di gestirla come parte della build principale.
# In questo caso CMake cerca il file CMakeLists.txt in sottocartella WickedEngine.
# Configura i target definiti in WickedEngine/CMakeLists.txt.
# I file risultanti (eseguibili, librerie, ecc.) saranno generati in build/WickedEngine.
add_subdirectory(WickedEngine)
# simile per sottocartella Content
add_subdirectory(Content)

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
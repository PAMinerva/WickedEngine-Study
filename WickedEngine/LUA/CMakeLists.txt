cmake_minimum_required(VERSION 3.19)
file(GLOB LUA_HEADERS CONFIGURE_DEPENDS *.h)
file(GLOB LUA_SOURCES CONFIGURE_DEPENDS *.c)
list(FILTER LUA_SOURCES EXCLUDE REGEX lua\.c)

add_library(LUA STATIC
    ${LUA_SOURCES}
    ${LUA_HEADERS}
)

if (UNIX)
    # for Linux, BSD, Solaris, Minix, Apple
    target_compile_definitions(LUA PUBLIC LUA_USE_POSIX=1)
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_compile_options(LUA PRIVATE
        -Wno-string-plus-int
    )
endif()

# La proprietà FOLDER è usata solo per organizzare i target in IDE con un'interfaccia grafica, come Visual Studio. 
# Non influisce sul processo di build. "ThirdParty" è il nome della cartella virtuale che sarà visibile nell'IDE.
# Ad esempio, in Visual Studio, il target LUA apparirà sotto una cartella chiamata ThirdParty nell'elenco dei progetti.
# La proprietà POSITION_INDEPENDENT_CODE (spesso abbreviata come PIC) viene utilizzata per indicare a CMake se 
# il codice generato per questo target è position-independent. 
# Questo significa che il codice prodotto può essere eseguito correttamente indipendentemente dalla posizione in 
# memoria in cui viene caricato (tramite offsets).
set_target_properties(LUA PROPERTIES
    FOLDER "ThirdParty"
    POSITION_INDEPENDENT_CODE ${WICKED_PIC}
)

install(FILES ${LUA_HEADERS}
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/WickedEngine/LUA/")

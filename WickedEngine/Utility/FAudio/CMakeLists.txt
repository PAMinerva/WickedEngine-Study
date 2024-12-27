# Based on upstream's CMakeLists.txt
# Changes compared to upstream:
# * cmake_minimum_required adjusted
# * unused options removed
# * windows related stuff removed
# * install rule added

# CMake Project for FAudio
# Written by @NeroBurner
cmake_minimum_required(VERSION 3.19)
project(FAudio C)

# Options
# Definisce un'opzione per specificare se costruire FAudio come libreria condivisa (shared) o statica. In questo caso è OFF di default, 
# quindi viene creata una libreria statica.
# L'utente può cambiarla passando -DBUILD_SHARED_LIBS=ON al comando cmake.
option(BUILD_SHARED_LIBS "Build shared library" OFF)
option(BUILD_SDL3 "Build against SDL 3.0" OFF)

# C99 Requirement
# Specifica che il codice deve essere compilato usando lo standard C99, 
# disabilitando estensioni proprietarie specifiche del compilatore
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_EXTENSIONS OFF)

# Version
SET(LIB_MAJOR_VERSION "0")
SET(LIB_MINOR_VERSION "25")
SET(LIB_REVISION "01")
SET(LIB_VERSION "${LIB_MAJOR_VERSION}.${LIB_MINOR_VERSION}.${LIB_REVISION}")

# Build Type
# Se l'utente non specifica un tipo di build (Debug, Release, etc.), imposta Release come default.
# Se l'utente esegue cmake .. -DCMAKE_BUILD_TYPE=Release[Debug], allora viene usato il tipo di build specificato.
# CACHE indica che la variabile deve essere memorizzata nella cache di CMake.
# STRING specifica il tipo di dato della variabile nella cache, seguito dalla stringa che descrive la variabile,
# che può essere così visualizzata in strumenti come cmake-gui.
# FORCE forza l'aggiornamento della variabile nella cache, anche se è già stata impostata in precedenza.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	# By default, we use Release
	message(STATUS "Setting build type to 'Release' as none was specified.")
	set(CMAKE_BUILD_TYPE "Release" CACHE
		STRING "Choose the type of build." FORCE
	)

	# Set the possible values of build type for cmake-gui
	# CACHE specifica che la proprietà deve essere impostata su una variabile nella cache di CMake.
	# CMAKE_BUILD_TYPE: Nome della variabile di cui si sta impostando la proprietà.
	# PROPERTY STRINGS: Specifica i valori possibili per la variabile.
	set_property(CACHE CMAKE_BUILD_TYPE PROPERTY
		STRINGS "Debug" "Release" "RelWithDebInfo"
	)
endif()

set(PUBLIC_HEADERS
	include/F3DAudio.h
	include/FACT3D.h
	include/FACT.h
	include/FAPOBase.h
	include/FAPOFX.h
	include/FAPO.h
	include/FAudioFX.h
	include/FAudio.h
	)

# Source lists
add_library(FAudio STATIC
	# Public Headers
	${PUBLIC_HEADERS}
	# Internal Headers
	src/FACT_internal.h
	src/FAudio_internal.h
	src/stb.h
	src/stb_vorbis.h
	src/qoa_decoder.h
	# Source Files
	src/F3DAudio.c
	src/FACT3D.c
	src/FACT.c
	src/FACT_internal.c
	src/FAPOBase.c
	src/FAPOFX.c
	src/FAPOFX_echo.c
	src/FAPOFX_eq.c
	src/FAPOFX_masteringlimiter.c
	src/FAPOFX_reverb.c
	src/FAudio.c
	src/FAudioFX_reverb.c
	src/FAudioFX_volumemeter.c
	src/FAudio_internal.c
	src/FAudio_internal_simd.c
	src/FAudio_operationset.c
	src/FAudio_platform_sdl2.c
	src/FAudio_platform_sdl3.c
	src/FAudio_platform_win32.c
	# Optional source files
	src/XNA_Song.c
)

install(FILES ${PUBLIC_HEADERS}
		DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/WickedEngine/FAudio/")

# FAudio folders as includes, for other targets to consume
# src/stb.h sopra non significa niente. Il compilatore dev e sapere esattamente il percorso completo degli headers.
# target_include_directories: Comando CMake utilizzato per specificare le directory di inclusione per un target.
# PUBLIC: Indica che le directory di inclusione sono necessarie sia per compilare il target FAudio sia per i target che dipendono da FAudio.
# $<BUILD_INTERFACE:...>: Generatore di espressioni CMake che specifica le directory di inclusione da utilizzare durante la fase di build.
# $<INSTALL_INTERFACE:...>: Generatore di espressioni CMake che specifica le directory di inclusione da utilizzare durante la fase di installazione.
# La directory di inclusione specificata durante la fase di installazione serve a garantire che i file di intestazione (header files) necessari 
# per utilizzare la libreria siano disponibili per altri progetti che desiderano utilizzare la libreria installata.
# Quando un progetto viene installato, i file di intestazione devono essere copiati in una posizione standard all'interno della 
# directory di installazione. Questo permette ad altri progetti di trovare e includere questi file di intestazione quando compilano 
# il loro codice che dipende dalla libreria installata.
# In questo caso, il comando target_include_directories sotto specifica che durante la fase di installazione, 
# i file di intestazione saranno disponibili nella directory `${CMAKE_INSTALL_INCLUDEDIR}/WickedEngine/FAudio`. 
# Questo significa che quando la libreria `FAudio` viene installata, i file di intestazione saranno copiati in questa directory, rendendoli 
# accessibili ad altri progetti che includono `FAudio`.
target_include_directories(FAudio PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/WickedEngine/FAudio>
)

# Soname
# impostando la proprietà OUTPUT_NAME di un target puoi modificare il nome dell'output impostato in precedenza da add_library.
# SOVERSION: Proprietà che specifica la versione sovraccaricata (soname) del target.
set_target_properties(FAudio PROPERTIES
	OUTPUT_NAME "FAudio"
	VERSION ${LIB_VERSION}
	SOVERSION ${LIB_MAJOR_VERSION}
)

# SDL Dependency
if (PLATFORM_WIN32)
	message(STATUS "not using SDL")
elseif (BUILD_SDL3)
	   if (DEFINED SDL3_INCLUDE_DIRS AND DEFINED SDL3_LIBRARIES)
			   message(STATUS "using pre-defined SDL3 variables SDL3_INCLUDE_DIRS and SDL3_LIBRARIES")
		target_include_directories(FAudio PUBLIC "$<BUILD_INTERFACE:${SDL3_INCLUDE_DIRS}>")
		target_link_libraries(FAudio PUBLIC ${SDL3_LIBRARIES})
			   if(INSTALL_MINGW_DEPENDENCIES)
					   install_shared_libs(${SDL3_LIBRARIES} DESTINATION bin NO_INSTALL_SYMLINKS)
			   endif()
	   else()
			   # Only try to autodetect if both SDL3 variables aren't explicitly set
			   find_package(SDL3 CONFIG)
			   if (TARGET SDL3::SDL3)
					   message(STATUS "using TARGET SDL3::SDL3")
					   target_link_libraries(FAudio PUBLIC SDL3::SDL3)
					   if(INSTALL_MINGW_DEPENDENCIES)
							   install_shared_libs(TARGETS SDL3::SDL3 DESTINATION bin NO_INSTALL_SYMLINKS REQUIRED)
					   endif()
			   elseif (TARGET SDL3)
					   message(STATUS "using TARGET SDL3")
					   target_link_libraries(FAudio PUBLIC SDL3)
					   if(INSTALL_MINGW_DEPENDENCIES)
							   install_shared_libs(TARGETS SDL3 DESTINATION bin NO_INSTALL_SYMLINKS REQUIRED)
					   endif()
			   else()
					   message(STATUS "no TARGET SDL3::SDL3, or SDL3, using variables")
			target_include_directories(FAudio PUBLIC "$<BUILD_INTERFACE:${SDL3_INCLUDE_DIRS}>")
			target_link_libraries(FAudio PUBLIC ${SDL3_LIBRARIES})
					   if(INSTALL_MINGW_DEPENDENCIES)
							   install_shared_libs(${SDL3_LIBRARIES} DESTINATION bin NO_INSTALL_SYMLINKS)
					   endif()
			   endif()
	   endif()
elseif (DEFINED SDL2_INCLUDE_DIRS AND DEFINED SDL2_LIBRARIES)
	message(STATUS "using pre-defined SDL2 variables SDL2_INCLUDE_DIRS and SDL2_LIBRARIES")
	target_include_directories(FAudio PUBLIC "$<BUILD_INTERFACE:${SDL2_INCLUDE_DIRS}>")
	target_link_libraries(FAudio PUBLIC ${SDL2_LIBRARIES})
	if(INSTALL_MINGW_DEPENDENCIES)
		install_shared_libs(${SDL2_LIBRARIES} DESTINATION bin NO_INSTALL_SYMLINKS)
	endif()
else()
	# Only try to autodetect if both SDL2 variables aren't explicitly set
	find_package(SDL2 CONFIG)
	if (TARGET SDL2::SDL2)
		message(STATUS "using TARGET SDL2::SDL2")
		target_link_libraries(FAudio PUBLIC SDL2::SDL2)
		if(INSTALL_MINGW_DEPENDENCIES)
			install_shared_libs(TARGETS SDL2::SDL2 DESTINATION bin NO_INSTALL_SYMLINKS REQUIRED)
		endif()
	elseif (TARGET SDL2)
		message(STATUS "using TARGET SDL2")
		target_link_libraries(FAudio PUBLIC SDL2)
		if(INSTALL_MINGW_DEPENDENCIES)
			install_shared_libs(TARGETS SDL2 DESTINATION bin NO_INSTALL_SYMLINKS REQUIRED)
		endif()
	else()
		message(STATUS "no TARGET SDL2::SDL2, or SDL2, using variables")
		target_include_directories(FAudio PUBLIC "$<BUILD_INTERFACE:${SDL2_INCLUDE_DIRS}>")
		target_link_libraries(FAudio PUBLIC ${SDL2_LIBRARIES})
		if(INSTALL_MINGW_DEPENDENCIES)
			install_shared_libs(${SDL2_LIBRARIES} DESTINATION bin NO_INSTALL_SYMLINKS)
		endif()
	endif()
endif()

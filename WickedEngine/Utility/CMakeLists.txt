cmake_minimum_required(VERSION 3.19)

# CMAKE_INSTALL_INCLUDEDIR è una variabile definita ed esposta dal modulo GNUInstallDirs
include(GNUInstallDirs)
if (PLATFORM MATCHES "SDL2")
	add_subdirectory(FAudio)
endif()

# GLOB_RECURSE: fa sì che il comando cerchi ricorsivamente i file che corrispondono ai pattern specificati (come *.h, *.hpp, e *.inl) 
# nella directory corrente e in tutte le sue sottodirectory.
# CONFIGURE_DEPENDS:Aggiorna automaticamente l'elenco dei file ogni volta che nuovi file vengono aggiunti alla directory 
# (durante la fase di configurazione di CMake).
#
# Scenario senza CONFIGURE_DEPENDS:
#
# Supponiamo di avere questo in un CMakeLists.txt:# 
# 	file(GLOB HEADER_FILES *.h)
# 	file(GLOB SOURCE_FILES *.cpp)
# 
# 	add_executable(MyApp ${HEADER_FILES} ${SOURCE_FILES})
#
# Se esegui:
# 
# 	cmake .
# 	make
#
# Cosa succede?
#
# Fase di configurazione (cmake .):
# 
# 	- file(GLOB *.h *.cpp) raccoglie i file .h e .cpp presenti nella directory al momento dell'esecuzione del comando.
# 	- Genera un build system (ad esempio un Makefile).
# 	- Modifica dopo la configurazione:
# 
# Supponiamo che tu crei un nuovo file newfile.cpp nella directory.
# Se esegui direttamente make, il nuovo file non sarà incluso nella compilazione. 
# Questo perché il build system non è stato rigenerato per includere il nuovo file.
# Soluzione senza CONFIGURE_DEPENDS: Devi eseguire manualmente di nuovo cmake . per aggiornare il build system.
#
# Scenario con CONFIGURE_DEPENDS:
#
# Modifichiamo il codice in questo modo:
# 
# 	file(GLOB HEADER_FILES CONFIGURE_DEPENDS *.h)
# 	file(GLOB SOURCE_FILES CONFIGURE_DEPENDS *.cpp)
# 
# 	add_executable(MyApp ${HEADER_FILES} ${SOURCE_FILES})
#
# Ora esegui:
# 
# 	cmake .
# 	make
#
# Cosa succede?
#
# Fase di configurazione (cmake .):
# 
# 	- CONFIGURE_DEPENDS dice a CMake di monitorare i file nella directory che corrispondono al pattern *.h e *.cpp.
# 	- Ogni volta che un file viene aggiunto, modificato o rimosso, CMake aggiornerà automaticamente il build system quando necessario.
#
# Modifica dopo la configurazione:
# 
# Supponiamo che tu aggiunga newfile.cpp nella directory.
# Ora, quando esegui make, il build system viene rigenerato automaticamente per includere newfile.cpp.
#
# Limiti di CONFIGURE_DEPENDS
# Non monitora l'intera directory: Solo i file corrispondenti al pattern specificato (ad esempio *.cpp) vengono monitorati.
# Non rimpiazza cmake . in tutti i casi: Se cambi il pattern o aggiungi una nuova istruzione file(GLOB ...), 
# devi comunque eseguire manualmente cmake .
file(GLOB_RECURSE HEADER_FILES CONFIGURE_DEPENDS
	*.h
	*.hpp
	*.inl
)
file(GLOB_RECURSE SOURCE_FILES CONFIGURE_DEPENDS
	*.c
	*.cpp
)

# workaround for source dir containing regex meta characters
string(REGEX REPLACE "([.+?*])" "\\\\\\1" SDIR "${CMAKE_CURRENT_SOURCE_DIR}")

# FAudio has their own CMakeLists
list(FILTER HEADER_FILES
	EXCLUDE REGEX "${SDIR}/FAudio/.*"
)
list(FILTER SOURCE_FILES
	EXCLUDE REGEX "${SDIR}/FAudio/.*"
)

install(FILES ${HEADER_FILES}
	DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/WickedEngine/Utility/")

# some *.cpp and *.c files are included by other files,
# so they are not source files that should be
# compiled on their own; hence we mark them as header files
# Seleziona i file che non devono essere compilati come unità indipendenti, perché sono inclusi direttamente in altri file sorgenti
file(GLOB_RECURSE DIRECTLY_INCLUDED_FILES CONFIGURE_DEPENDS
	# included by utility_common.cpp
	mikktspace.c
	zstd.c

	# included by volk.h
	volk.c
)

# set_property serve per impostare una proprietà su uno o più elementi, come target, sorgenti, file, directory, ecc. 
# SOURCE: specifica che la proprietà che segue verrà applicata ai file sorgenti elencati in ${DIRECTLY_INCLUDED_FILES}.
# La proprietà HEADER_FILE_ONLY è una proprietà speciale per i file sorgenti in CMake.
# Se è impostata su ON, dice a CMake di non compilare direttamente questi file come parte di un target, 
# anche se sono inclusi nell'elenco dei file sorgenti.
set_property(
	SOURCE ${DIRECTLY_INCLUDED_FILES}
	PROPERTY HEADER_FILE_ONLY ON
)

# linux doesn't support it, and windows includes
# it directly in wiGraphicsDevice_DX12.cpp
list(REMOVE_ITEM SOURCE_FILES
	${CMAKE_CURRENT_SOURCE_DIR}/D3D12MemAlloc.cpp
)
add_library(Utility STATIC
	${SOURCE_FILES}
	${HEADER_FILES}
)

set_target_properties("Utility" PROPERTIES
	FOLDER "ThirdParty"
	POSITION_INDEPENDENT_CODE ${WICKED_PIC}
)

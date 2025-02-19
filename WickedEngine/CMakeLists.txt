cmake_minimum_required(VERSION 3.12)

set(MIN_OpenImageDenoise_VERSION "2.0")

# Il modulo CMakePackageConfigHelpers viene utilizzato per creare file di configurazione che rendono una libreria o un progetto esportabile.
# Questo è particolarmente utile per consentire ad altri progetti di trovare e utilizzare il tuo progetto o libreria tramite CMake.
# I file di configurazione (come WickedEngineConfig.cmake) consentono a chi utilizza il tuo progetto/libreria di:
#    - Trovare il progetto/libreria: specificando dove si trovano i file binari, header, librerie statiche/dinamiche.
#    - Collegare correttamente il progetto: configurando automaticamente le dipendenze necessarie, le directory d’inclusione e altre impostazioni.
# Ad esempio, immagina di voler usare WickedEngine in un altro progetto. Con i file di configurazione corretti, puoi semplicemente scrivere nel tuo CMakeLists.txt:
#    find_package(WickedEngine REQUIRED)
#    target_link_libraries(MyTarget PRIVATE WickedEngine)
# Il file WickedEngineConfig.cmake dice a CMake dove trovare WickedEngine e come collegarlo.
# Riepilogo delle differenze:
# CMakeLists.txt (in sottocartella WickedEngine):
#	Gestisce la configurazione, la costruzione e l'installazione del progetto WickedEngine
#	Contiene istruzioni per compilare e installare il progetto
#	È eseguito durante la configurazione e la compilazione del progetto	
# WickedEngineConfig.cmake (in build/cmake o install)
#	Configura il pacchetto per altri progetti che desiderano usare WickedEngine
#	Contiene variabili e informazioni per "trovare" WickedEngine (con find_package)
#	Viene utilizzato per integrare WickedEngine in altri progetti CMake
# Il modulo GNUInstallDirs, che definisce directory standard di installazione, come:
#    CMAKE_INSTALL_BINDIR: directory per eseguibili (bin).
#    CMAKE_INSTALL_LIBDIR: directory per librerie (lib).
#    CMAKE_INSTALL_INCLUDEDIR: directory per header (include).
# Interessante solo se user esegue cmake --install build. Altrimenti, nulla di ciò che è definito nelle istruzioni install() di questo file avverrà. 
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)
set(INSTALL_LOCAL_CONFIGDIR "${CMAKE_BINARY_DIR}/cmake") # CMAKE_BINARY_DIR: directory di build del progetto (cartella build se si usa cmake . -B build).
set(INSTALL_CONFIGDIR "${CMAKE_INSTALL_LIBDIR}/cmake/WickedEngine")

if (WICKED_DYNAMIC_LIBRARY)
	set(WICKED_LIBRARY_TYPE SHARED)
	message(STATUS "Building WickedEngine as a shared library")
else()
	set(WICKED_LIBRARY_TYPE STATIC)
	message(STATUS "Building WickedEngine as a static library")
endif()

if (WIN32)
	# TODO: Choose whether to use SDL2 on windows as well
	set(TARGET_NAME WickedEngine_Windows)
else ()
	set(TARGET_NAME WickedEngine_Linux)
	# Cerca e configura librerie esterne come SDL2, OpenImageDenoise e Threads.
    # QUIET evita errori visibili se la libreria non viene trovata.
	find_package(SDL2 REQUIRED)
	find_package(OpenImageDenoise "${MIN_OpenImageDenoise_VERSION}" QUIET)
	find_package(Threads REQUIRED)
	if(NOT ${OpenImageDenoise_FOUND})
		message("OpenImageDenoise>=${MIN_OpenImageDenoise_VERSION} not found, it will be disabled.")
	else()
		message("OpenImageDenoise ${OpenImageDenoise_VERSION} Found.")
	endif()

	# CMake incoraggia l'uso di target namespaces, in cui i target appartenenti a un pacchetto o libreria sono prefissati dal nome del pacchetto. Ad esempio:
	# Il pacchetto SDL2 può avere target:
	#	SDL2::SDL2: per la libreria principale di SDL2.
	#	SDL2::SDL2main: per eventuali componenti aggiuntivi.
	#	SDL2::SDL2test: per un'ipotetica libreria di test.
	# Questo approccio aiuta a:
	#	Evitare conflitti di nomi. Ad esempio, potrebbe esistere un target generico chiamato SDL2 in un altro progetto.
	#	Rendere esplicito il legame tra il target e il pacchetto o libreria a cui appartiene.
	# La maggior parte delle librerie moderne, come SDL2, usa il sistema di target esportabili. 
	# Quando utilizzi find_package(SDL2), CMake cerca file di configurazione (es. SDL2Config.cmake) che descrivono i target esportabili definiti dalla libreria stessa. 
	# Questi target seguono la convenzione dei namespaces.
	# Esempio: Il file SDL2Config.cmake esporta un target chiamato SDL2::SDL2. Se hai un progetto che usa SDL2, puoi semplicemente fare:
	#	find_package(SDL2 REQUIRED)
	#	target_link_libraries(MyTarget PRIVATE SDL2::SDL2)
	# Tuttavia, le versioni più vecchie di SDL2 non supportano il moderno approccio di CMake con target predefiniti (ad es. SDL2::SDL2).
	# Questo significa che se SDL2 utilizza un sistema CMake più vecchio, senza target definiti, e si debba creare manualmente un target per garantire il corretto funzionamento.
	# Questo controllo IF verifica se il target SDL2::SDL2 è stato creato da find_package(SDL2 REQUIRED). 
	# Se il target esiste, significa che SDL2 è stata trovata e configurata correttamente. In caso contrario, è necessario creare manualmente il target.
	if(NOT TARGET SDL2::SDL2)
		# using old SDL2 cmake, lets create a SDL2 target ourselves
		# find_library cerca la libreria SDL2 nel sistema. Il risultato è salvato nella variabile temporanea SDL2_LIBRARY_FILE_LOCATION, 
		# che contiene il percorso completo della libreria trovata (es. /usr/lib/libSDL2.so su Linux)
		find_library(SDL2_LIBRARY_FILE_LOCATION SDL2 REQUIRED)

		# add_library crea un target SDL2::SDL2 di tipo imported, indicativo del fatto che non stiamo costruendo SDL2 ma vogliamo usare la libreria già compilata.
		add_library(SDL2::SDL2 SHARED IMPORTED)

		# set_target_properties imposta le proprietà del target SDL2::SDL2. In questo caso, stiamo impostando le directory include e quella delle librerie da collegare.
		# Le variabili come SDL2_INCLUDE_DIRS, SDL2_LIBRARIES e SDL2_LIBRARY_FILE_LOCATION sono tipicamente definite in un file di configurazione fornito dalla libreria
		# SDL2 stessa (ad esempio SDL2Config.cmake). In questo caso, il comando find_package(SDL2) importa le variabili necessarie direttamente da quel file.
		set_target_properties(SDL2::SDL2 PROPERTIES
			INTERFACE_INCLUDE_DIRECTORIES ${SDL2_INCLUDE_DIRS}
			INTERFACE_LINK_LIBRARIES ${SDL2_LIBRARIES}
			IMPORTED_LOCATION ${SDL2_LIBRARY_FILE_LOCATION}
		)
	    # La variabile temporanea viene rimossa per evitare conflitti futuri.
		unset(SDL2_LIBRARY_FILE_LOCATION)
	endif()

	# minimum sdl version is 2.0.14 for controller LED support
	# Controlla che la versione minima di SDL2 sia 2.0.14 (necessaria per alcune funzionalità avanzate come il supporto ai LED del controller).
	# Se la versione soddisfa i requisiti, definisce il macro SDL2_FEATURE_CONTROLLER_LED, che può essere utilizzata nel codice C++ per abilitare funzionalità specifiche.
	if((${SDL_VERSION_MAJOR} GREATER_EQUAL 2) AND (${SDL2_VERSION_MINOR} GREATER_EQUAL 0) AND (${SDL2_VERSION_PATCH} GREATER_EQUAL 14))
		add_compile_definitions(SDL2_FEATURE_CONTROLLER_LED=1)
	endif()
endif()

# Aggiunge le sottodirectory LUA e Utility al progetto, permettendo di gestirle come parte della build principale.
add_subdirectory(LUA)
add_subdirectory(Utility)

# add_compile_definitions definisce una macro di precompiplazione per tutti i file sorgente del progetto.
add_compile_definitions(JPH_DEBUG_RENDERER=1)
set(PHYSICS_REPO_ROOT ${CMAKE_CURRENT_SOURCE_DIR})

# Include il file Jolt.cmake che si trova nella directory Jolt all'interno della directory PHYSICS_REPO_ROOT.
# In altre parole, questa riga esegue le istruzioni definite in Jolt.cmake, che contiene le istruzioni necessarie per configurare, 
# includere o costruire la libreria Jolt Physics.
include(${PHYSICS_REPO_ROOT}/Jolt/Jolt.cmake)
set_target_properties(Jolt PROPERTIES
	POSITION_INDEPENDENT_CODE ${WICKED_PIC}
)

# Utilizza file(GLOB ...) per raccogliere tutti i file .h e .cpp nella directory corrente.
# Rimuove esplicitamente offlineshadercompiler.cpp dalla lista dei file sorgente.
# CONFIGURE_DEPENDS indica a CMake di monitorare i file specificati nel pattern per eventuali modifiche durante la configurazione. 
# Se CMake rileva che i file sono cambiati (per esempio se un file è stato aggiunto, modificato o rimosso), eseguirà nuovamente la 
# configurazione quando esegui di nuovo il progetto (ad esempio, se esegui cmake .) senza bisogno di forzare una rigenerazione manuale del progetto.
# La fase di configurazione in CMake è il momento in cui il sistema di build (CMake) raccoglie informazioni sul progetto, 
# sul sistema in cui si sta costruendo il progetto e su come dovrebbe essere organizzato il processo di compilazione. 
# È il passo in cui CMake esamina l'ambiente di sviluppo (come le versioni dei compilatori, la presenza di librerie esterne, le opzioni di configurazione, ecc.) 
# e genera i file necessari per il sistema di build scelto (come Makefile, Visual Studio project files, ecc.).
#Fasi principali di CMake
# Configurazione (Configuration):
#	- Rilevamento delle variabili di sistema: CMake controlla il sistema per determinare caratteristiche come:
#   	- Il compilatore C++ e le sue versioni.
#		- La presenza di librerie di terze parti, come SDL2 o OpenGL.
#		- La configurazione dell'ambiente di esecuzione (ad esempio se stai compilando su Windows, Linux, macOS).
#	- Esecuzione dei CMakeLists.txt: CMake esegue i comandi nei file CMakeLists.txt per raccogliere informazioni su come dovrebbe essere costruito il progetto. 
#	  Questo include la raccolta di file di sorgente, la definizione dei target di compilazione, e la gestione delle dipendenze.
#	- Generazione dei file di build: In base alle informazioni raccolte e alle opzioni configurate, CMake genera i file di build (come Makefile, Visual Studio project files, ecc.).
# Compilazione (Build):
# 	Una volta completata la fase di configurazione, il sistema di build (ad esempio make o Visual Studio) è quindi utilizzato per compilare il codice sorgente in eseguibili, 
# 	librerie o altri artefatti.
#Installazione (Install):
# 	Dopo che il progetto è stato costruito, CMake può essere utilizzato per installare il software, ovvero copiare i file necessari nelle giuste directory del sistema. 
# 	Ad esempio, il binario, le librerie o gli header possono essere copiati in una directory di sistema come /usr/local/ su Linux, o C:/Program Files/ su Windows.
file(GLOB HEADER_FILES CONFIGURE_DEPENDS *.h)
file(GLOB SOURCE_FILES CONFIGURE_DEPENDS *.cpp)
list(REMOVE_ITEM SOURCE_FILES offlineshadercompiler.cpp)

# Crea una libreria (STATIC o SHARED, a seconda del tipo configurato) con i file sorgenti e header definiti.
add_library(${TARGET_NAME} ${WICKED_LIBRARY_TYPE}
	${SOURCE_FILES}
	${HEADER_FILES}
)
# MSVC è un identificatore in CMake che indica che il sistema di build sta usando il compilatore Microsoft Visual C++ (MSVC)
# CMAKE_CXX_FLAGS è una variabile che contiene le opzioni di compilazione per il compilatore C++. È usata da CMake per passare argomenti al compilatore.
if(MSVC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
endif()

# Crea un alias per il target di CMake già esistente, specificato dalla variabile ${TARGET_NAME}. 
# In pratica, WickedEngine diventa un altro nome per il target definito in ${TARGET_NAME}, che è WickedEngine_Windows.
add_library(WickedEngine ALIAS ${TARGET_NAME})
# set_target_properties è un comando di CMake che consente di impostare proprietà su un target.
# PUBLIC_HEADER è una proprietà speciale che viene utilizzata per definire quali file header dovrebbero essere considerati pubblici e 
# quindi inclusi nell'installazione del target. Questi file header vengono installati insieme alla libreria o all'eseguibile, e 
# quindi gli altri progetti che utilizzano questa libreria possono trovarli.
set_target_properties(${TARGET_NAME} PROPERTIES
	PUBLIC_HEADER "${HEADER_FILES}"
	POSITION_INDEPENDENT_CODE ${WICKED_PIC}
)

# target_include_directories è un comando di CMake utilizzato per specificare le directory dove il compilatore dovrebbe cercare i file di intestazione (.h) per il target ${TARGET_NAME}.
# SYSTEM indica che le directory specificate sono considerate come directory di sistema, il che significa che il compilatore non dovrebbe generare warning per gli include in queste directory.
# $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> indica che il compilatore dovrebbe cercare i file di intestazione nella directory corrente (CMAKE_CURRENT_SOURCE_DIR) durante la fase di build.
# $<BUILD_INTERFACE:...> è un generatore espressione di CMake. È un modo per specificare che una determinata directory di inclusione deve essere utilizzata solo durante la fase 
# di costruzione del progetto (build), non durante l'installazione.
# $<INSTALL_INTERFACE:...> è un'altra espressione di CMake che viene utilizzata per definire cosa accade durante la fase di installazione.
target_include_directories(${TARGET_NAME} SYSTEM PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/WickedEngine>
)

set(WICKEDENGINE_STATIC_LIBRARIES
	Jolt
	LUA
	Utility
)

if (WIN32)
	# target_compile_definitions definisce una macro di precompiplazione durante la compilazione del codice sorgente associato al target ${TARGET_NAME}.
	# PUBLIC specifica che la proprietà è applicata sia al target stesso sia a tutti i target che linkano a questo target.
	target_compile_definitions(${TARGET_NAME} PUBLIC
		UNICODE _UNICODE
	)

	# target_compile_options viene utilizzato per specificare le opzioni di compilazione (flag del compilatore) che devono essere applicate a uno specifico target in CMake.
	# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3") imposterebbe il flag /W3 per il compilatore C++ a livello globale.
	# PRIVATE specifica che l'opzione è applicata solo al target ${TARGET_NAME} e non ai target che dipendono da esso.
	target_compile_options(${TARGET_NAME} PRIVATE /W3)

	set(LIBDXCOMPILER "dxcompiler.dll")
else ()
	# `ska::flat_hash_map` has issues on linux because of the hash function being identity
	# in same cases. Use `robin_hood::unordered_flat_map` instead
	target_compile_definitions(${TARGET_NAME} PUBLIC WI_UNORDERED_MAP_TYPE=2)

	# target_link_libraries collega le librerie specificate (Threads::Threads, SDL2::SDL2, ecc.) al target ${TARGET_NAME}. 
	# Il qualificatore PUBLIC indica che queste librerie:
	#	- Saranno usate per compilare il target ${TARGET_NAME}.
	#	- Saranno propagate a tutti i target che dipendono da ${TARGET_NAME}.
	# $<$<BOOL:${OpenImageDenoise_FOUND}>:OpenImageDenoise> è una espressione generativa di CMake (generator expression), 
	# che include OpenImageDenoise solo se è stato trovato.
	# BOOL indica di controllare se la variabile OpenImageDenoise_FOUND è vera o falsa.
	# Se la condizione è vera, aggiunge OpenImageDenoise alla lista delle librerie collegate.
	# In pratica, se find_package(OpenImageDenoise ...) ha trovato la libreria, questa sarà collegata; altrimenti, sarà ignorata.
	target_link_libraries(${TARGET_NAME} PUBLIC
		Threads::Threads
		SDL2::SDL2
		$<$<BOOL:${OpenImageDenoise_FOUND}>:OpenImageDenoise> # links OpenImageDenoise only if it's found
	)

	# Imposta la proprietà POSITION_INDEPENDENT_CODE del target FAudio a ${WICKED_PIC}.
	set_target_properties(FAudio PROPERTIES
		POSITION_INDEPENDENT_CODE ${WICKED_PIC}
	)
	# Aggiunge la libreria FAudio alla lista delle librerie statiche di WickedEngine.
	set(WICKEDENGINE_STATIC_LIBRARIES ${WICKEDENGINE_STATIC_LIBRARIES} FAudio)

	if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
	target_compile_options(${TARGET_NAME} PRIVATE
			-Wuninitialized
			#-Wwrite-strings
			#-Winit-self
			#-Wreturn-type
			#-Wreorder
			#-Werror=delete-non-virtual-dtor
			#-Werror
			#uncomment this to stop the compilation at the first error
			# -Wfatal-errors
			)
	elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	# add some warnings and set them as errors
	# read more details here: https://gcc.gnu.org/onlinedocs/gcc/Warning-Options.html
	target_compile_options(${TARGET_NAME} PRIVATE
			-Wuninitialized
			-Wmaybe-uninitialized
			-Wwrite-strings
			-Winit-self
			-Wreturn-type
			-Wreorder
			-Werror=delete-non-virtual-dtor
			-Werror
			#uncomment this to stop the compilation at the first error
			# -Wfatal-errors
			)
	endif()

	target_link_libraries(${TARGET_NAME} PRIVATE dl)

	set(LIBDXCOMPILER "libdxcompiler.so")
endif()

target_link_libraries(${TARGET_NAME} PUBLIC ${WICKEDENGINE_STATIC_LIBRARIES})

if (PLATFORM MATCHES "SDL2")
	target_compile_definitions(${TARGET_NAME} PUBLIC SDL2=1)
endif()

# crea un target chiamato dxcompiler. Poiché è dichiarato come IMPORTED, il target rappresenta una libreria preesistente, 
# non una che verrà compilata dal progetto e niente verrà messo in build/WickedEngine/Debug.
add_library(dxcompiler SHARED IMPORTED)
# Questa riga imposta la posizione della libreria condivisa associata al target dxcompiler.
# TARGET dxcompiler: Si riferisce al target dxcompiler appena creato.
# PROPERTY IMPORTED_LOCATION: Specifica il percorso fisico del file .dll (o .so) che corrisponde al target. 
# Questo permette a CMake di sapere dove trovare il file quando è necessario collegarlo.
set_property(TARGET dxcompiler PROPERTY
		IMPORTED_LOCATION ${WICKED_ROOT_DIR}/WickedEngine/${LIBDXCOMPILER} )


# Copy the shader library next to the executable
# add_custom_command è un comando di CMake che consente di eseguire comandi personalizzati durante la fase di build.
# COMMAND ${CMAKE_COMMAND}: Specifica un comando da eseguire. CMAKE_COMMAND è un'espressione di CMake che rappresenta il comando cmake.
# POST_BUILD: Indica che il comando verrà eseguito dopo la fase di build del target specificato.
# -E copy_if_different: Specifica di copiare il file dxcompiler.dll (o libdxcompiler.so) nella directory CMAKE_CURRENT_BINARY_DIR, solo se è diverso (o non esiste già nella destinazione).
# CMAKE_CURRENT_BINARY_DIR è la directory di build corrente, dove CMake genera i file intermedi per il target attualmente configurato. In questo caso build\WickedEngine.
# (se ti serve la directory dove vengono generati gli eseguibili usa RUNTIME_OUTPUT_DIRECTORY).
# -E create_symlink: Crea un collegamento simbolico (symlink) da una directory sorgente a una destinazione.
# Su Windows, create_symlink implica che cmake userà l'API di sistema CreateSymbolicLink per creare un collegamento simbolico, il che richiede privilegi particolari
# che un normale utente potrebbe non avere. Per evitare di dover eseguire CMake come amministratore, si userà il comando mklink di Windows per creare il collegamento simbolico.
# Su Windows, create_symlink implica che cmake userà l'API di sistema CreateSymbolicLink per creare un collegamento simbolico, il che richiede privilegi particolari
# che un normale utente potrebbe non avere. Per evitare di dover eseguire CMake come amministratore, si userà il comando mklink di Windows per creare il collegamento simbolico.
if (WIN32)
    add_custom_command(
        TARGET ${TARGET_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${WICKED_ROOT_DIR}/WickedEngine/${LIBDXCOMPILER} ${CMAKE_CURRENT_BINARY_DIR}
		COMMAND if not exist "${CMAKE_CURRENT_BINARY_DIR}/shaders" \( mklink /J "${CMAKE_CURRENT_BINARY_DIR}/shaders" "${WICKED_ROOT_DIR}/WickedEngine/shaders" \)
    )
else()
    add_custom_command(
        TARGET ${TARGET_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${WICKED_ROOT_DIR}/WickedEngine/${LIBDXCOMPILER} ${CMAKE_CURRENT_BINARY_DIR}
        COMMAND ${CMAKE_COMMAND} -E create_symlink ${WICKED_ROOT_DIR}/WickedEngine/shaders ${CMAKE_CURRENT_BINARY_DIR}/shaders
    )
endif()


# OFFLINE SHADER COMPILER
# add_executable: Crea un nuovo target di tipo eseguibile (e lo mette in cartella dei file generati: build/WickedEngine/Debug).
# Il primo argomento (offlineshadercompiler) è il nome del target eseguibile.
# Il secondo argomento (offlineshadercompiler.cpp) è il file sorgente che verrà compilato per produrre l'eseguibile.
# Questo comando dice a CMake di generare un file eseguibile (offlineshadercompiler.exe su Windows o semplicemente offlineshadercompiler su sistemi Unix) 
# compilando il file specificato.
add_executable(offlineshadercompiler
		offlineshadercompiler.cpp
)

# collega la libreria WickedEngine_Windows (o WickedEngine_Linux) al target offlineshadercompiler
# così l'eseguibile offlineshadercompiler avrà accesso a tutte le funzionalità implementate in WickedEngine necessarie per compilare gli shader.
# PUBLIC significa che tutte le proprietà di compilazione (come include directories o definizioni) ereditate dalla libreria collegata saranno propagate 
# anche ai target che dipendono da offlineshadercompiler
target_link_libraries(offlineshadercompiler
		PUBLIC ${TARGET_NAME})


# comandi install eseguiti solo se l'utente esegue cmake --install build. 
# Altrimenti, nulla di ciò che è definito nelle istruzioni install() di questo file avverrà.
# In questo caso si vuole installare il target offlineshadercompiler nel sistema.
# TARGETS offlineshadercompiler: Specifica che solo il target offlineshadercompiler, 
# che è stato creato precedentemente come un eseguibile (con add_executable), deve essere installato.
# RUNTIME indica che solo i file runtime (eseguibili) devono essere installati.
install(TARGETS offlineshadercompiler
		RUNTIME DESTINATION "${CMAKE_INSTALL_LIBDIR}/WickedEngine")

# Installa (copia) la directory WickedEngine/Content in lib/WickedEngine (o usr/lib/WickedEngine in linux).
install(DIRECTORY "${WICKED_ROOT_DIR}/Content"
		DESTINATION "${CMAKE_INSTALL_LIBDIR}/WickedEngine")

#### LOCAL FOLDER INSTALLATION
# Create importable target here

#set_target_properties(${TARGET_NAME} PROPERTIES PUBLIC_HEADER "${HEADERS}")
# EXPORT_NAME è una proprietà di un target in CMake che definisce il nome con cui quel target sarà esportato nei file 
# di configurazione o esportazione generati da CMake. Questo nome può differire dal nome del target originale all'interno del progetto.
set_target_properties("${TARGET_NAME}" PROPERTIES EXPORT_NAME "${TARGET_NAME}")

# TARGETS "${TARGET_NAME}" ${WICKEDENGINE_STATIC_LIBRARIES} specifica che i target sono ${TARGET_NAME} e tutte le librerie statiche
# specificate in WICKEDENGINE_STATIC_LIBRARIES.
# EXPORT Engine-Targets: Aggiunge i target installati a un set di esportazione chiamato Engine-Targets.
# Utilizzato in comandi export per generare file di configurazione CMake, come WickedEngineTargets.cmake, per permettere l'uso di questi 
# target in progetti esterni con find_package (vedi comando export più avanti).
# LIBRARY: Specifica dove installare le librerie condivise (shared libraries).
# ARCHIVE: Specifica dove installare le librerie statiche (static libraries).
# PUBLIC_HEADER: Specifica dove installare i file header pubblici associati in precedenza al target usando la proprietà PUBLIC_HEADER.
install(TARGETS "${TARGET_NAME}" ${WICKEDENGINE_STATIC_LIBRARIES}
		EXPORT Engine-Targets
		LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
		ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}/WickedEngine"
		PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/WickedEngine")

# FILES: specifica i file da installare. In questo caso, il file dxcompiler.dll (o libdxcompiler.so)
install(FILES ${LIBDXCOMPILER}
		DESTINATION "${CMAKE_INSTALL_LIBDIR}/WickedEngine")

# export viene utilizzato per creare un file .cmake contenente informazioni su target esportabili. 
# Questo file può essere usato da altri progetti per configurare facilmente il proprio ambiente di build e includere le librerie esportate.
# EXPORT Engine-Targets: Specifica il gruppo di target da esportare.
# Questo gruppo di target deve essere stato già definito in una chiamata a install(TARGETS ...), come visto in precedenza.
# FILE: Specifica il percorso in cui salvare il file .cmake generato
# NAMESPACE: Specifica un namespace per i target esportati. 
# Ad esempio, il target Jolt diventa WickedEngine::Jolt.
# Questo evita conflitti di nomi, rendendo più chiaro a quale libreria appartengono i target.
# Il file WickedEngineTargets.cmake contiene:
#	Percorsi delle librerie e dei file associati ai target.
# 	Proprietà come flag di compilazione, directory di include, ecc.
#	Relazioni tra target (es., dipendenze tra di essi).
export(EXPORT Engine-Targets
		FILE "${CMAKE_BINARY_DIR}/cmake/WickedEngineTargets.cmake"
		NAMESPACE WickedEngine::
		)

# Simile al comando `export` sopra, ma crea il file durante la fase di installazione
# (eseguito quando l'utente utilizza il comando `cmake --install`).
# export(EXPORT ...) genera un file .cmake durante la fase di configurazione del progetto (quando esegui cmake ..), utile per il debug o per integrare il progetto in fase di sviluppo.
# install(EXPORT ...) genera il file .cmake solo durante l'installazione (cmake --install), rendendolo parte del pacchetto installato per gli utenti finali.
install(EXPORT Engine-Targets
		FILE WickedEngineTargets.cmake
		NAMESPACE WickedEngine::
		DESTINATION ${INSTALL_CONFIGDIR})

set(_CONFIG_INSTALL_DIR_INCLUDE "${WICKED_ROOT_DIR}")
set(_CONFIG_INSTALL_DIR_LIB "${WICKED_ROOT_DIR}")
# configure_package_config_file è un comando di CMake che prende un template di configurazione (un file .in) e lo elabora per generare un file di configurazione 
# del pacchetto (ad esempio, .cmake). Questo comando sostituisce le variabili nel template con i valori effettivi durante la fase di configurazione.
# ${CMAKE_CURRENT_SOURCE_DIR}/WickedEngineConfig.cmake.in: È il template di configurazione. È un file che contiene variabili come @_CONFIG_INSTALL_DIR_INCLUDE@, 
# che CMake sostituirà con i valori corrispondenti (impostati con set; vedi poco sopra).
# ${CMAKE_BINARY_DIR}/cmake/WickedEngineConfig.cmake: È il percorso dove il file configurato verrà salvato. Questo file verrà poi usato per la configurazione 
# durante la fase di installazione o durante l'inclusione del pacchetto.
# INSTALL_DESTINATION ${INSTALL_LOCAL_CONFIGDIR} a CMake un'indicazione su dove installare il file generato, nel contesto di un'installazione automatica 
# (cioè, quando si usa cmake --install per copiare i file in una destinazione finale). La vera installazione del file avviene con install(FILES ...), 
# che copia fisicamente il file dalla cartella di build alla destinazione finale.
configure_package_config_file(
		${CMAKE_CURRENT_SOURCE_DIR}/WickedEngineConfig.cmake.in
		${CMAKE_BINARY_DIR}/cmake/WickedEngineConfig.cmake
		INSTALL_DESTINATION ${INSTALL_LOCAL_CONFIGDIR}
)

# message(STATUS "STATUS: CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
# message(STATUS "STATUS: CMAKE_INSTALL_INCLUDEDIR: ${CMAKE_INSTALL_INCLUDEDIR}")
# message(STATUS "STATUS: CMAKE_INSTALL_LIBDIR: ${CMAKE_INSTALL_LIBDIR}")
# CMAKE_INSTALL_PREFIX: "C:/Program Files (x86)/WickedEngine" su Windows, "/usr/local/WickedEngine" su Linux.
# CMAKE_INSTALL_INCLUDEDIR: include
# CMAKE_INSTALL_LIBDIR: lib
set(_CONFIG_INSTALL_DIR_INCLUDE "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}")
set(_CONFIG_INSTALL_DIR_LIB "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
configure_package_config_file(
		${CMAKE_CURRENT_SOURCE_DIR}/WickedEngineConfig.cmake.in
		${CMAKE_BINARY_DIR}/cmake/install/WickedEngineConfig.cmake
		INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)

# copia il file WickedEngineConfig.cmake generato nella cartella di installazione finale.
install(FILES ${CMAKE_BINARY_DIR}/cmake/install/WickedEngineConfig.cmake
		DESTINATION ${INSTALL_CONFIGDIR}
)

# copia directory degli Shaders, escludendo le sottocartelle spirv e hlsl6 ed i file *.vcxitems*
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/shaders
		DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/WickedEngine"
		PATTERN "spirv" EXCLUDE
		PATTERN "hlsl6" EXCLUDE
		PATTERN "*.vcxitems*" EXCLUDE
		)

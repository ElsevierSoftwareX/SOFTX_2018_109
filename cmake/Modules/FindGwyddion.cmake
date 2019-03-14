# - Try to find Gwyddion libs
# Once done this will define
#  GWY_FOUND - System has Gwy libs
#  GWY_INCLUDE_DIRS - The Gwy libs include directories
#  GWY_LIBRARIES - The libraries needed to use Gwy libs
#  GWY_DEFINITIONS - Compiler switches required for using Gwy libs

# Steps:
# 1) try to use user-provided paths
# 2) try pkg-config
# 3) download and compile


## components: nogui -> gui = all
#if(NOT GWY_FIND_COMPONENTS)
#    set(GWY_FIND_COMPONENTS all)
#endif()

macro(check_gwy_dirs INCDIR LIBDIR)
  
  message("Looking for Gwyddion includes in ${INCDIR} ...")
 
  # gwyddion.h chosen as a representative of Gwy headers presence
  find_path(CHECKED_GWY_INCLUDE_DIR libgwyddion/gwyddion.h
    HINTS ${INCDIR}
    PATH_SUFFIXES gwyddion)
  if(NOT (CHECKED_GWY_INCLUDE_DIR STREQUAL "CHECKED_GWY_INCLUDE_DIR-NOTFOUND"))
    message("... includes found in ${CHECKED_GWY_INCLUDE_DIR}")
    set(GWY_INCLUDES_FOUND 1)
  else()
    message("... includes not found")
  endif()

  message("Looking for Gwyddion libraries in ${LIBDIR} ...")
  
  find_library(GWY_LIBRARY_GWYAPP gwyapp2 HINTS ${LIBDIR})
  find_library(GWY_LIBRARY_GWYDDION gwyddion2 HINTS ${LIBDIR})
  find_library(GWY_LIBRARY_GWYDGETS gwydgets2 HINTS ${LIBDIR})
  find_library(GWY_LIBRARY_GWYDRAW gwydraw2 HINTS ${LIBDIR})
  find_library(GWY_LIBRARY_MODULE gwymodule2 HINTS ${LIBDIR})
  find_library(GWY_LIBRARY_GWYPROCESS gwyprocess2 HINTS ${LIBDIR})

  # libgwyddion chosen as a representative of Gwy libraries presence
  if(NOT (GWY_LIBRARY_GWYDDION STREQUAL "GWY_LIBRARY_GWYDDION-NOTFOUND"))
    message("... libgwyddion in ${GWY_LIBRARY_GWYDDION}")
    set(GWY_LIBS_FOUND 1)
  else()
    message("... libgwyddion not found")
  endif()

  # in newer deb packaging, gwyconfig.h moved from lib to include
  find_path(CHECKED_GWYCONFIG_INCLUDE_DIR gwyconfig.h
    HINTS ${LIBDIR} ${INCDIR}
    PATH_SUFFIXES gwyddion gwyddion/include)

  if(NOT (CHECKED_GWYCONFIG_INCLUDE_DIR STREQUAL "CHECKED_GWYCONFIG_INCLUDE_DIR-NOTFOUND"))
    message("gwyconfig.h found in ${CHECKED_GWYCONFIG_INCLUDE_DIR}")
    set(GWY_GWYCONFIG_FOUND 1)
  else()
    message("gwyconfig.h not found")
  endif()
  
endmacro(check_gwy_dirs)

### 0) check for GTK2
find_package(PkgConfig REQUIRED)
if(WIN32 AND MINGW)
  # for some reason PKG_CONFIG_PATH is formatted as C:\\mingw32\\..., so change it to something pkg-config can use
  if(CMAKE_SIZEOF_VOID_P EQUAL 4) # 32-bit environment
    set(ENV{PKG_CONFIG_PATH} /mingw32/lib/pkgconfig;/mingw32/share/pkgconfig)
  elseif(CMAKE_SIZEOF_VOID_P EQUAL 8) # 64-bit environment
    set(ENV{PKG_CONFIG_PATH} /mingw64/lib/pkgconfig;/mingw64/share/pkgconfig)
  else()
    message(WARNING "Could not determine the system bitness, configuration will probably fail")
  endif()
endif()

pkg_check_modules(GTK2 REQUIRED gtk+-2.0>=2.24)

### 1) user-provided paths
if(GWY_INCLUDE_DIR AND GWY_LIBDIR)
  
  message(STATUS "Using paths to Gwyddion headers and libraries provided by the user ...")
  check_gwy_dirs(${GWY_INCLUDE_DIR} ${GWY_LIBDIR})

  if(GWY_INCLUDES_FOUND AND GWY_LIBS_FOUND AND GWY_GWYCONFIG_FOUND)
    message(STATUS "... success!")

    if(GTK2_FOUND)
      # Gwyddion might have been built with GtkGLExt
      # try to locate it (if we find it, then Gwyddion should have found and used it as well)
      pkg_check_modules(GTKGLEXT gtkglext-1.0)
    endif(GTK2_FOUND)
    
    set(GWY_INCLUDE_DIRS ${CHECKED_GWY_INCLUDE_DIR} ${CHECKED_GWYCONFIG_INCLUDE_DIR} ${GTKGLEXT_INCLUDE_DIRS})
    set(USER_GWY_FOUND 1)

    # try to make version detection and extraction cleaner
    file (STRINGS ${CHECKED_GWY_INCLUDE_DIR}/libgwyddion/gwyversion.h VERSION_LINE REGEX "#define GWY_VERSION_STRING \".*\"")
    string(REPLACE "#define GWY_VERSION_STRING \"" "" VERSION_LINE ${VERSION_LINE})
    string(REPLACE "\"" "" GWY_VERSION ${VERSION_LINE})
  else()
    message(STATUS "... Gwyddion headers and libraries provided by the user not found.")
  endif()
  
endif(GWY_INCLUDE_DIR AND GWY_LIBDIR)


### 2) if user-provided paths failed, try pkg-config
if(NOT USER_GWY_FOUND)
  
  message(STATUS "Trying to find Gwyddion headers and libraries using pkg-config ...")
  pkg_check_modules(PC_GWY gwyddion)

  if(PC_GWY_FOUND)
    message(STATUS "... success!")
    set(GWY_DEFINITIONS ${PC_GWY_CFLAGS_OTHER})

    # only Gwy libraries are provided
    find_library(GWY_LIBRARY_GWYDDION gwyddion2 HINTS ${PC_GWY_LIBDIR} ${PC_GWY_LIBRARY_DIRS})
    find_library(GWY_LIBRARY_GWYPROCESS gwyprocess2 HINTS ${PC_GWY_LIBDIR} ${PC_GWY_LIBRARY_DIRS})

    if(GTK2_FOUND)
      find_library(GWY_LIBRARY_GWYAPP gwyapp2 HINTS ${PC_GWY_LIBDIR} ${PC_GWY_LIBRARY_DIRS})
      find_library(GWY_LIBRARY_GWYDGETS gwydgets2 HINTS ${PC_GWY_LIBDIR} ${PC_GWY_LIBRARY_DIRS})
      find_library(GWY_LIBRARY_GWYDRAW gwydraw2 HINTS ${PC_GWY_LIBDIR} ${PC_GWY_LIBRARY_DIRS})
      find_library(GWY_LIBRARY_MODULE gwymodule2 HINTS ${PC_GWY_LIBDIR} ${PC_GWY_LIBRARY_DIRS})
    endif(GTK2_FOUND)

    # at the moment use include dirs as provided by pkg-config
    set(GWY_INCLUDE_DIRS ${PC_GWY_INCLUDE_DIR} ${PC_GWY_INCLUDE_DIRS})
    set(GWY_VERSION ${PC_GWY_VERSION})
    
  else(PC_GWY_FOUND)
    message(STATUS "... not found.")
  endif(PC_GWY_FOUND)
endif(NOT USER_GWY_FOUND)


### 3) download and build from source (check for existing custom build first)
if(NOT PC_GWY_FOUND)
  
  message(STATUS "Trying to use our own Gwyddion (downloaded and built from source) ...")

  # workaround for msys2
  # use CMAKE_BINDIR in execute_process as command parameters; CMAKE_BINARY_DIR when provided by cmake
  if(WIN32 AND MINGW)
    string(REPLACE "C:/msys64" "" CMAKE_BINDIR ${CMAKE_BINARY_DIR})
  else()
    set(CMAKE_BINDIR ${CMAKE_BINARY_DIR})
  endif()
  
  set(OWN_GWY_VER 2.52)
  set(OWN_GWY_SRC_FILE "gwyddion-${OWN_GWY_VER}.tar.gz")
  set(OWN_GWY_URL http://gwyddion.net/download/${OWN_GWY_VER}/${OWN_GWY_SRC_FILE})
  set(OWN_GWY_SRC_DIR "${CMAKE_BINARY_DIR}/gwyddion-${OWN_GWY_VER}")
  set(OWN_GWY_INST_DIR "${CMAKE_BINARY_DIR}/gwyddion-${OWN_GWY_VER}-inst")
  set(OWN_GWY_INCLUDE_DIR ${OWN_GWY_INST_DIR}/include/gwyddion)
  set(OWN_GWY_LIBDIR ${OWN_GWY_INST_DIR}/lib)

  # first check if an older build is present
  check_gwy_dirs(${OWN_GWY_INCLUDE_DIR} ${OWN_GWY_LIBDIR})

  if(NOT (GWY_INCLUDES_FOUND AND GWY_LIBS_FOUND AND GWY_GWYCONFIG_FOUND))

    if(NOT EXISTS ${CMAKE_BINARY_DIR}/${OWN_GWY_SRC_FILE})
      message(STATUS "Downloading Gwyddion sources from ${OWN_GWY_URL} ...")
      file(DOWNLOAD ${OWN_GWY_URL} ${CMAKE_BINARY_DIR}/${OWN_GWY_SRC_FILE} STATUS DLSTATUS)
      message(STATUS "Download finished with status: ${DLSTATUS}")
    endif()

    if(EXISTS ${CMAKE_BINARY_DIR}/${OWN_GWY_SRC_FILE}) # download succeeded or already present
      message(STATUS "Extracting Gwyddion sources from ${CMAKE_BINDIR}/${OWN_GWY_SRC_FILE} ...")

      # --force-local applies on Windows (MSVC generator) when the path starts with C:/...
      execute_process(COMMAND tar -x --force-local -f ${CMAKE_BINDIR}/${OWN_GWY_SRC_FILE} WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
      
      message(STATUS "Patching Gwyddion to compile just the libraries ...")
      # --forward: do not complain when the patch has already been applied
      if(GTK2_FOUND)
	execute_process(COMMAND patch --forward -i ${CMAKE_SOURCE_DIR}/gwyddion/Makefile.in.patch WORKING_DIRECTORY ${OWN_GWY_SRC_DIR})
      elseif()
	execute_process(COMMAND patch --forward -i ${CMAKE_SOURCE_DIR}/gwyddion/Makefile.in.nogui.patch WORKING_DIRECTORY ${OWN_GWY_SRC_DIR})
      endif()

      message(STATUS "Configuring Gwyddion (this will take some time) ...")
      execute_process(COMMAND sh ./configure --prefix=${OWN_GWY_INST_DIR}
	--disable-gtk-doc --disable-dependency-tracking 
	--disable-schemas-install --disable-desktop-file-update --disable-pygwy
	--without-pascal --without-perl --without-python --without-ruby
	--without-kde4-thumbnailer --without-unique
	# --without-libiconv-prefix --without-libintl-prefix
	--without-x --without-gtksourceview --without-gl --without-zlib --without-bzip2 --without-libxml2
	WORKING_DIRECTORY ${OWN_GWY_SRC_DIR}
	RESULT_VARIABLE RES
	OUTPUT_FILE ${CMAKE_BINARY_DIR}/gwyddion.configure.out.log
	ERROR_FILE ${CMAKE_BINARY_DIR}/gwyddion.configure.err.log)
      message(STATUS "Configure result: ${RES}")

      # try to determine the number of CPUs available to speed up Gwy compilation
      include(ProcessorCount)
      ProcessorCount(NCPUS)
      if(NCPUS EQUAL 0)
        set(NCPUS 1)
      endif()

      message(STATUS "Compiling Gwyddion libraries, using ${NCPUS} CPUs (this will take even more time) ...")
      execute_process(COMMAND make -j${NCPUS} CFLAGS=-Ofast
	WORKING_DIRECTORY ${OWN_GWY_SRC_DIR}
	RESULT_VARIABLE RES
	OUTPUT_FILE ${CMAKE_BINARY_DIR}/gwyddion.make.out.log
	ERROR_FILE ${CMAKE_BINARY_DIR}/gwyddion.make.err.log)
      message(STATUS "Compile result: ${RES}")
      
      message(STATUS "Installing Gwyddion libraries ...")
      execute_process(COMMAND make install
	WORKING_DIRECTORY ${OWN_GWY_SRC_DIR}
	RESULT_VARIABLE RES
	OUTPUT_FILE ${CMAKE_BINARY_DIR}/gwyddion.install.out.log
	ERROR_FILE ${CMAKE_BINARY_DIR}/gwyddion.install.err.log)
      message(STATUS "Install result: ${RES}")

      check_gwy_dirs(${OWN_GWY_INCLUDE_DIR} ${OWN_GWY_LIBDIR})

    else(EXISTS ${CMAKE_BINARY_DIR}/${OWN_GWY_SRC_FILE})
      message(SEND_ERROR "Gwyddion source package ${CMAKE_BINARY_DIR}/${GWY_SRC_FILE} not found, cannot continue")
    endif(EXISTS ${CMAKE_BINARY_DIR}/${OWN_GWY_SRC_FILE})

  endif(NOT (GWY_INCLUDES_FOUND AND GWY_LIBS_FOUND AND GWY_GWYCONFIG_FOUND))

  if(GWY_INCLUDES_FOUND AND GWY_LIBS_FOUND AND GWY_GWYCONFIG_FOUND)
    message(STATUS "... all the required files are present, success!")
    
    set(GWY_INCLUDE_DIRS ${CHECKED_GWY_INCLUDE_DIR} ${CHECKED_GWYCONFIG_INCLUDE_DIR})
    set(GWY_VERSION ${OWN_GWY_VER})
    set(OWN_GWY_FOUND 1)
  else()
    message(SEND_ERROR "... something went wrong, some of the required files are missing!")
  endif()

endif(NOT PC_GWY_FOUND)


### 4) finish
if(USER_GWY_FOUND OR PC_GWY_FOUND OR OWN_GWY_FOUND)
  mark_as_advanced(GWY_INCLUDE_DIR GWYCONFIG_INCLUDE_DIR
    GWY_LIBRARY_GWYAPP GWY_LIBRARY_GWYDDION GWY_LIBRARY_GWYDGETS
    GWY_LIBRARY_GWYDRAW GWY_LIBRARY_MODULE GWY_LIBRARY_GWYPROCESS)

  set(GWY_LIBRARIES
    ${GWY_LIBRARY_GWYAPP} ${GWY_LIBRARY_GWYDDION} ${GWY_LIBRARY_GWYDGETS}
    ${GWY_LIBRARY_GWYDRAW} ${GWY_LIBRARY_MODULE} ${GWY_LIBRARY_GWYPROCESS}
    ${GTKGLEXT_LIBRARIES})

  # handle the QUIETLY and REQUIRED arguments and set GWY_FOUND to TRUE if all listed variables are TRUE
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(GWY
    REQUIRED_VARS GWY_INCLUDE_DIRS GWY_LIBRARY_GWYDDION
    VERSION_VAR PC_GWY_VERSION
    HANDLE_COMPONENTS)
endif(USER_GWY_FOUND OR PC_GWY_FOUND OR OWN_GWY_FOUND)

# https://stackoverflow.com/questions/20746936/cmake-of-what-use-is-find-package-if-you-need-to-specify-cmake-module-path-an

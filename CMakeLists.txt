# zjistit spravnou min. verzi
cmake_minimum_required(VERSION 2.8)

message("Niget CMake configuration script")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

project(Niget LANGUAGES C)

include(CheckLanguage)
check_language(Fortran)
if(CMAKE_Fortran_COMPILER)
  enable_language(Fortran)
  set(FORTRAN ON)
  get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
  if(Fortran_COMPILER_NAME MATCHES "ifort.*")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -nofor_main")
  endif()
else()
  message(STATUS "No Fortran support, ODR functionality will not be available")
  add_definitions(-DNOFORTRAN)
  set(FORTRAN OFF)
endif()

# poresit, aby se tyhle/nektere promenne necachovaly (takze pri naslednem cmake .. se pouzily vychozi hodnoty)
option(DEBUG "Debug build" ON)
option(VERBOSE "Increase Niget's verbosity for debugging purposes" OFF)

if(DEBUG)
  set(CMAKE_BUILD_TYPE Debug)
  message(STATUS "Build type: DEBUG")
else()
  set(CMAKE_BUILD_TYPE Release)
  message(STATUS "Build type: RELEASE")
endif()

if(VERBOSE)
  add_definitions(-DSAVEMC)
  add_definitions(-DSAVEMC2)
endif()

if(MSVC)
  add_definitions(-D_USE_MATH_DEFINES)
endif()


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


### 1) look for GLib & GObject, REQUIRED
# napsat vlastni hledac GLibu, nebo neco, co by fungovalo i na Windows/MSVC?
# do not look for GLib explicitly - Gwyddion requires GObject, which brings in GLib as a dependency
# pkg_check_modules(GLIB REQUIRED glib-2.0)
pkg_check_modules(GOBJECT REQUIRED gobject-2.0)


### 2) look for GTK2
pkg_check_modules(GTK2 gtk+-2.0>=2.24)


### 3) look for corresponding Gwyddion components, REQUIRED
# Gwy components not implemented yet in FindGwyddion.cmake !
if(GTK2_FOUND)
  find_package(Gwyddion REQUIRED COMPONENTS all)
else()
  find_package(Gwyddion REQUIRED COMPONENTS nogui)
endif()


### 4) compile Fortran code if desired
if(FORTRAN)
  set(SRCS_ODRPACK real_precision.f lpkbls.f odr.f)
  add_library(odrpack OBJECT ${SRCS_ODRPACK})

  set(SRCS_ODRPACK_WRAPPERS odrwrap.f90 odrwrap_power_3param.f90 odrwrap_power_2param.f90 odrwrap_linear.f90 odrwrap_hertz_radial.f90)
  add_library(odrpack_wrappers OBJECT ${SRCS_ODRPACK_WRAPPERS})
  add_dependencies(odrpack_wrappers odrpack)
endif()


### 5) nonGUI parts
set(SRCS_COMMON_NOGUI fddata.c file-utils.c fit-utils.c niget-common.c instrument.c mc-utils.c settings.c)

# getline is not part of the standard C library on Windows; use a standalone implementation
if(WIN32)
  # a deprecated autotools symbol, prevents a warning on conflicting types for malloc and realloc
  add_definitions(-DSTDC_HEADERS)
  list(APPEND SRCS_COMMON_NOGUI getline.c)
endif()

add_library(common_nogui OBJECT ${SRCS_COMMON_NOGUI})
target_include_directories(common_nogui PUBLIC ${GOBJECT_INCLUDE_DIRS} ${GWY_INCLUDE_DIRS})

set(SRCS_TOOLS_NOGUI
  tool_op.c tool_op_unc.c
  tool_op_odr.c tool_op_odr_unc.c
  tool_tangent.c tool_tangent_unc.c
  tool_hertz.c tool_hertz_unc.c
  tool_hertz_odr.c tool_hertz_odr_unc.c
  tool_twoslopes.c tool_twoslopes_unc.c
  tool_stiffness.c tool_stiffness_unc.c
  tool_Ph2.c tool_apopins.c
  tool_epwork.c tool_epwork_unc.c)

add_library(tools_nogui OBJECT ${SRCS_TOOLS_NOGUI})
target_include_directories(tools_nogui PUBLIC ${GOBJECT_INCLUDE_DIRS} ${GWY_INCLUDE_DIRS})

set(SRCS_NIGET_CMD niget-cmd.c $<TARGET_OBJECTS:common_nogui> $<TARGET_OBJECTS:tools_nogui>)
if(FORTRAN)
  list(APPEND SRCS_NIGET_CMD $<TARGET_OBJECTS:odrpack> $<TARGET_OBJECTS:odrpack_wrappers>)
endif()

add_executable(niget-cmd ${SRCS_NIGET_CMD})
target_compile_options(niget-cmd PUBLIC -Wall)
target_include_directories(niget-cmd PUBLIC ${GOBJECT_INCLUDE_DIRS} ${GWY_INCLUDE_DIRS})
target_link_libraries(niget-cmd m ${GOBJECT_LIBRARIES}
  ${GWY_LIBRARY_GWYPROCESS} ${GWY_LIBRARY_GWYDDION})


### 6) GUI parts
if(GTK2_FOUND)
  set(SRCS_COMMON_GUI gui-utils.c)

  add_library(common_gui OBJECT ${SRCS_COMMON_GUI})
  target_include_directories(common_gui PUBLIC ${GTK2_INCLUDE_DIRS} ${GWY_INCLUDE_DIRS})

  set(SRCS_TOOLS_GUI
    tool_op_gui.c tool_op_unc_gui.c
    tool_op_odr_gui.c tool_op_odr_unc_gui.c
    tool_tangent_gui.c tool_tangent_unc_gui.c
    tool_hertz_gui.c tool_hertz_unc_gui.c
    tool_hertz_odr_gui.c tool_hertz_odr_unc_gui.c
    tool_twoslopes_gui.c tool_twoslopes_unc_gui.c
    tool_stiffness_gui.c tool_stiffness_unc_gui.c
    tool_Ph2_gui.c tool_apopins_gui.c
    tool_epwork_gui.c tool_epwork_unc_gui.c)
  add_library(tools_gui OBJECT ${SRCS_TOOLS_GUI})
  target_include_directories(tools_gui PUBLIC ${GTK2_INCLUDE_DIRS} ${GWY_INCLUDE_DIRS})
  
  set(SRCS_NIGET niget-gtk.c
    $<TARGET_OBJECTS:common_nogui> $<TARGET_OBJECTS:common_gui>
    $<TARGET_OBJECTS:tools_nogui> $<TARGET_OBJECTS:tools_gui>
    tool_contact.c dialog_about.c dialog_area_calibration.c)
  if(FORTRAN)
    list(APPEND SRCS_NIGET $<TARGET_OBJECTS:odrpack> $<TARGET_OBJECTS:odrpack_wrappers>)
  endif()
  add_executable(niget ${SRCS_NIGET})
  target_compile_options(niget PUBLIC -Wall)
  target_include_directories(niget PUBLIC ${GTK2_INCLUDE_DIRS} ${GWY_INCLUDE_DIRS})
  target_link_libraries(niget m ${GTK2_LIBRARIES}
    ${GWY_LIBRARY_GWYPROCESS} ${GWY_LIBRARY_GWYDDION} ${GWY_LIBRARY_GWYDRAW} ${GWY_LIBRARY_GWYDGETS} ${GWY_LIBRARY_GWYAPP} ${GWY_LIBRARY_MODULE})

else(GTK2_FOUND)
  message("GTK 2.24 or newer not found, will build only commandline version")
endif(GTK2_FOUND)

### 7) documentation

set(DOC_SRC_DIR ${CMAKE_SOURCE_DIR}/doc)

### 7a) PDF

find_package(LATEX COMPONENTS PDFLATEX BIBTEX)
if(LATEX_PDFLATEX_FOUND AND LATEX_BIBTEX_FOUND)
  message(STATUS "pdfLaTeX and BibTeX found, enabling PDF documentation build")
  add_custom_target(pdfdoc
    TEXINPUTS=${DOC_SRC_DIR}: pdflatex ${DOC_SRC_DIR}/niget.tex
    COMMAND BIBINPUTS=${DOC_SRC_DIR} openout_any=a bibtex ${CMAKE_BINARY_DIR}/niget # alternativne TMPDIR=.
    COMMAND TEXINPUTS=${DOC_SRC_DIR}: pdflatex ${DOC_SRC_DIR}/niget.tex
    COMMAND TEXINPUTS=${DOC_SRC_DIR}: pdflatex ${DOC_SRC_DIR}/niget.tex
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
else(LATEX_PDFLATEX_FOUND AND LATEX_BIBTEX_FOUND)
  message("pdfLaTeX or BibTeX not found, cannot build PDF documentation")
endif(LATEX_PDFLATEX_FOUND AND LATEX_BIBTEX_FOUND)

### 7b) HTML

find_program(LATEXML latexml)
if(LATEXML)
  message(STATUS "LaTeXML found, enabling HTML documentation build")
  add_custom_target(htmldoc
    latexml --destination=niget.xml  ${DOC_SRC_DIR}/niget.tex
    COMMAND latexmlpost --destination=html4/index.html --format=html4 --split niget.xml
    COMMAND latexmlpost --destination=html5/index.html --format=html5 --split niget.xml
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
else(LATEXML)
  message("STATUS LaTeXML not found, cannot build HTML documentation")
endif(LATEXML)


# https://cmake.org/Wiki/CMake:How_To_Find_Libraries

# zvazit neco jako configure_file - vygeneruje .h soubor s #defines

#cmake -G "Your Generator" -D CMAKE_C_COMPILER=gcc-4.2 -D CMAKE_CXX_COMPILER=g++-4.2 path/to/your/source

#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#  message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

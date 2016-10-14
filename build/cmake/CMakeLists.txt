cmake_minimum_required(VERSION 2.8.12)

project(appplugin)

cmake_policy(SET CMP0020 NEW)

if(${CMAKE_MAJOR_VERSION} GREATER 2)
  cmake_policy(SET CMP0007 NEW)
  cmake_policy(SET CMP0026 OLD)
  cmake_policy(SET CMP0042 NEW)  
  cmake_policy(SET CMP0043 NEW)  
  cmake_policy(SET CMP0053 OLD)  
endif()

set(APPPLUGIN_ABSOLUTE_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../..)

set_directory_properties(PROPERTIES CLEAN_NO_CUSTOM 1)

if(MSVC)
  set(CMAKE_DEBUG_POSTFIX "d")
endif()

#
# main variables for the build script
#

# variables that control the system versioning
set(APPPLUGIN_VERSION_MAJOR 1)
set(APPPLUGIN_VERSION_MINOR 0)
set(APPPLUGIN_VERSION_PATCH 0)
set(APPPLUGIN_VERSION_STATUS "rc")
set(APPPLUGIN_VERSION_STRING "${APPPLUGIN_VERSION_MAJOR}.${APPPLUGIN_VERSION_MINOR}.${APPPLUGIN_VERSION_PATCH}-${APPPLUGIN_VERSION_STATUS}")
set(APPPLUGIN_VERSION "0x050000")
set(APPPLUGIN_PKG_VERSION ${APPPLUGIN_VERSION_MAJOR}.${APPPLUGIN_VERSION_MINOR}.${APPPLUGIN_VERSION_PATCH})

find_package(terralib REQUIRED)

find_package(Boost REQUIRED)

find_package(Qt5 5.1 REQUIRED COMPONENTS Core Gui Widgets PrintSupport)

find_package(Qt5LinguistTools REQUIRED COMPONENTS lupdate lrelease)

add_definitions(-D_SCL_SECURE_NO_WARNINGS -DBOOST_ALL_NO_LIB -DBOOST_FILESYSTEM_VERSION=3 -DBOOST_UBLAS_TYPE_CHECK=0)

set(CMAKE_AUTOMOC ON)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

set(APPPLUGIN_BASE_DESTINATION_DIR "")

set (TE_3RDPARTY "$ENV{TE_3RDPARTY}" CACHE PATH "Where to find TerraLib 3rdParty?")

if(NOT DEFINED APPPLUGIN_DESTINATION_LIBRARY)
  set(APPPLUGIN_DESTINATION_LIBRARY "${APPPLUGIN_BASE_DESTINATION_DIR}lib")
endif()

if(NOT DEFINED APPPLUGIN_DESTINATION_RUNTIME)
  set(APPPLUGIN_DESTINATION_RUNTIME "${APPPLUGIN_BASE_DESTINATION_DIR}bin")
endif()

if(NOT DEFINED APPPLUGIN_DESTINATION_ARCHIVE)
  set(APPPLUGIN_DESTINATION_ARCHIVE "${APPPLUGIN_BASE_DESTINATION_DIR}lib")
endif()

if(NOT DEFINED APPPLUGIN_DESTINATION_HEADERS)
  set(APPPLUGIN_DESTINATION_HEADERS "${APPPLUGIN_BASE_DESTINATION_DIR}include")
endif()

if(NOT DEFINED APPPLUGIN_DESTINATION_SHARE)
  set(APPPLUGIN_DESTINATION_SHARE "${APPPLUGIN_BASE_DESTINATION_DIR}share")
endif()

add_subdirectory(appplugin)

TeInstallQt5Plugins()

if(WIN32)
  get_filename_component(_appShare ${APPPLUGIN_ABSOLUTE_ROOT_DIR}/share/terralib DIRECTORY)
  get_filename_component(_teShare ${terralib_INCLUDE_DIRS}/../share/terralib DIRECTORY) 
  get_filename_component(_qtPlugins ${terralib_INCLUDE_DIRS}/../qtplugins/platforms DIRECTORY) 

  install (
    FILES ${_teplugins} ${_dep}
    DESTINATION bin
    COMPONENT runtime
  )  
  
  install(
    DIRECTORY ${_appShare}
    DESTINATION .
    COMPONENT runtime
  )
  
  install(
    DIRECTORY ${_teShare}
    DESTINATION .
    COMPONENT runtime
  )
   
  file(GLOB dll_list ${terralib_DIR}/../../../bin/*.dll)
  file(GLOB exe_list ${terralib_DIR}/../../../bin/*.exe)
  
  set(_teBins ${dll_list} ${exe_list})
  
  install(
    FILES ${_teBins}
    DESTINATION bin
    COMPONENT runtime
  )
endif()


include(InstallRequiredSystemLibraries)

# CPack
set(CPACK_GENERATOR NSIS)
set(CPACK_PACKAGE_NAME "TerraView APP")
set(CPACK_PACKAGE_VENDOR "Terralib3rdParty")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "TerraView APP")
set(CPACK_PACKAGE_VERSION "1.0.0")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "TerraView APP")
SET(CPACK_NSIS_MODIFY_PATH OFF)

#for the sortcuts
set(shortCutName "TerraView APP ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}${bitName}") 
set(shortCutNameTerraview "TerraView 5.2.0")

set(shortcutAPP_startmenu "CreateShortCut \\\"$SMPROGRAMS\\\\$STARTMENU_FOLDER/${shortCutName}.lnk\\\" \\\"$INSTDIRb\\\\bin\\\\terraview.exe\\\" \\\"\\\" \\\"$INSTDIR\\\\share\\\\terralib\\\\images\\\\ico\\\\terralib-globe.ico\\\"")
set(shortcutAPP_desktop "CreateShortCut \\\"$DESKTOP\\\\${shortCutName}.lnk\\\" \\\"$INSTDIR\\\\bin\\\\terraview.exe\\\" \\\"\\\" \\\"$INSTDIR\\\\share\\\\terralib\\\\images\\\\ico\\\\terralib-globe.ico\\\"")
set(shortcutTV_startmenu "CreateShortCut \\\"$SMPROGRAMS\\\\$STARTMENU_FOLDER/${shortCutNameTerraview}.lnk\\\" \\\"$INSTDIRb\\\\bin\\\\terraview.exe\\\" \\\"\\\" \\\"$INSTDIR\\\\share\\\\terralib\\\\images\\\\ico\\\\terralib-globe.ico\\\"")
set(shortcutTV_desktop "CreateShortCut \\\"$DESKTOP\\\\${shortCutNameTerraview}.lnk\\\" \\\"$INSTDIR\\\\bin\\\\terraview.exe\\\" \\\"\\\" \\\"$INSTDIR\\\\share\\\\terralib\\\\images\\\\ico\\\\terralib-globe.ico\\\"")
set(full_shortcut_list "${shortcutAPP_startmenu};${shortcutAPP_desktop};${shortcutTV_startmenu};${shortcutTV_desktop}")

set(CPACK_PACKAGE_EXECUTABLES "terraview;${shortCutName}" "terraview;${shortCutNameTerraview}")
set(CPACK_NSIS_CREATE_ICONS_EXTRA "SetOutPath \\\"$INSTDIR\\\\bin\\\";${full_shortcut_list}")

list (APPEND CPACK_NSIS_DELETE_ICONS_EXTRA "Delete \\\"$SMPROGRAMS\\\\$START_MENU\\\\${shortCutName}.lnk\\\";Delete \\\"$DESKTOP\\\\${shortCutName}.lnk\\\"")
list (APPEND CPACK_NSIS_DELETE_ICONS_EXTRA ";Delete \\\"$SMPROGRAMS\\\\$START_MENU\\\\${shortCutNameTerraview}.lnk\\\";Delete \\\"$DESKTOP\\\\${shortCutNameTerraview}.lnk\\\"")

INCLUDE(CPack)


if(WIN32)
  add_definitions(-DAPPPLUGINDLL)
endif()

set(APPPLG_ABSOLUTE_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../..)

include_directories(
  ${APPPLG_ABSOLUTE_ROOT_DIR}/src
  ${terralib_INCLUDE_DIRS}
  ${terralib_DIR}
  ${Boost_INCLUDE_DIR}
)

file(GLOB APPPLG_SRC_FILES ${APPPLG_ABSOLUTE_ROOT_DIR}/src/*.cpp)
file(GLOB APPPLG_HDR_FILES ${APPPLG_ABSOLUTE_ROOT_DIR}/src/*.h)

source_group("Source Files"  FILES ${APPPLG_SRC_FILES})
source_group("Header Files"  FILES ${APPPLG_HDR_FILES})

#core
file(GLOB APPPLG_CORE_SRC_FILES ${APPPLG_ABSOLUTE_ROOT_DIR}/src/core/*.cpp)
file(GLOB APPPLG_CORE_HDR_FILES ${APPPLG_ABSOLUTE_ROOT_DIR}/src/core/*.h)

source_group("Source Files\\core"  FILES ${APPPLG_CORE_SRC_FILES})
source_group("Header Files\\core"  FILES ${APPPLG_CORE_HDR_FILES})

#qt
file(GLOB APPPLG_QT_SRC_FILES ${APPPLG_ABSOLUTE_ROOT_DIR}/src/qt/*.cpp)
file(GLOB APPPLG_QT_HDR_FILES ${APPPLG_ABSOLUTE_ROOT_DIR}/src/qt/*.h)

source_group("Source Files\\qt"  FILES ${APPPLG_QT_SRC_FILES})
source_group("Header Files\\qt"  FILES ${APPPLG_QT_HDR_FILES})

file(GLOB APPPLG_QT_UI_FILES  ${APPPLG_ABSOLUTE_ROOT_DIR}/src/qt/ui/*.ui)

source_group("UI Files"  FILES ${APPPLG_QT_UI_FILES})

# uic'ing
QT5_WRAP_UI(APPPLG_GEN_HDR_FILES ${APPPLG_QT_UI_FILES})

set(APPPLG_FILES ${APPPLG_SRC_FILES} 
				   ${APPPLG_HDR_FILES} 
				   ${APPPLG_CORE_SRC_FILES} 
				   ${APPPLG_CORE_HDR_FILES} 
				   ${APPPLG_QT_SRC_FILES} 
				   ${APPPLG_QT_HDR_FILES}
				   ${APPPLG_GEN_HDR_FILES})
				   
set(APPPLG_HDR_FILES ${APPPLG_HDR_FILES} 
				     ${APPPLG_CORE_HDR_FILES} 
				     ${APPPLG_QT_HDR_FILES}
				     ${APPPLG_GEN_HDR_FILES})
				   
add_library(app_plugin SHARED ${APPPLG_FILES})

target_link_libraries(app_plugin terralib_mod_plugin terralib_mod_qt_apf)

set_target_properties(app_plugin
                      PROPERTIES VERSION ${APPPLUGIN_VERSION_MAJOR}.${APPPLUGIN_VERSION_MINOR}
                                 SOVERSION ${APPPLUGIN_VERSION_MAJOR}.${APPPLUGIN_VERSION_MINOR})

qt5_use_modules(app_plugin Widgets)

install(TARGETS app_plugin
        EXPORT app-targets
        RUNTIME DESTINATION ${APPPLUGIN_DESTINATION_RUNTIME} COMPONENT runtime
        LIBRARY DESTINATION ${APPPLUGIN_DESTINATION_LIBRARY} COMPONENT runtime
        ARCHIVE DESTINATION ${APPPLUGIN_DESTINATION_ARCHIVE} COMPONENT runtime)

install(FILES ${APPPLG_HDR_FILES}
        DESTINATION ${APPPLUGIN_DESTINATION_HEADERS}/app_plugin COMPONENT devel)

export(TARGETS app_plugin APPEND FILE ${CMAKE_BINARY_DIR}/app-exports.cmake)

#add_custom_command(
#    TARGET app_plugin
#    POST_BUILD
#	COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_BINARY_DIR}/bin/Debug/app_plugind.dll" "${terralib_DIR}/debug"
#    COMMENT "Copying plugin files..."
#    VERBATIM
#)

PROJECT(GoTopology)


# Include directories

INCLUDE_DIRECTORIES(
  ${GoTopology_SOURCE_DIR}/include
  ${GoToolsCore_SOURCE_DIR}/include
  ${GoTools_COMMON_INCLUDE_DIRS}
  )


# Linked in libraries

SET(DEPLIBS
  GoToolsCore
  sisl
  )


# Make the GoTopology library
FILE(GLOB_RECURSE GoTopology_SRCS src/*.C include/*.h)
if (BUILD_AS_SHARED_LIBRARY)
    ADD_LIBRARY(GoTopology SHARED ${GoTopology_SRCS})
else (BUILD_AS_SHARED_LIBRARY)
    ADD_LIBRARY(GoTopology ${GoTopology_SRCS})
endif (BUILD_AS_SHARED_LIBRARY)
TARGET_LINK_LIBRARIES(GoTopology ${DEPLIBS})
SET_PROPERTY(TARGET GoTopology
  PROPERTY FOLDER "GoTopology/Libs")
SET_TARGET_PROPERTIES(GoTopology PROPERTIES SOVERSION ${GoTools_ABI_VERSION})


# Apps and tests
MACRO(ADD_APPS SUBDIR PROPERTY_FOLDER IS_TEST)
  FILE(GLOB_RECURSE GoTopology_APPS ${SUBDIR}/*.C)
  FOREACH(app ${GoTopology_APPS})
    GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
    ADD_EXECUTABLE(${appname} ${app})
    TARGET_LINK_LIBRARIES(${appname} GoTopology ${DEPLIBS})
    SET_TARGET_PROPERTIES(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${SUBDIR})
    SET_PROPERTY(TARGET ${appname}
      PROPERTY FOLDER "GoTopology/${PROPERTY_FOLDER}")
    IF(${IS_TEST})
      ADD_TEST(${appname} ${SUBDIR}/${appname}
		--log_format=XML --log_level=all --log_sink=../Testing/${appname}.xml)
      SET_TESTS_PROPERTIES( ${appname} PROPERTIES LABELS "${SUBDIR}" )
    ENDIF(${IS_TEST})
  ENDFOREACH(app)
ENDMACRO(ADD_APPS)

IF(GoTools_COMPILE_APPS)
  ADD_APPS(app "Apps" FALSE)
  ADD_APPS(examples "Examples" FALSE)
ENDIF(GoTools_COMPILE_APPS)

IF(GoTools_COMPILE_TESTS)
  SET(DEPLIBS ${DEPLIBS} ${Boost_LIBRARIES})
  ADD_APPS(test/unit "Unit Tests" TRUE)
  ADD_APPS(test/integration "Integration Tests" TRUE)
  ADD_APPS(test/acceptance "Acceptance Tests" TRUE)
ENDIF(GoTools_COMPILE_TESTS)


# 'install' target

IF(WIN32)
  # Windows
  INSTALL(TARGETS GoTopology DESTINATION ${GoTools_INSTALL_PREFIX}/lib)
  INSTALL(DIRECTORY include/GoTools/topology
    DESTINATION ${GoTools_INSTALL_PREFIX}/include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ELSE(WIN32)
  # Linux
  INSTALL(TARGETS GoTopology DESTINATION lib COMPONENT topology)
  INSTALL(DIRECTORY include/GoTools/topology
    DESTINATION include/GoTools
    COMPONENT topology-dev
    FILES_MATCHING PATTERN "*.h"
    PATTERN "*.hpp"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)

SET(CPACK_STRIP_FILES ${CPACK_STRIP_FILES} libGoTopology.so)

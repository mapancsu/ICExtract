include(GenerateExportHeader)
include_directories(${PYTHON_INCLUDE_DIRS})

add_library(ICExtract SHARED PIC.cpp)
target_link_libraries(ICExtract ${PYTHON_LIBRARIES} OpenMS)
set_target_properties(ICExtract PROPERTIES COMPILE_FLAGS "-DICEXTRACT_EXPORTS")
install(TARGETS ICExtract DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../python)
generate_export_header(ICExtract)

add_executable(testICExtract main.cpp PIC.cpp)
set_target_properties(testICExtract PROPERTIES COMPILE_FLAGS "-DICEXTRACT_STATIC_DEFINE")
target_link_libraries(testICExtract ${PYTHON_LIBRARIES} OpenMS)
include(GenerateExportHeader)
include_directories(${PYTHON_INCLUDE_DIRS})

add_library(decDIA SHARED PIC.cpp) #utils.cpp 
#target_include_directories(decDIA PRIVATE ${TBB_INCLUDE_DIRS})
target_link_libraries(decDIA ${PYTHON_LIBRARIES} OpenMS)
set_target_properties(decDIA PROPERTIES COMPILE_FLAGS "-DDECDIA_EXPORTS")
install(TARGETS decDIA DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../python)
generate_export_header(decDIA)

add_executable(testdecDIA main.cpp PIC.cpp)
set_target_properties(testdecDIA PROPERTIES COMPILE_FLAGS "-DDECDIA_STATIC_DEFINE")
#target_include_directories(testMZXML PRIVATE ${TBB_INCLUDE_DIRS})
target_link_libraries(testdecDIA ${PYTHON_LIBRARIES} OpenMS)
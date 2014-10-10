file(GLOB EMTG_FILES *.h *.cpp)


cmake_policy(SET CMP0003 NEW)

add_executable(emtg ${EMTG_FILES})
target_link_libraries(emtg ${EMTG_LIBRARIES})
install(TARGETS emtg RUNTIME DESTINATION "${PROJECT_SOURCE_DIR}/bin")


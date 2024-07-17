set(LEBEDEV_SOURCE_DIR ${PROJECT_SOURCE_DIR}/share/lebedev CACHE STRING "Path to SAD basis and density files")
set(LEBEDEV_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/share/${PROJECT_NAME}/lebedev)
install(
  DIRECTORY share/lebedev
  DESTINATION share/${PROJECT_NAME}
  )


set(PROGS
  sc
)

foreach(PROG ${PROGS})
  set(PROG_TARGET ${PROG})
  add_executable(${PROG_TARGET} ${PROG}.cpp)
  target_include_directories(${PROG_TARGET} PRIVATE ../include)
  target_link_libraries(${PROG_TARGET} PRIVATE
    ZLIB::ZLIB
    cxxopts::cxxopts
  )
  install(TARGETS ${PROG_TARGET} RUNTIME DESTINATION bin COMPONENT ${PROG_TARGET})
endforeach()
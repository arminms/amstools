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

macro(do_test prog arg result)
    add_test(test-${prog}${arg} ${prog} ${arg})
    set_tests_properties(test-${prog}${arg} PROPERTIES PASS_REGULAR_EXPRESSION ${result})
endmacro(do_test)
 
do_test(sc -p "Option ‘p’ does not exist")

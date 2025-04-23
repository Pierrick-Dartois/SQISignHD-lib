set(SOURCE_FILES_SQISIGNDIM2_GENERIC_REF_TESTS
    ${SQISIGNDIM2X_DIR}/test/test_sqisigndim2.c
)
add_executable(sqisign_test_sqisigndim2_${SVARIANT_LOWER} ${SOURCE_FILES_SQISIGNDIM2_GENERIC_REF_TESTS})
target_link_libraries(sqisign_test_sqisigndim2_${SVARIANT_LOWER} ${LIB_SQISIGNDIM2_${SVARIANT_UPPER}} ${LIB_DIM2ID2ISO_${SVARIANT_UPPER}} ${LIB_HD_${SVARIANT_UPPER}} ${LIB_ID2ISO_${SVARIANT_UPPER}} ${LIB_EC} ${LIB_KLPT_${SVARIANT_UPPER}} ${LIB_QUATERNION} ${LIB_PRECOMP_${SVARIANT_UPPER}} ${LIB_INTBIG} ${LIB_PUBLIC} ${LIB_GF_${SVARIANT_UPPER}} ${LIB_EC_${SVARIANT_UPPER}} ${GMP} sqisign_common_sys )
target_include_directories(sqisign_test_sqisigndim2_${SVARIANT_LOWER} PRIVATE ${INC_PUBLIC} ${INC_INTBIG} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_QUATERNION} ${INC_KLPT} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_EC} ${INC_COMMON} ${INC_ID2ISO} ${INC_HD} ${INC_DIM2ID2ISO} ${INC_SQISIGNDIM2} ./include/)

add_test(sqisign_test_sqisigndim2_${SVARIANT_LOWER} sqisign_test_sqisigndim2_${SVARIANT_LOWER})

add_custom_command(
  TARGET sqisign_test_sqisigndim2_${SVARIANT_LOWER}
  POST_BUILD
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy $<TARGET_FILE:sqisign_test_sqisigndim2_${SVARIANT_LOWER}> ${CMAKE_BINARY_DIR}/test/sqisign2d_${SVARIANT_LOWER}
)

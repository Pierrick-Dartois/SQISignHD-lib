add_executable(sqisign_test_gf_${SVARIANT_LOWER}_fp ${GFX_DIR}/test/test_fp.c ${GFX_DIR}/test/test_utils.c ${GFX_DIR}/test/sha3.c)
target_link_libraries(sqisign_test_gf_${SVARIANT_LOWER}_fp ${LIB_GF_${SVARIANT_UPPER}})
target_include_directories(sqisign_test_gf_${SVARIANT_LOWER}_fp PRIVATE ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_COMMON} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_PUBLIC})
add_test(sqisign_test_gf_${SVARIANT_LOWER}_fp sqisign_test_gf_${SVARIANT_LOWER}_fp test ${SQISIGN_TEST_REPS})

add_executable(sqisign_test_gf_${SVARIANT_LOWER}_fp2 ${GFX_DIR}/test/test_fp2.c ${GFX_DIR}/test/test_utils.c ${GFX_DIR}/test/sha3.c)
target_link_libraries(sqisign_test_gf_${SVARIANT_LOWER}_fp2 ${LIB_GF_${SVARIANT_UPPER}})
target_include_directories(sqisign_test_gf_${SVARIANT_LOWER}_fp2 PRIVATE ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_COMMON} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_PUBLIC})
add_test(sqisign_test_gf_${SVARIANT_LOWER}_fp2 sqisign_test_gf_${SVARIANT_LOWER}_fp2 test ${SQISIGN_TEST_REPS})

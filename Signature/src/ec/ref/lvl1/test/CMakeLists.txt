add_executable(biextension.test_${SVARIANT_LOWER} ${ECX_DIR}/test/biextension-test.c)
	target_include_directories(biextension.test_${SVARIANT_LOWER} PUBLIC ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_INTBIG} ${INC_PRECOMP_${SVARIANT_UPPER}} ${PROJECT_SOURCE_DIR}/include  ../include ${INC_EC} ${INC_COMMON} ${INC_QUATERNION} .)
	target_link_libraries(biextension.test_${SVARIANT_LOWER} ${LIB_INTBIG} ${LIB_EC_${SVARIANT_UPPER}} ${LIB_QUATERNION} sqisign_common_sys)

add_executable(basis-gen.test_${SVARIANT_LOWER} ${ECX_DIR}/test/basis-gen-test.c)
	target_include_directories(basis-gen.test_${SVARIANT_LOWER} PUBLIC ${ECX_DIR}/test ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_INTBIG} ${INC_PRECOMP_${SVARIANT_UPPER}} ${PROJECT_SOURCE_DIR}/include  ../include ${INC_EC} ${INC_COMMON} .)
	target_link_libraries(basis-gen.test_${SVARIANT_LOWER} ${LIB_INTBIG} ${LIB_EC_${SVARIANT_UPPER}})

add_test(ec_biextension.test_${SVARIANT_LOWER} biextension.test_${SVARIANT_LOWER} ${SQISIGN_TEST_REPS})
add_test(ec_basis_gen.test_${SVARIANT_LOWER} basis-gen.test_${SVARIANT_LOWER} ${SQISIGN_TEST_REPS})


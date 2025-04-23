
set(SOURCE_FILES_GF_${SVARIANT_UPPER}_BROADWELL
    gf65376.c
    fp.c
    ${GFX_DIR}/fp2.c
    ${GFX_DIR}/mp.c
)

add_library(${LIB_GF_${SVARIANT_UPPER}} ${SOURCE_FILES_GF_${SVARIANT_UPPER}_BROADWELL})
target_include_directories(${LIB_GF_${SVARIANT_UPPER}} PRIVATE common ${INC_COMMON} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} include ${PROJECT_SOURCE_DIR}/include ${INC_COMMON})
target_compile_options(${LIB_GF_${SVARIANT_UPPER}} PRIVATE ${C_OPT_FLAGS})

add_subdirectory(test)

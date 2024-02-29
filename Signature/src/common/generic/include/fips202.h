// SPDX-License-Identifier: Apache-2.0

#ifndef FIPS202_H
#define FIPS202_H

#include <stddef.h>

int SHAKE128(unsigned char *output, size_t outputByteLen, const unsigned char *input, size_t inputByteLen);
int SHAKE256(unsigned char *output, size_t outputByteLen, const unsigned char *input, size_t inputByteLen);

#endif

#ifndef TEST_ARITH_DATA_H
#define TEST_ARITH_DATA_H

#include "../include/ec.h"

digit_t k[NWORDS_FIELD]={0xddc5adf6dab467a3,0x765b68a791f364,0x423401757c27ace3,0xe835f7ae776d33e};
ec_point_t AC={0}, P={0}, twoP={0}, PpQ={0}, PmQ={0}, kP={0}, PpkQ={0};

AC.z.re[0]=0x1;

P.x.re[0]=0xe2cd67b2c5245429;
P.x.re[1]=0xf7c2f1c59c03503a;
P.x.re[2]=0xba30a2d8574d4963;
P.x.re[3]=0x8bbbfdf34f96e7d;
P.x.im[0]=0xfbfea27033366c29;
P.x.im[1]=0x8c1967b16b1d2b7d;
P.x.im[2]=0x1dbe1ae1d962bce1;
P.x.im[3]=0x3fdfb47b31d7964;
P.z.re[0]=0x1;

Q.x.re[0]=0x261d00e4f8c77a8;
Q.x.re[1]=0x60fb6729eee7f55f;
Q.x.re[2]=0xebc99f1fb441b116;
Q.x.re[3]=0x10b9d46c05930a59;
Q.x.im[0]=0x6a0c25bec4d0bf6a;
Q.x.im[1]=0x6361e9dcae04927e;
Q.x.im[2]=0x281c11479c11894e;
Q.x.im[3]=0x1c5ddc5bcd1f60dc;
Q.z.re[0]=0x1;

twoP.x.re[0]=0xe7def443c13470a9;
twoP.x.re[1]=0x1e8512ae900af666;
twoP.x.re[2]=0xb654db27e0fc7de4;
twoP.x.re[3]=0xb3c8cd89da3b27a;
twoP.x.im[0]=0xfb5ffd4bf8405c7a;
twoP.x.im[1]=0x1d0dcf2acc096468;
twoP.x.im[2]=0x7ff1c5a22300ac86;
twoP.x.im[3]=0x194dde3c6eec4b81;
twoP.z.re[0]=0x1;

PpQ.x.re[0]=0x3f6ff3a991476f75;
PpQ.x.re[1]=0x8668c072f41ebeb;
PpQ.x.re[2]=0xbed4f9ef2f9d0332;
PpQ.x.re[3]=0x14f395622432d70c;
PpQ.x.im[0]=0xd714ff5a75f80c55;
PpQ.x.im[1]=0xfff9b4f35854987a;
PpQ.x.im[2]=0x46fbcc8b103f5abb;
PpQ.x.im[3]=0x15677e93177269e5;
PpQ.z.re[0]=0x1;

PmQ.x.re[0]=0x6454dcfb3244a25b;
PmQ.x.re[1]=0x8e2b965ac3fbf149;
PmQ.x.re[2]=0xeb1d0dc4e276059d;
PmQ.x.re[3]=0x159d9738a730509b;
PmQ.x.im[0]=0xccb1e0d511f7d01e;
PmQ.x.im[1]=0x42042af4ecc87f49;
PmQ.x.im[2]=0x9d1ea2740e883e6c;
PmQ.x.im[3]=0x138e8b8eeb53b28e;
PmQ.z.re[0]=0x1;

kP.x.re[0]=0x6cb1d857bf44ccdb;
kP.x.re[1]=0xe5f51c019a5aaf8d;
kP.x.re[2]=0x1441de9ac02f1620;
kP.x.re[3]=0x1effe8b4f4039ae6;
kP.x.im[0]=0x4a89079091765d8a;
kP.x.im[1]=0x430a2c0edd87dc71;
kP.x.im[2]=0x37c3eaf7f906ef1a;
kP.x.im[3]=0x1f48811e16f8ab4d;
kP.z.re[0]=0x1;

PpkQ.x.re[0]=0x4e2741218dc7da3d;
PpkQ.x.re[1]=0x1542024cbd537b72;
PpkQ.x.re[2]=0xc9c0fa0faf9493f3;
PpkQ.x.re[3]=0x21edd3a306404eed;
PpkQ.x.im[0]=0xc47d2394ad258719;
PpkQ.x.im[1]=0x977f96076f104230;
PpkQ.x.im[2]=0xbb2f5ef2c13d43a9;
PpkQ.x.im[3]=0x25b081e443d94125;
PpkQ.z.re[0]=0x1;

#endif
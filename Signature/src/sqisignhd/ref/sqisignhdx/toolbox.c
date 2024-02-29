
#include <stdio.h>
#include <time.h> 

#include "toolbox.h"

static clock_t global_timer;

clock_t
tic()
{
    global_timer = clock();
    return global_timer;
}

float
tac()
{
    float ms = (1000. * (float) (clock() - global_timer) / CLOCKS_PER_SEC);
    return ms;
}

float
TAC(const char *str)
{
    float ms = (1000. * (float) (clock() - global_timer) / CLOCKS_PER_SEC);
    #ifndef NDEBUG
    printf("%s [%d ms]\n", str, (int) ms);
    #endif
    return ms;
}

float
toc(const clock_t t)
{
    float ms = (1000. * (float) (clock() - t) / CLOCKS_PER_SEC);
    return ms;
}

float
TOC(const clock_t t, const char *str)
{
    float ms = (1000. * (float) (clock() - t) / CLOCKS_PER_SEC);
    #ifndef NDEBUG
    printf("%s [%d ms]\n", str, (int) ms);
    #endif
    return ms;
}





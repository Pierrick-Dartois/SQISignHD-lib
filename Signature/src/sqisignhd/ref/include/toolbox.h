
#ifndef TOOLBOX_H
#define TOOLBOX_H

#include <time.h> 

clock_t tic();
float tac(); /* time in ms since last tic */
float TAC(const char *str); /* same, but prints it with label 'str' */
float toc(const clock_t t); /* time in ms since t */
float TOC(const clock_t t, const char *str); /* same, but prints it with label 'str' */

#endif





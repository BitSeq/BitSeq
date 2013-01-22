#ifndef COMMON_H
#define COMMON_H

#include<string>

using namespace std;

#ifdef BIOC_BUILD

#include <R.h>
#include <R_ext/Utils.h>

#define message(...) Rprintf(__VA_ARGS__)

#else

#include<cstdio>

#define message(...) printf(__VA_ARGS__)
#define warning(...) {fprintf(stderr,"WARNING: ");fprintf(stderr, __VA_ARGS__);}
#define error(...) {fprintf(stderr,"ERROR: ");fprintf(stderr, __VA_ARGS__);}

#endif

void buildTime(char *argv0, string compileDate, string compileTime);

bool progressLog(long cur,long outOf, long steps = 10);

#endif

#ifndef COMMON_H
#define COMMON_H

#include<string>

using std::string;

const char bitseq_version[] = BS_VERSION;

#ifdef BIOC_BUILD

#include <R.h>
#include <R_ext/Utils.h>

#define R_INTERUPT R_CheckUserInterrupt()

#define message(...) Rprintf(__VA_ARGS__)

const long samplesAtOnce = 50;

#else

#include<cstdio>

#define R_INTERUPT

#define message(...) printf(__VA_ARGS__)
#define warning(...) {fprintf(stderr,"WARNING: ");fprintf(stderr, __VA_ARGS__);}
#define error(...) {fprintf(stderr,"ERROR: ");fprintf(stderr, __VA_ARGS__);}

#endif

void buildTime(char *argv0, string compileDate, string compileTime, const char *version = bitseq_version);

bool progressLog(long cur,long outOf, long steps = 10, char nl = '\n');

#endif

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
char * _strdup_printf(char * fmt, ...) {
  char zero[2];
  va_list va;
  va_start(va, fmt);
  int s = vsnprintf(zero, 0, fmt, va);
  va_end(va);
  char * rt = malloc(s);
  va_start(va, fmt);
  vsnprintf(rt, s, fmt, va);
  va_end(va);
  return rt;
}

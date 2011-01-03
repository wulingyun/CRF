#include <R.h>
#include <Rdefines.h>

#define setValues(r, c, v) for (int i = 0; i < length(r); i++) c[i] = v

SEXP getListElement(SEXP list, const char *tag);
void setListElement(SEXP list, int i, const char *tag, SEXP value);

void setDim2(SEXP array, int x1, int x2);
void setDim3(SEXP array, int x1, int x2, int x3);

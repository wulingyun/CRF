#include <R.h>
#include <Rdefines.h>

#define set_zeros(rarray, carray) for (int i = 0; i < length(rarray); i++) carray[i] = 0

SEXP getListElement(SEXP list, const char *tag);
void setListElement(SEXP list, int i, const char *tag, SEXP value);

void setDim2(SEXP array, int x1, int x2);
void setDim3(SEXP array, int x1, int x2, int x3);

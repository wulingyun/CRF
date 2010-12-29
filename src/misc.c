#include "CRF.h"

/* get the list element named tag, or return NULL */

SEXP getListElement(SEXP list, const char *tag)
{
	SEXP value = R_NilValue, names = getAttrib(list, R_NamesSymbol);
	for (int i = 0; i < length(list); i++)
		if (strcmp(CHAR(STRING_ELT(names, i)), tag) == 0)
		{
			value = VECTOR_ELT(list, i);
			break;
		}
	return value;
}

/* set the list element */

void setListElement(SEXP list, int i, const char *tag, SEXP value)
{
	SEXP names = getAttrib(list, R_NamesSymbol);
	if (names == R_NilValue)
	{
		PROTECT(names = NEW_STRING(length(list)));
		SET_STRING_ELT(names, i, mkChar(tag));
		setAttrib(list, R_NamesSymbol, names);
		UNPROTECT(1);
	}
	else
		SET_STRING_ELT(names, i, mkChar(tag));
	SET_VECTOR_ELT(list, i, value);
}

/* set dim of array */

void setDim2(SEXP array, int x1, int x2)
{
	SEXP _dim;
	PROTECT(_dim = NEW_INTEGER(2));
	INTEGER_POINTER(_dim)[0] = x1;
	INTEGER_POINTER(_dim)[1] = x2;
	SET_DIM(array, _dim);
	UNPROTECT(1);
}

void setDim3(SEXP array, int x1, int x2, int x3)
{
	SEXP _dim;
	PROTECT(_dim = NEW_INTEGER(3));
	INTEGER_POINTER(_dim)[0] = x1;
	INTEGER_POINTER(_dim)[1] = x2;
	INTEGER_POINTER(_dim)[2] = x3;
	SET_DIM(array, _dim);
	UNPROTECT(1);
}

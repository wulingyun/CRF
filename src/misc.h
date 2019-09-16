#include <R.h>
#include <Rdefines.h>

/* macros */

#define min(a, b) a < b ? a : b;
#define max(a, b) a > b ? a : b;

/* swap variables */

template <class T>
inline void swap2(T &a, T &b)
{
  T temp = a;
  a = b;
  b = temp;
}

/* initialize the list */

template <class T>
inline void SetValues(SEXP r, T *c, T v)
{
  for (int i = 0; i < length(r); i++)
    c[i] = v;
}

/* get the variable */

inline SEXP GetVar(SEXP env, const char *name)
{
  return findVar(install(name), env);
}

inline SEXP GetVarAsInteger(SEXP env, const char *name)
{
  SEXP temp, var;
  PROTECT(temp = GetVar(env, name));
  var = AS_INTEGER(temp);
  UNPROTECT(1);
  return(var);
}

inline SEXP GetVarAsNumeric(SEXP env, const char *name)
{
  SEXP temp, var;
  PROTECT(temp = GetVar(env, name));
  var = AS_NUMERIC(temp);
  UNPROTECT(1);
  return(var);
}

inline SEXP GetVarAsList(SEXP env, const char *name)
{
  SEXP temp, var;
  PROTECT(temp = GetVar(env, name));
  var = AS_LIST(temp);
  UNPROTECT(1);
  return(var);
}

/* set the variable */

inline void SetVar(SEXP env, const char *name, SEXP var)
{
  defineVar(install(name), var, env);
}

/* get the list element */

inline SEXP GetListElement(SEXP list, int i)
{
  if (i >= 0 && i < length(list)) return VECTOR_ELT(list, i);
  return R_NilValue;
}

/* get the list element named tag, or return NULL */

inline SEXP GetListElement(SEXP list, const char *tag)
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

inline void SetListElement(SEXP list, int i, const char *tag, SEXP value)
{
  SEXP _names = getAttrib(list, R_NamesSymbol);
  if (_names == R_NilValue)
  {
    PROTECT(_names = NEW_STRING(length(list)));
    SET_STRING_ELT(_names, i, mkChar(tag));
    setAttrib(list, R_NamesSymbol, _names);
    UNPROTECT(1);
  }
  else
    SET_STRING_ELT(_names, i, mkChar(tag));
  SET_VECTOR_ELT(list, i, value);
}

/* set dim of array */

inline void SetDim2(SEXP array, int x1, int x2)
{
  SEXP _dim;
  PROTECT(_dim = NEW_INTEGER(2));
  INTEGER_POINTER(_dim)[0] = x1;
  INTEGER_POINTER(_dim)[1] = x2;
  SET_DIM(array, _dim);
  UNPROTECT(1);
}

inline void SetDim3(SEXP array, int x1, int x2, int x3)
{
  SEXP _dim;
  PROTECT(_dim = NEW_INTEGER(3));
  INTEGER_POINTER(_dim)[0] = x1;
  INTEGER_POINTER(_dim)[1] = x2;
  INTEGER_POINTER(_dim)[2] = x3;
  SET_DIM(array, _dim);
  UNPROTECT(1);
}

inline void SetDim4(SEXP array, int x1, int x2, int x3, int x4)
{
  SEXP _dim;
  PROTECT(_dim = NEW_INTEGER(4));
  INTEGER_POINTER(_dim)[0] = x1;
  INTEGER_POINTER(_dim)[1] = x2;
  INTEGER_POINTER(_dim)[2] = x3;
  INTEGER_POINTER(_dim)[3] = x4;
  SET_DIM(array, _dim);
  UNPROTECT(1);
}

/* allocate vector */

template <class T>
inline T *R_allocVector(int n)
{
  T *vector = (T *) R_alloc(n, sizeof(T));
  return vector;
}

/* allocate vector (Calloc version) */

template <class T>
inline T *C_allocVector(int n)
{
  T *vector = (T *) Calloc(n, T);
  return vector;
}

template <class T>
inline void C_freeVector(T *vector)
{
  Free(vector);
}

/* allocate multidimensional array */

template <class T, int n>
inline T **R_allocArray(int dim[n])
{
  int size1, size2;
  T **array, **sub1, **sub2, **tmp;
  size1 = dim[0];
  array = sub1 = (T **) R_alloc(size1, sizeof(T *));
  for (int i = 1; i < n-1; i++)
  {
    size2 = size1 * dim[i];
    tmp = sub2 = (T **) R_alloc(size2, sizeof(T *));
    for (int j = 0; j < size1; j++)
    {
      sub1[j] = (T *) tmp;
      tmp += dim[i];
    }
    size1 = size2;
    sub1 = sub2;
  }
  T *block = (T *) R_alloc(size1 * dim[n-1], sizeof(T));
  for (int i = 0; i < size1; i++)
  {
    sub1[i] = block;
    block += dim[n-1];
  }
  return array;
}

template <class T>
inline T **R_allocArray(int dim1, int dim2)
{
  int dim[] = {dim1, dim2};
  return R_allocArray<T, 2>(dim);
}

template <class T>
inline T **R_allocArray(int dim1, int dim2, int dim3)
{
  int dim[] = {dim1, dim2, dim3};
  return R_allocArray<T, 3>(dim);
}

/* allocate multidimensional array (Calloc version) */

template <class T, int n>
inline T **C_allocArray(int dim[n])
{
  int size1, size2;
  T **array, **sub1, **sub2, **tmp;
  size1 = dim[0];
  array = sub1 = (T **) Calloc(size1, T *);
  for (int i = 1; i < n-1; i++)
  {
    size2 = size1 * dim[i];
    tmp = sub2 = (T **) Calloc(size2, T *);
    for (int j = 0; j < size1; j++)
    {
      sub1[j] = (T *) tmp;
      tmp += dim[i];
    }
    size1 = size2;
    sub1 = sub2;
  }
  T *block = (T *) Calloc(size1 * dim[n-1], T);
  for (int i = 0; i < size1; i++)
  {
    sub1[i] = block;
    block += dim[n-1];
  }
  return array;
}

template <class T>
inline T **C_allocArray(int dim1, int dim2)
{
  int dim[] = {dim1, dim2};
  return C_allocArray<T, 2>(dim);
}

template <class T>
inline T **C_allocArray(int dim1, int dim2, int dim3)
{
  int dim[] = {dim1, dim2, dim3};
  return C_allocArray<T, 3>(dim);
}

template <class T, int n>
inline void C_freeArray(T **array)
{
  void *p;
  for (int i = 0; i < n-1; i++)
  {
    p = array;
    array = (T **) array[0];
    Free(p);
  }
  Free(array);
}

/* allocate 2D array with varied dim2 */

template <class T>
inline T **R_allocArray2(int dim1, int *dim2)
{
  T *block, **array;
  int array_size;
  array_size = 0;
  for (int i = 0; i < dim1; i++)
  {
    if (dim2[i] >= 0)
      array_size += dim2[i];
    else
    {
      array_size = -1;
      break;
    }
  }
  block = (T *) R_alloc(array_size, sizeof(T));
  array = (T **) R_alloc(dim1, sizeof(T *));
  for (int i = 0; i < dim1; i++)
  {
    array[i] = block;
    block += dim2[i];
  }
  return array;
}

/* allocate 2D array with varied dim2 (Calloc version) */

template <class T>
inline T **C_allocArray2(int dim1, int *dim2)
{
  T *block, **array;
  int array_size;
  array_size = 0;
  for (int i = 0; i < dim1; i++)
  {
    if (dim2[i] >= 0)
      array_size += dim2[i];
    else
    {
      array_size = -1;
      break;
    }
  }
  block = (T *) Calloc(array_size, T);
  array = (T **) Calloc(dim1, T *);
  for (int i = 0; i < dim1; i++)
  {
    array[i] = block;
    block += dim2[i];
  }
  return array;
}

template <class T>
inline void C_freeArray2(T **array)
{
  Free(array[0]);
  Free(array);
}

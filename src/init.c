#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "interp_c.h"

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

static R_NativePrimitiveArgType biliip_t[10] = {
  REALSXP, /* X0,  */
  REALSXP, /* Y0,  */
  REALSXP, /* Z0,  */
  INTSXP,  /* N0,  */
  REALSXP, /* X,   */
  REALSXP, /* Y,   */
  REALSXP, /* Z,   */
  INTSXP,  /* NX,  */
  INTSXP,  /* NY   */
  INTSXP   /* IER  */
};

/* .Call calls */
extern SEXP _interp_aSpline(SEXP, SEXP, SEXP, SEXP,SEXP);
extern SEXP _interp_inHull(SEXP, SEXP, SEXP, SEXP);
extern SEXP _interp_interpDeltri(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _interp_interpShull(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                                SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _interp_left(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _interp_nearestNeighbours(SEXP, SEXP);
extern SEXP _interp_on(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _interp_onHull(SEXP, SEXP, SEXP, SEXP);
extern SEXP _interp_partDerivGrid(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _interp_partDerivPoints(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _interp_shullDeltri(SEXP, SEXP, SEXP);
extern SEXP _interp_triFind(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _interp_circum(SEXP, SEXP);
extern SEXP _interp_BiLinear(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _interp_ConvexHull(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_interp_aSpline",           (DL_FUNC) &_interp_aSpline,            5},
    {"_interp_inHull",            (DL_FUNC) &_interp_inHull,             4},
    {"_interp_BiLinear",          (DL_FUNC) &_interp_BiLinear,           5},
    {"_interp_interpDeltri",      (DL_FUNC) &_interp_interpDeltri,       6},
    {"_interp_interpShull",       (DL_FUNC) &_interp_interpShull,       18},
    {"_interp_left",              (DL_FUNC) &_interp_left,               7},
    {"_interp_nearestNeighbours", (DL_FUNC) &_interp_nearestNeighbours,  2},
    {"_interp_on",                (DL_FUNC) &_interp_on,                 7},
    {"_interp_onHull",            (DL_FUNC) &_interp_onHull,             4},
    {"_interp_partDerivGrid",     (DL_FUNC) &_interp_partDerivGrid,     12},
    {"_interp_partDerivPoints",   (DL_FUNC) &_interp_partDerivPoints,   12},
    {"_interp_shullDeltri",       (DL_FUNC) &_interp_shullDeltri,        3},
    {"_interp_triFind",           (DL_FUNC) &_interp_triFind,            8},
    {"_interp_circum",            (DL_FUNC) &_interp_circum,             2},
    {"_interp_ConvexHull",        (DL_FUNC) &_interp_ConvexHull,         2},
    {NULL, NULL, 0}
};

static R_FortranMethodDef fortranMethods[] = {
  {"biliip", (DL_FUNC) &F77_SUB(biliip), 10, biliip_t}, /* bilinear    */
  {NULL, NULL, 0}
};

void R_init_interp(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, fortranMethods, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

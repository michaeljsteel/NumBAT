//------------------------------------------------------------------------------
// UMFPACK/Demo/umf4_f77zwrapper: Fortran interface for UMFPACK
//------------------------------------------------------------------------------

// UMFPACK, Copyright (c) 2005-2022, Timothy A. Davis, All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0+

//------------------------------------------------------------------------------

/* FORTRAN interface for the C-callable UMFPACK library (complex / int version
 * only and complex / int64_t versions only).  This is HIGHLY
 * non-portable.  You will need to modify this depending on how your FORTRAN
 * and C compilers behave.
 *
 * See umf4_f77wrapper.c for more information.
 *
 * The complex values are provided in two separate arrays.  Ax contains the
 * real part and Az contains the imaginary part.  The solution vector is in
 * x (the real part) and xz (the imaginary part.  b is the real part of the
 * right-hand-side and bz is the imaginary part.  Does not support the
 * packed complex type.
 */

#include "inttypes.h"
#include "umfpack.h"
#include <ctype.h>
#include <stdio.h>
#ifdef NULL
#undef NULL
#endif
#define NULL 0
#define LEN 200

/* -------------------------------------------------------------------------- */
/* integer type: int32_t or int64_t */
/* -------------------------------------------------------------------------- */

#if defined (ZLONG)

#define Int int64_t
#define UMFPACK_defaults	 umfpack_zl_defaults
#define UMFPACK_free_numeric	 umfpack_zl_free_numeric
#define UMFPACK_free_symbolic	 umfpack_zl_free_symbolic
#define UMFPACK_numeric		 umfpack_zl_numeric
#define UMFPACK_report_control	 umfpack_zl_report_control
#define UMFPACK_report_info	 umfpack_zl_report_info
#define UMFPACK_save_numeric	 umfpack_zl_save_numeric
#define UMFPACK_save_symbolic	 umfpack_zl_save_symbolic
#define UMFPACK_load_numeric	 umfpack_zl_load_numeric
#define UMFPACK_load_symbolic	 umfpack_zl_load_symbolic
#define UMFPACK_scale		 umfpack_zl_scale
#define UMFPACK_solve		 umfpack_zl_solve
#define UMFPACK_symbolic	 umfpack_zl_symbolic

#else

#define Int int
#define UMFPACK_defaults	 umfpack_zi_defaults
#define UMFPACK_free_numeric	 umfpack_zi_free_numeric
#define UMFPACK_free_symbolic	 umfpack_zi_free_symbolic
#define UMFPACK_numeric		 umfpack_zi_numeric
#define UMFPACK_report_control	 umfpack_zi_report_control
#define UMFPACK_report_info	 umfpack_zi_report_info
#define UMFPACK_save_numeric	 umfpack_zi_save_numeric
#define UMFPACK_save_symbolic	 umfpack_zi_save_symbolic
#define UMFPACK_load_numeric	 umfpack_zi_load_numeric
#define UMFPACK_load_symbolic	 umfpack_zi_load_symbolic
#define UMFPACK_scale		 umfpack_zi_scale
#define UMFPACK_solve		 umfpack_zi_solve
#define UMFPACK_symbolic	 umfpack_zi_symbolic

#endif

/* -------------------------------------------------------------------------- */
/* construct a file name from a file number (not user-callable) */
/* -------------------------------------------------------------------------- */

static void make_filename (Int filenum, char *prefix, char *filename)
{
    char *psrc, *pdst ;
#ifdef ZLONG
    sprintf (filename, "%s%"PRId64".umf", prefix, filenum) ;
#else
    sprintf (filename, "%s%d.umf", prefix, filenum) ;
#endif
    /* remove any spaces in the filename */
    pdst = filename ;
    for (psrc = filename ; *psrc ; psrc++)
    {
	if (!isspace (*psrc)) *pdst++ = *psrc ;
    }
    *pdst = '\0' ;
}

/* ========================================================================== */
/* === with underscore ====================================================== */
/* ========================================================================== */

/* Solaris, Linux, and SGI IRIX.  Probably Compaq Alpha as well. */

#define umf4zdef_FORTRAN  UMF4ZDEF
#define umf4zpcon_FORTRAN UMF4ZPCON
#define umf4zsym_FORTRAN  UMF4ZSYM
#define umf4znum_FORTRAN  UMF4ZNUM
#define umf4zsolr_FORTRAN UMF4ZSOLR
#define umf4zsol_FORTRAN  UMF4ZSOL
#define umf4zscal_FORTRAN UMF4ZSCAL
#define umf4zpinf_FORTRAN UMF4ZPINF
#define umf4zfnum_FORTRAN UMF4ZFNUM
#define umf4zfsym_FORTRAN UMF4ZFSYM
#define umf4zsnum_FORTRAN UMF4ZSNUM
#define umf4zssym_FORTRAN UMF4ZSSYM
#define umf4zlnum_FORTRAN UMF4ZLNUM
#define umf4zlsym_FORTRAN UMF4ZLSYM

#include "umf4_f77zwrapper_methods.h"



#ifndef numgrid_EXPORT_H
#define numgrid_EXPORT_H

#ifdef numgrid_STATIC_DEFINE
#  define numgrid_EXPORT
#  define numgrid_NO_EXPORT
#else
#  ifndef numgrid_EXPORT
#    ifdef numgrid_EXPORTS
        /* We are building this library */
#      define numgrid_EXPORT __attribute__((visibility("default")))
#    else
        /* We are using this library */
#      define numgrid_EXPORT __attribute__((visibility("default")))
#    endif
#  endif

#  ifndef numgrid_NO_EXPORT
#    define numgrid_NO_EXPORT __attribute__((visibility("hidden")))
#  endif
#endif

#ifndef numgrid_DEPRECATED
#  define numgrid_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef numgrid_DEPRECATED_EXPORT
#  define numgrid_DEPRECATED_EXPORT numgrid_EXPORT numgrid_DEPRECATED
#endif

#ifndef numgrid_DEPRECATED_NO_EXPORT
#  define numgrid_DEPRECATED_NO_EXPORT numgrid_NO_EXPORT numgrid_DEPRECATED
#endif

#if 1 /* DEFINE_NO_DEPRECATED */
#  ifndef numgrid_NO_DEPRECATED
#    define numgrid_NO_DEPRECATED
#  endif
#endif

#endif

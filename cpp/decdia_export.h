
#ifndef DECDIA_EXPORT_H
#define DECDIA_EXPORT_H

#ifdef DECDIA_STATIC_DEFINE
#  define DECDIA_EXPORT
#  define DECDIA_NO_EXPORT
#else
#  ifndef DECDIA_EXPORT
#    ifdef decDIA_EXPORTS
        /* We are building this library */
#      ifdef _MSC_VER
#      		define DECDIA_EXPORT __declspec(dllexport)
#      else
#			define DECDIA_EXPORT
#      endif
#    else
        /* We are using this library */
#      ifdef _MSC_VER
#      define DECDIA_EXPORT __declspec(dllimport)
#      else
#			define DECDIA_EXPORT
#      endif
#    endif
#  endif

#  ifndef DECDIA_NO_EXPORT
#    define DECDIA_NO_EXPORT 
#  endif
#endif

#ifndef DECDIA_DEPRECATED
#  define DECDIA_DEPRECATED __declspec(deprecated)
#endif

#ifndef DECDIA_DEPRECATED_EXPORT
#  define DECDIA_DEPRECATED_EXPORT DECDIA_EXPORT DECDIA_DEPRECATED
#endif

#ifndef DECDIA_DEPRECATED_NO_EXPORT
#  define DECDIA_DEPRECATED_NO_EXPORT DECDIA_NO_EXPORT DECDIA_DEPRECATED
#endif

#define DEFINE_NO_DEPRECATED 0
#if DEFINE_NO_DEPRECATED
# define DECDIA_NO_DEPRECATED
#endif

#endif

#pragma once

#if !defined(OPTIMIZATION_STATIC_DEFINE)
#  if defined(_WIN32)
#    if defined(OPTIMIZATION_BUILDING_LIBRARY)
#      define OPTIMIZATION_API __declspec(dllexport)
#    else
#      define OPTIMIZATION_API __declspec(dllimport)
#    endif
#  elif defined(__GNUC__) || defined(__clang__)
#    define OPTIMIZATION_API __attribute__((visibility("default")))
#  else
#    define OPTIMIZATION_API
#  endif
#else
#  define OPTIMIZATION_API
#endif

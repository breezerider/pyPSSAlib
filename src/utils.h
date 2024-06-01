#include <boost/config.hpp>


#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

// inspired by https://stackoverflow.com/questions/30670909/boostpython-extract-c-value-from-numpy-ndarray-return-value
#define PY_ERRMSG(msg) \
  { PyErr_SetString(PyExc_TypeError, (boost::format("PY_ERRMSG(%1%:%2%): %3%") % (__FILE__) % (__LINE__) % (msg)).str().c_str()); };

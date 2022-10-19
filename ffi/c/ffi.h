#ifndef EPA_FFI_H
#define EPA_FFI_H

// Closure
struct epa_function {
  void (*function)();
  void* data;
  void (*destructor)(void*);

#ifdef __cplusplus
  ~epa_function() { if (destructor) destructor(data); };
#endif
};

typedef struct epa_function epa_function;

// XXX: non-standard features: statement expressions and __VA_OPT__.
#define epa_call(type, func, ...) \
  ({ \
     epa_function* _epa_func = (func); \
     ((type)_epa_func->function)(__VA_ARGS__ __VA_OPT__(,) _epa_func->data); \
  })

#ifdef __cplusplus
extern "C" {
#endif

void* epa_get_error();
int epa_get_error_layer();
void epa_set_error(void* error, int layer);
void epa_clear_error();
const char* epa_cpp_error_message(void* error);

const int epa_cpp_error_layer = 0x002b2b43; // "C++"

epa_function* epa_make_function(
    void (*function)(),
    void* data,
    void (*destructor)(void*)
);

void epa_destroy_function(epa_function*);

#ifdef __cplusplus
};
#endif

#endif // EPA_FFI_H

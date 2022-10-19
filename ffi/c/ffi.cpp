#include <exception>

#include "ffi.hpp"

namespace epa {
namespace ffi {

Error error { nullptr, 0 };

void set_error() {
  error.error = new std::exception_ptr(std::current_exception());
  error.layer = epa_cpp_error_layer;
};

extern "C" void* epa_get_error() {
  return error.error;
};

extern "C" int epa_get_error_layer() {
  return error.layer;
};

extern "C" void epa_set_error(void* err, int layer) {
  epa_clear_error();
  error.error = err;
  error.layer = layer;
};

extern "C" void epa_clear_error() {
  if (!error.error) return;
  if (error.layer == epa_cpp_error_layer)
    delete reinterpret_cast<std::exception_ptr*>(error.error);
  error.error = nullptr;
  error.layer = 0;
};

extern "C" const char* epa_cpp_error_message(void* error) {
  try {
    std::rethrow_exception(*reinterpret_cast<std::exception_ptr*>(error));
  } catch (std::exception& e) {
    return e.what();
  };
};

extern "C" Function* epa_make_function(
    void (*function)(),
    void* data,
    void (*destructor)(void*)
) {
  return new Function { function, data, destructor };
};

extern "C" void epa_destroy_function(Function* function) {
  delete function;
};

}; // namespace ffi
}; // namespace epa

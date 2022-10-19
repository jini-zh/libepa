#include <functional>
#include <stdexcept>
#include <type_traits>

#include "ffi.h"

#define FFI_CATCH_R(result) \
  catch (ForeignError&) { \
    return result; \
  } catch (std::exception&) { \
    set_error(); \
    return result; \
  }

#define FFI_CATCH FFI_CATCH_R(nullptr)

namespace epa {
namespace ffi {

using Function = epa_function;

struct Error {
  void* error;
  int layer;
};

extern Error error;

void set_error();

class ForeignError: public std::exception {
  public:
    ForeignError(): error_(epa::ffi::error) {};

    const Error& error() const { return error_; };
  private:
    Error error_;
};

template <typename>
struct is_std_function_ {
  static const bool value = false;
};

template <typename Signature>
struct is_std_function_<std::function<Signature>> {
  static const bool value = true;
};

template <typename T>
inline constexpr bool is_std_function
  = is_std_function_<std::remove_cvref_t<T>>::value;

template <typename>
struct is_function_ {
  static const bool value = false;
};

template <typename Result, typename... Args>
struct is_function_<Result (Args...)> {
  static const bool value = true;
};

template <typename T>
inline constexpr bool is_function = is_function_<std::remove_cvref_t<T>>::value;

template <typename Signature> struct std_function_signature_;

template <typename Signature>
struct std_function_signature_<std::function<Signature>> {
  using signature = Signature;
};

template <typename F>
using std_function_signature = std_function_signature_<F>::signature;

template <typename T>
using lower_t = std::conditional_t<is_std_function<T>, Function*, T>;

template <typename> struct lower_;

template <typename T>
inline
std::enable_if_t<!std::is_same_v<std::remove_cvref_t<T>, Function*>, T&>
lower(T& x) {
  return x;
};

template <typename Signature>
inline
std::enable_if_t<is_function<Signature>, std::function<Signature>>
lower(Function* f) {
  return lower_<Signature>(f);
};

template <typename F>
inline
std::enable_if_t<is_std_function<F>, F>
lower(Function* f) {
  return lower<std_function_signature<F>>(f);
};

template <typename Signature>
void destructor(std::function<Signature>* function) {
  delete function;
};

template <typename Result, typename... Args>
lower_t<Result>
trampoline(lower_t<Args>... args, std::function<Result (Args...)>* f);

template <typename T>
inline
std::enable_if_t<!is_std_function<T>, T&&>
lift(T&& x) {
  return std::forward<T>(x);
};

template <typename Result, typename... Args>
Function*
lift(std::function<Result (Args...)>&& f) {
  try {
    auto pf = new std::function<Result (Args...)>(std::move(f));
    try {
      return new Function {
        reinterpret_cast<void (*)()>(trampoline<Result, Args...>),
        pf,
        reinterpret_cast<void (*)(void*)>(destructor<Result (Args...)>)
      };
    } catch (...) {
      delete pf;
      throw;
    };
  } catch (std::exception&) {
    set_error();
    return nullptr;
  };
};

template <typename T> struct lift_on_stack_ {
  T& x;
  lift_on_stack_(T& x): x(x) {};
  operator T&() { return x; };
};

template <typename Result, typename... Args>
struct lift_on_stack_<const std::function<Result (Args...)>&> {
  Function f;

  lift_on_stack_(const std::function<Result (Args...)>& function) {
    f.function   = reinterpret_cast<void (*)()>(trampoline<Result, Args...>);
    f.data       = &const_cast<std::function<Result (Args...)>&>(function);
    f.destructor = nullptr;
  };

  operator Function*() { return &f; };
};

template <typename Result, typename... Args>
lower_t<Result>
trampoline(lower_t<Args>... args, std::function<Result (Args...)>* f) {
  try {
    return lift((*f)(lower<std::remove_cvref_t<Args>>(args)...));
  } catch (ForeignError& e) {
    error = e.error();
    return lower_t<Result>();
  } catch (std::exception&) {
    set_error();
    return lower_t<Result>();
  };
};

template <typename Result, typename... Args>
struct lower_<Result (Args...)> {
  using F = std::function<Result (Args...)>;
  using T = lower_t<Result> (*)(lower_t<Args>..., F*);

  Function* f;

  lower_(Function* f): f(f) {};

  operator std::function<Result (Args...)>() const {
    if (reinterpret_cast<T>(f->function) == trampoline<Result, Args...>)
      return *reinterpret_cast<F*>(f->data);
    return [f = f](Args... args) -> Result {
      lower_t<Result> result = reinterpret_cast<T>(f->function)(
          lift_on_stack_<Args>(args)...,
          static_cast<F*>(f->data)
      );
      if (error.error) throw ForeignError();
      return lower<Result>(result);
    };
  };
};

}; // namespace epa
}; // namespace ffi

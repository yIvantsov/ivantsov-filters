#pragma once

#include <array>
#include <cmath>
#include <concepts>
#include <numbers>
#include <type_traits>
#include <utility>

namespace ivantsov {
template<typename Sample>
using sample_type = typename decltype([] {
     if constexpr (requires { typename Sample::value_type; }) {
          return std::type_identity<typename Sample::value_type> {};
     } else {
          return std::type_identity<Sample> {};
     }
}())::type;

enum struct Correction { None, Omega, Sigma };

namespace Linear {
enum struct Type { HighPass, LowPass, AllPass, HighShelf, LowShelf, BandPass, BandStop, BandShelf };
enum struct Param { Fc, Sigma, Gain, Zeta };

using enum Type;
using enum Param;

template<auto Type>
concept is_band = Type == BandPass || Type == BandStop || Type == BandShelf;

template<auto Type>
concept is_shelf = Type == HighShelf || Type == LowShelf || Type == BandShelf;

template<auto Type, auto Order, typename Sample, auto Correction = Correction::None>
     requires std::floating_point<sample_type<Sample>> && (Order > 0) && (!is_band<Type> || Order % 2 == 0)
struct Filter {
     using T = sample_type<Sample>;

     auto initialise(const T x) {
          fs = x;
          reset();
     }

     template<auto N = 0>
     auto processed(const Sample x) {
          if constexpr (N == 0 && Order % 2 != 0) { return processed<N + 1>(processed6(x)); }
          if constexpr (N < (Order + 1) / 2) { return processed<N + 1>(processed12<N>(x)); }
          return x;
     }

     auto reset() {
          z.fill({Sample {T {}}, Sample {T {}}});
     }

     template<auto... Param>
     auto set(const decltype(Param, T {})... x)
          requires(... && (Param == Fc || (Param == Sigma && Correction == Correction::None) ||
                           (Param == Gain && is_shelf<Type>) || (Param == Zeta && Order > 1))) {
          (..., [&] {
               if constexpr (Param == Fc) { w = fs / (std::numbers::pi_v<T> * x); }
               if constexpr (Param == Sigma) { sigma = x; }
               if constexpr (Param == Gain) { gain = std::pow(x, T {1} / Order); }
               if constexpr (Param == Zeta) { zeta = x; }
               if constexpr (Param == Fc && Correction != Correction::None) { correct_cutoff(); }
          }());
          update();
     }

     template<auto N = 0>
     auto magnitude(const T x) const {
          if constexpr (N == 0 && Order % 2 != 0) { return magnitude6(x) * magnitude<N + 1>(x); }
          if constexpr (N < (Order + 1) / 2) { return magnitude12<N>(x) * magnitude<N + 1>(x); }
          return T {1};
     }

 private:
     auto processed6(const Sample x) {
          const auto theta {(x - z[0][0]) * b[0][0]};
          const auto y {Type == HighPass ? theta * b[0][1] : theta * b[0][1] + z[0][0]};
          z[0][0] = z[0][0] + theta;
          return Type == LowShelf ? y * gain : y;
     }

     template<auto N>
     auto processed12(const Sample x) {
          const auto theta {(x - z[N][0] - z[N][1] * b[N][1]) * b[N][0]};
          const auto y {(Type == HighPass || Type == BandPass) ? theta * b[N][3] + z[N][1] * b[N][2]
                                                               : theta * b[N][3] + z[N][1] * b[N][2] + z[N][0]};
          z[N][0] = z[N][0] + theta;
          z[N][1] = -z[N][1] - theta * b[N][1];
          return Type == LowShelf ? y * (gain * gain) : y;
     }

     auto magnitude6(const T x) const {
          const auto sq {[](const auto t) { return t * t; }};
          const auto d {1 - x};
          const auto F {[x, d, sq](const auto t) { return 1 + x + sq(t + t - 1) * d; }};
          const auto num {Type == HighPass ? 4 * d * sq(b[0][1]) : F(b[0][1])};
          return std::sqrt(num / F(1 / b[0][0])) * (Type == LowShelf ? gain : 1);
     }

     template<auto N>
     auto magnitude12(const T x) const {
          const auto sq {[](const auto t) { return t * t; }};
          const auto d {1 - x};
          const auto F {
          [x, d, sq](const auto t, const auto u) { return sq(1 + x - u * d) + sq(t + t - u - 1) * d * (1 + x); }};
          const T num {Type == HighPass ? 4 * sq(d * b[N][3])
                       : Type == BandPass
                         ? d * (sq(b[N][1] * b[N][2]) * d + sq(2 * b[N][3] - b[N][1] * b[N][2]) * (1 + x))
                         : F(b[N][3], b[N][1] * b[N][2])};
          return std::sqrt(num / F(1 / b[N][0], sq(b[N][1]))) * (Type == LowShelf ? gain * gain : 1);
     }

     template<auto N = 0>
     auto update() {
          if constexpr (N == 0 && Order % 2 != 0) { update6(); }
          if constexpr (N >= Order % 2) { update12<N>(); }
          if constexpr (N + 1 < (Order + 1) / 2) { update<N + 1>(); }
     }

     auto update6() {
          const auto v {[this](const auto t) { return std::sqrt(t + sigma * sigma); }};
          b[0][0] = Type == HighShelf ? 2 / (1 + v(w * w / gain))
                  : Type == LowShelf  ? 2 / (1 + v(w * w * gain))
                                      : 2 / (1 + v(w * w));
          b[0][1] = Type == HighShelf ? (1 + v(w * w * gain)) / 2
                  : Type == LowShelf  ? (1 + v(w * w / gain)) / 2
                  : Type == AllPass   ? (1 - v(w * w)) / 2
                  : Type == LowPass   ? (1 + sigma) / 2
                                      : w / 2;
     }

     template<auto N>
     auto section(const T x) const {
          if constexpr (Order == 2 || (is_band<Type> && 2 * N + 1 == Order / 2)) { return std::pair {w, x}; }
          if constexpr (is_band<Type>) {
               const auto [s, d] {std::array {x * std::sin(std::numbers::pi_v<T> * (2 * N + 1) / Order), 1 + x * x}};
               const auto u {std::sqrt((d + std::sqrt(d * d - 4 * s * s)) / 2)};
               const auto radius {u + std::sqrt(u * u - 1)};
               return std::pair {N < Order / 4 ? w / radius : w * radius, s / u};
          }
          const auto zeta_N {std::cos(std::numbers::pi_v<T> * (2 * N + 1 - Order % 2) / (2 * Order))};
          return std::pair {w, zeta_N * (N + 1 == (Order + 1) / 2 ? x * std::numbers::sqrt2_v<T> : 1)};
     }

     template<auto N>
     auto update12() {
          const auto [w_N, zeta_N] {section<N>(Type == BandShelf ? zeta / gain : zeta)};
          const auto [w_sq, sigma_sq, zeta_sq] {std::array {w_N * w_N, sigma * sigma, zeta_N * zeta_N}};
          const auto vk {[sigma_sq](const auto t, const auto u) {
               const auto k {t * (u + u) - t};
               return std::pair {std::sqrt(t * t + sigma_sq * (k + k + sigma_sq)), k + sigma_sq};
          }};
          const auto [v, k] = Type == LowShelf  ? vk(w_sq * gain, zeta_sq)
                            : Type == HighShelf ? vk(w_sq / gain, zeta_sq)
                                                : vk(w_sq, zeta_sq);
          b[N][0] = 2 / (1 + std::sqrt(2 * (k + v)) + v);
          b[N][1] = std::sqrt(v);
          const auto [w_a, zeta_a] {Type == BandShelf ? section<N>(zeta * gain) : std::pair {w_N, zeta_N}};
          const auto [v_a, k_a] = Type == LowShelf  ? vk(w_sq / gain, zeta_sq)
                                : Type == HighShelf ? vk(w_sq * gain, zeta_sq)
                                : Type == BandShelf ? vk(w_a * w_a, zeta_a * zeta_a)
                                                    : std::pair {v, k};
          b[N][2] = Type == HighPass ? w_sq / b[N][1]
                  : Type == BandPass ? 2 * w_N * zeta * sigma / b[N][1]
                  : Type == LowPass  ? sigma_sq / b[N][1]
                  : Type == BandStop ? (w * w - sigma_sq) / b[N][1]
                  : Type == AllPass  ? b[N][1]
                                     : v_a / b[N][1];
          b[N][3] = Type == HighPass ? w_sq / 2
                  : Type == BandPass ? w_N * zeta * (1 + sigma)
                  : Type == LowPass  ? (1 + sigma) * (1 + sigma) / 2
                  : Type == BandStop ? (1 + w * w - sigma_sq) / 2
                  : Type == AllPass  ? (1 - std::sqrt(2 * (k + v)) + v) / 2
                                     : (1 + std::sqrt(2 * (k_a + v_a)) + v_a) / 2;
     }

     auto correct_cutoff() {
          const std::common_type_t<T, double> v {w}; // Sigma correction can fail on float for small fc.
          const auto ctg {1 / std::tan(1 / v)};
          if constexpr (Correction == Correction::Omega) { w = std::sqrt(sigma * sigma + ctg * ctg); }
          if constexpr (Correction == Correction::Sigma) { sigma = std::sqrt(v * v - ctg * ctg); }
     }

     T fs {};
     T w {};
     T zeta {1 / std::numbers::sqrt2_v<T>};
     T sigma {2 / std::numbers::pi_v<T>};
     T gain {1};
     std::array<std::array<T, 4>, (Order + 1) / 2> b {};
     std::array<std::array<Sample, 2>, (Order + 1) / 2> z {};
};
} // namespace Linear
} // namespace ivantsov

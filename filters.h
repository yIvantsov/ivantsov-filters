#pragma once

#include <array>
#include <cmath>
#include <numbers>

namespace ivantsov
{
    using namespace std::numbers;

    template<typename T, typename S>
    concept is_algebraic = requires {
        { (((S {} * T {} + S {}) - S {}) + T {}) - T {} } -> std::same_as<S>;
    };

    enum struct Warp { None, Sigma };

    namespace Linear::FirstOrder::Details
    {
        enum struct Type { HighPass, LowPass, AllPass, HighShelf, LowShelf };
        using enum Type;

        template<Type Type, std::floating_point T, typename Sample, Warp Warp>
        requires is_algebraic<T, Sample>
        struct StateSpace
        {
            static constexpr bool is_shelf {Type == HighShelf || Type == LowShelf};

            auto initialize(const auto sample_rate) { fs = sample_rate; }

            auto processed(const auto x)
            {
                const auto theta {(x - z) * b[0]};
                const auto y {Type == HighPass ? theta * b[1] : theta * b[1] + z};
                z = z + theta;
                return Type == LowShelf ? y * b[2] : y;
            }

            auto reset() { z = {}; }

            auto set_fc(const auto x)
            {
                w = fs / (T {2} * pi_v<T> * x);
                if (Warp == Warp::Sigma && w > inv_pi_v<T>) warp_sigma();
                update_blocks();
            }

            auto set_g(const auto x) requires is_shelf
            {
                b[2] = x;
                update_blocks();
            }

            auto set_sigma(const auto x) requires(Warp == Warp::None)
            {
                sigma = x / (T {2} * pi_v<T>);
                update_blocks();
            }

        private:
            auto update_blocks()
            {
                const auto v {[this](const auto x) { return std::sqrt(x + sigma * sigma); }};
                b[0] = Type == HighShelf ? T {1} / (T {0.5} + v(w * w / b[2]))
                     : Type == LowShelf  ? T {1} / (T {0.5} + v(w * w * b[2]))
                                         : T {1} / (T {0.5} + v(w * w));
                b[1] = Type == HighShelf ? T {0.5} + v(w * w * b[2])
                     : Type == LowShelf  ? T {0.5} + v(w * w / b[2])
                     : Type == AllPass   ? T {0.5} - v(w * w)
                     : Type == LowPass   ? T {0.5} + sigma
                                         : w;
            }

            // scaled approximation of [paper1, eq (2.11)]. Max error: 1%.
            auto warp_sigma() { sigma = T {0.40824999} * (T {0.05843357} - w * w) / (T {0.04593294} - w * w); }

            T sigma {inv_pi_v<T>};
            T fs {};
            T w {};
            std::array<T, 2 + is_shelf> b {T {}, T {0.5} + sigma};
            Sample z {};
        };
    } // namespace Linear::FirstOrder::Details

    namespace Linear::SecondOrder::Details
    {
        enum struct Type { HighPass, BandPass, LowPass, AllPass, Notch, HighShelf, LowShelf, MidShelf };
        using enum Type;

        template<Type Type, std::floating_point T, typename Sample, Warp Warp>
        requires is_algebraic<T, Sample>
        struct StateSpace
        {
            static constexpr auto is_shelf {Type == LowShelf || Type == HighShelf || Type == MidShelf};

            auto initialize(const auto sample_rate) { fs = sample_rate; }

            auto processed(const auto x)
            {
                const auto theta {(x - z[0] - z[1] * b[1]) * b[0]};
                const auto y {(Type == HighPass || Type == BandPass) ? theta * b[3] + z[1] * b[2]
                                                                     : theta * b[3] + z[1] * b[2] + z[0]};
                z = {z[0] + theta, Sample {} - z[1] - theta * b[1]};
                return Type == LowShelf ? y * b[4] : y;
            }

            auto reset() { z = {}; }

            auto set_fc(const auto x)
            {
                w = fs / (sqrt2_v<T> * pi_v<T> * x);
                if (Warp == Warp::Sigma && w > inv_pi_v<T> * sqrt2_v<T>) warp_sigma();
                update_blocks();
            }

            auto set_damping(const auto x)
            {
                zeta = x;
                update_blocks();
            }

            auto set_g(const auto x) requires is_shelf
            {
                b[4] = x;
                update_blocks();
            }

            auto set_sigma(const auto x) requires(Warp == Warp::None)
            {
                sigma = x / (sqrt2_v<T> * pi_v<T>);
                update_blocks();
            }

        private:
            auto update_blocks()
            {
                const auto w_sq {w * w};
                const auto sigma_sq {sigma * sigma};
                const auto zeta_sq {zeta * zeta};
                const auto vk {[sigma_sq](const auto x, const auto y) {
                    const auto t {x * (y + y) - x};
                    return std::pair {std::sqrt(x * x + sigma_sq * (t + t + sigma_sq)), t + sigma_sq};
                }};
                const auto [v, k] {Type == LowShelf    ? vk(w_sq * std::sqrt(b[4]), zeta_sq)
                                   : Type == HighShelf ? vk(w_sq / std::sqrt(b[4]), zeta_sq)
                                   : Type == MidShelf  ? vk(w_sq, zeta_sq / b[4])
                                                       : vk(w_sq, zeta_sq)};
                const auto [v_a, k_a] {Type == LowShelf    ? vk(w_sq / std::sqrt(b[4]), zeta_sq)
                                       : Type == HighShelf ? vk(w_sq * std::sqrt(b[4]), zeta_sq)
                                                           : vk(w_sq, zeta_sq * b[4])};
                b[0] = T {1} / (v + std::sqrt(v + k) + T {0.5});
                b[1] = std::sqrt(v + v);
                b[2] = Type == HighPass ? T {2} * w_sq / b[1]
                     : Type == BandPass ? T {4} * w * zeta * sigma / b[1]
                     : Type == LowPass  ? T {2} * sigma_sq / b[1]
                     : Type == Notch    ? T {2} * (w_sq - sigma_sq) / b[1]
                     : Type == AllPass  ? b[1]
                                        : T {2} * v_a / b[1];
                b[3] = Type == HighPass ? w_sq
                     : Type == BandPass ? T {2} * w * zeta * (sigma + T {1} / sqrt2_v<T>)
                     : Type == LowPass  ? T {0.5} + sigma_sq + sqrt2_v<T> * sigma
                     : Type == Notch    ? T {0.5} + w_sq - sigma_sq
                     : Type == AllPass  ? T {0.5} + v - std::sqrt(v + k)
                                        : T {0.5} + v_a + std::sqrt(v_a + k_a);
            }

            // scaled approximation of [paper1, eq (2.11)]. Max error: 1%.
            auto warp_sigma() { sigma = T {0.57735268} * (T {0.11686715} - w * w) / (T {0.09186588} - w * w); }

            T fs {};
            T w {};
            T zeta {};
            T sigma {sqrt2_v<T> * inv_pi_v<T>};
            std::array<T, 4 + is_shelf> b {};
            std::array<Sample, 2> z {};
        };
    } // namespace Linear::SecondOrder::Details

    namespace Linear::FirstOrder
    {
        template<typename T, typename Sample = T, Warp Warp = Warp::None>
        using LowPass = Details::StateSpace<Details::LowPass, T, Sample, Warp>;
        template<typename T, typename Sample = T, Warp Warp = Warp::None>
        using HighPass = Details::StateSpace<Details::HighPass, T, Sample, Warp>;
        template<typename T, typename Sample = T, Warp Warp = Warp::None>
        using AllPass = Details::StateSpace<Details::AllPass, T, Sample, Warp>;
        template<typename T, typename Sample = T, Warp Warp = Warp::None>
        using LowShelf = Details::StateSpace<Details::LowShelf, T, Sample, Warp>;
        template<typename T, typename Sample = T, Warp Warp = Warp::None>
        using HighShelf = Details::StateSpace<Details::HighShelf, T, Sample, Warp>;
    } // namespace Linear::FirstOrder

    namespace Linear::SecondOrder
    {
        template<typename T, typename Sample = T, Warp Warp = Warp::None>
        using LowPass = Details::StateSpace<Details::LowPass, T, Sample, Warp>;
        template<typename T, typename Sample = T, Warp Warp = Warp::None>
        using BandPass = Details::StateSpace<Details::BandPass, T, Sample, Warp>;
        template<typename T, typename Sample = T, Warp Warp = Warp::None>
        using HighPass = Details::StateSpace<Details::HighPass, T, Sample, Warp>;
        template<typename T, typename Sample = T, Warp Warp = Warp::None>
        using Notch = Details::StateSpace<Details::Notch, T, Sample, Warp>;
        template<typename T, typename Sample = T, Warp Warp = Warp::None>
        using AllPass = Details::StateSpace<Details::AllPass, T, Sample, Warp>;
        template<typename T, typename Sample = T, Warp Warp = Warp::None>
        using LowShelf = Details::StateSpace<Details::LowShelf, T, Sample, Warp>;
        template<typename T, typename Sample = T, Warp Warp = Warp::None>
        using MidShelf = Details::StateSpace<Details::MidShelf, T, Sample, Warp>;
        template<typename T, typename Sample = T, Warp Warp = Warp::None>
        using HighShelf = Details::StateSpace<Details::HighShelf, T, Sample, Warp>;
    } // namespace Linear::SecondOrder
} // namespace ivantsov

#pragma once

#include <array>
#include <cmath>
#include <numbers>

namespace ivantsov
{
    template<typename T, typename S>
    concept is_algebraic = requires {
        { (((S {} * T {} + S {}) - S {}) + T {}) - T {} } -> std::same_as<S>;
    };

    enum struct Warp { None, Sigma };

    namespace Linear::FirstOrder::Details
    {
        enum struct Type { HighPass, LowPass, AllPass, HighShelf, LowShelf };

        template<Type Type, std::floating_point T, typename Sample, Warp Warp>
        requires is_algebraic<T, Sample>
        struct StateSpace
        {
            static constexpr bool is_shelf {Type == Type::HighShelf || Type == Type::LowShelf};

            auto initialize(const auto sample_rate) { fs = sample_rate; }

            auto processed(const auto x)
            {
                const auto theta {(x - z) * b[0]};
                auto y {theta * b[1]};
                if constexpr (Type != Type::HighPass) y = y + z;
                z = z + theta;
                if constexpr (Type == Type::LowShelf) return y * b[2];
                else return y;
            }

            auto reset() { z = {}; }

            auto set_fc(const auto x)
            {
                w = fs / (T {2} * std::numbers::pi_v<T> * x);
                if constexpr (Warp == Warp::Sigma)
                {
                    if (w > std::numbers::inv_pi_v<T>)
                    {
                        sigma = (T {0.519355108} * w * w - T {0.2251345831}) / (w * w - T {0.5514625429});
                    }
                }
                update_blocks();
            }

            auto set_g(const auto x) requires is_shelf
            {
                b[2] = x;
                update_blocks();
            }

            auto set_sigma(const auto x) requires(Warp == Warp::None)
            {
                sigma = x / (T {2} * std::numbers::pi_v<T>);
                update_blocks();
            }

        private:
            auto update_blocks()
            {
                const auto v {[this](const auto x) { return std::sqrt(x + sigma * sigma); }};
                if constexpr (!is_shelf) b[0] = T {1} / (T {0.5} + v(w * w));
                if constexpr (Type == Type::HighPass) b[1] = w;
                if constexpr (Type == Type::LowPass) b[1] = T {0.5} + sigma;
                if constexpr (Type == Type::AllPass) b[1] = T {0.5} - v(w * w);
                if constexpr (Type == Type::LowShelf)
                {
                    b[0] = T {1} / (T {0.5} + v(w * w * b[2]));
                    b[1] = T {0.5} + v(w * w / b[2]);
                }
                if constexpr (Type == Type::HighShelf)
                {
                    b[0] = T {1} / (T {0.5} + v(w * w / b[2]));
                    b[1] = T {0.5} + v(w * w * b[2]);
                }
            }

            T sigma {std::numbers::inv_pi_v<T>};
            T fs {};
            T w {};
            std::array<T, 2 + is_shelf> b {T {}, T {0.5} + sigma};
            Sample z {};
        };
    }// namespace Linear::FirstOrder::Details

    namespace Linear::SecondOrder::Details
    {
        enum struct Type { HighPass, BandPass, LowPass, AllPass, Notch, HighShelf, LowShelf, MidShelf };

        template<Type Type, std::floating_point T, typename Sample, Warp Warp>
        requires is_algebraic<T, Sample>
        struct StateSpace
        {
            static constexpr bool is_shelf {Type == Type::LowShelf || Type == Type::HighShelf ||
                                            Type == Type::MidShelf};

            auto initialize(const auto sample_rate) { fs = sample_rate; }

            auto processed(const auto x)
            {
                const auto theta {(x - z[0] - z[1] * b[1]) * b[0]};
                auto y {theta * b[3] + z[1] * b[2]};
                if constexpr (Type != Type::HighPass && Type != Type::BandPass) y = y + z[0];
                z = {z[0] + theta, Sample {} - z[1] - theta * b[1]};
                if constexpr (Type == Type::LowShelf) return y * b[4];
                else return y;
            }

            auto reset() { z = {}; }

            auto set_fc(const auto x)
            {
                w = fs / (std::numbers::sqrt2_v<T> * std::numbers::pi_v<T> * x);
                if constexpr (Warp == Warp::Sigma)
                {
                    if (w > T {1} / (std::numbers::sqrt2_v<T> * std::numbers::pi_v<T>) )
                    {
                        sigma = (T {0.7344790372} * w * w - T {0.3183883807}) / (w * w - T {0.5514625429});
                    }
                }
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
                sigma = x / (std::numbers::sqrt2_v<T> * std::numbers::pi_v<T>);
                update_blocks();
            }

        private:
            auto update_blocks()
            {
                const auto w_sq {w * w};
                const auto sigma_sq {sigma * sigma};
                const auto l_vk {[&](const auto x, const auto y) {
                    const auto t {x * (y + y - T {1})};
                    return std::pair {std::sqrt(x * x + T {2} * t * sigma_sq + sigma_sq * sigma_sq), t + sigma_sq};
                }};
                const auto [v, k] {[&] {
                    if constexpr (Type == Type::LowShelf) return l_vk(w_sq * std::sqrt(b[4]), zeta * zeta);
                    if constexpr (Type == Type::HighShelf) return l_vk(w_sq / std::sqrt(b[4]), zeta * zeta);
                    if constexpr (Type == Type::MidShelf) return l_vk(w_sq, zeta * zeta / b[4]);
                    else return l_vk(w_sq, zeta * zeta);
                }()};
                b[0] = T {1} / (v + std::sqrt(v + k) + T {0.5});
                b[1] = std::sqrt(v + v);
                if constexpr (Type == Type::HighPass)
                {
                    b[2] = T {2} * w_sq / b[1];
                    b[3] = w_sq;
                }
                if constexpr (Type == Type::BandPass)
                {
                    b[2] = T {4} * w * zeta * sigma / b[1];
                    b[3] = T {2} * w * zeta * (sigma + T {1} / std::numbers::sqrt2_v<T>);
                }
                if constexpr (Type == Type::LowPass)
                {
                    b[2] = T {2} * sigma_sq / b[1];
                    b[3] = sigma_sq + std::numbers::sqrt2_v<T> * sigma + T {0.5};
                }
                if constexpr (Type == Type::Notch)
                {
                    b[2] = T {2} * (w_sq - sigma_sq) / b[1];
                    b[3] = w_sq - sigma_sq + T {0.5};
                }
                if constexpr (Type == Type::AllPass)
                {
                    b[2] = b[1];
                    b[3] = v - std::sqrt(v + k) + T {0.5};
                }
                if constexpr (is_shelf)
                {
                    const auto [v_b, k_b] {[&] {
                        if constexpr (Type == Type::LowShelf) return l_vk(w_sq / std::sqrt(b[4]), zeta * zeta);
                        if constexpr (Type == Type::HighShelf) return l_vk(w_sq * std::sqrt(b[4]), zeta * zeta);
                        else return l_vk(w_sq, zeta * zeta * b[4]);
                    }()};
                    b[2] = T {2} * v_b / b[1];
                    b[3] = v_b + std::sqrt(v_b + k_b) + T {0.5};
                }
            }

            T fs {};
            T w {};
            T zeta {};
            T sigma {std::numbers::sqrt2_v<T> * std::numbers::inv_pi_v<T>};
            std::array<T, 4 + is_shelf> b {};
            std::array<Sample, 2> z {};
        };
    }// namespace Linear::SecondOrder::Details

    namespace Linear::FirstOrder
    {
        template<typename T, typename Sample = T, Warp Warp = Warp::None>
        using LowPass = Details::StateSpace<Details::Type::LowPass, T, Sample, Warp>;
        template<typename T, typename Sample = T, Warp Warp = Warp::None>
        using HighPass = Details::StateSpace<Details::Type::HighPass, T, Sample, Warp>;
        template<typename T, typename Sample = T, Warp Warp = Warp::None>
        using AllPass = Details::StateSpace<Details::Type::AllPass, T, Sample, Warp>;
        template<typename T, typename Sample = T, Warp Warp = Warp::None>
        using LowShelf = Details::StateSpace<Details::Type::LowShelf, T, Sample, Warp>;
        template<typename T, typename Sample = T, Warp Warp = Warp::None>
        using HighShelf = Details::StateSpace<Details::Type::HighShelf, T, Sample, Warp>;
    }// namespace Linear::FirstOrder

    namespace Linear::SecondOrder
    {
        template<typename T, typename Sample = T, Warp Warp = Warp::None>
        using LowPass = Details::StateSpace<Details::Type::LowPass, T, Sample, Warp>;
        template<typename T, typename Sample = T, Warp Warp = Warp::None>
        using BandPass = Details::StateSpace<Details::Type::BandPass, T, Sample, Warp>;
        template<typename T, typename Sample = T, Warp Warp = Warp::None>
        using HighPass = Details::StateSpace<Details::Type::HighPass, T, Sample, Warp>;
        template<typename T, typename Sample = T, Warp Warp = Warp::None>
        using Notch = Details::StateSpace<Details::Type::Notch, T, Sample, Warp>;
        template<typename T, typename Sample = T, Warp Warp = Warp::None>
        using AllPass = Details::StateSpace<Details::Type::AllPass, T, Sample, Warp>;
        template<typename T, typename Sample = T, Warp Warp = Warp::None>
        using LowShelf = Details::StateSpace<Details::Type::LowShelf, T, Sample, Warp>;
        template<typename T, typename Sample = T, Warp Warp = Warp::None>
        using MidShelf = Details::StateSpace<Details::Type::MidShelf, T, Sample, Warp>;
        template<typename T, typename Sample = T, Warp Warp = Warp::None>
        using HighShelf = Details::StateSpace<Details::Type::HighShelf, T, Sample, Warp>;
    }// namespace Linear::SecondOrder
}// namespace ivantsov

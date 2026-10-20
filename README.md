# ivantsov-filters

C++20 decramped IIR filters with configurable order and audio rate modulation.

## Example

Use scalar or SIMD sample types.

```cpp
#include "filters.h"
#include <cmath>
#include <numbers>
#include <xsimd/xsimd.hpp> // For the xsimd example.

auto main() -> int {
     using namespace ivantsov;
     using namespace ivantsov::Linear;
     using Sample = double; // Or xsimd::batch<double>.

     Filter<HighShelf, 2, Sample, Correction::Sigma> high_shelf {};
     high_shelf.initialise(44'100.0);
     high_shelf.set<Fc, Zeta, Gain>(1'000.0, 0.71, 2.0);

     const auto y {high_shelf.processed(Sample {1.0})};
     const auto magnitude {high_shelf.magnitude(std::cos(2 * std::numbers::pi_v<double> * 1'000.0 / 44'100.0))};

     high_shelf.reset();
}
```

## Parameters

`Fc` is the cutoff frequency in Hz. `Gain` is the shelf gain (2 ≈ +6 dB). For orders above 1, `Zeta` sets damping as `1 / (2 * Q)`.

`magnitude(x)` returns the magnitude response and takes `x = cos(2 * pi * f / fs)`, where `f` is the frequency in Hz and `fs` is the sample rate in Hz.

`BandPass`, `BandStop`, and `BandShelf` require an even order.

## References

Yuriy Ivantsov, [On the ideal bilinear and biquadratic digital filter](https://ivantsovy.com/research/paper1.pdf).

Yuriy Ivantsov, [On the state space of a linear digital filter](https://ivantsovy.com/research/paper2.pdf).

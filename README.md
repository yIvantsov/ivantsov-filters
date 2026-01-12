# ivantsov-filters

Fast first and second order decramped IIR filters suited for audio rate modulation. 
Requires C++20. 

## Example

```cpp
#include "filters.h"

template<std::floating_point T>
struct StereoSample
{
    T left;
    T right;

    // operator overloading required for the is_algebraic concept
    StereoSample operator+(const StereoSample x) const { return {left + x.left, right + x.right}; }
    StereoSample operator-(const StereoSample x) const { return {left - x.left, right - x.right}; }
    StereoSample operator+(const T x) const { return {left + x, right + x}; }
    StereoSample operator-(const T x) const { return {left - x, right - x}; }
    StereoSample operator*(const T x) const { return {left * x, right * x}; }
};

int main()
{
    using namespace ivantsov::Linear; 
    SecondOrder::HighShelf<float, StereoSample<float>, ivantsov::Warp::Sigma> high_shelf {};
    high_shelf.initialize(44'100.0f); // sets the sampling rate in Hz
    high_shelf.set_fc(100.0f); // sets the cutoff frequeny in Hz
    high_shelf.set_damping(0.71f); // sets the damping (1.0f / (2.0f * Q))
    high_shelf.set_g(2.0f); // sets the gain (2.0f == 6dB) 
    const auto [left, right] {high_shelf.processed(StereoSample {1.0f, 1.0f})};
    high_shelf.reset(); // resets the states
    
    // computes the magnitude response at 10Hz. 
    const auto magnitude {high_shelf.magnitude(std::cos(2.0 * pi * 10.0 / 44'100.0))}
}
```

## References
1) Yuriy Ivantsov, “On the ideal bilinear and biquadratic digital filter.” [Online]. Available: https://ivantsovy.com/research/paper1.pdf
2) Yuriy Ivantsov, “On the state space of a linear digital filter.” [Online]. Available: https://ivantsovy.com/research/paper2.pdf

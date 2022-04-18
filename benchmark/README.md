A couple of comments on the benchmark results:
- Starting from March, 20th. the `max_time` is 3600., `atol` is 10^-5 and `rtol` is 10^-6. Before it was 120., 10^-5 and 10^-5.
- April 15th, Uniform inexact Newton:
```
    cgatol = min(maxtol, ξ * gNorm2^(1 + ζ))
    cgatol = max(mintol, cgatol) # add some feasible limit
    cgrtol = min(maxtol, ξ * gNorm2^ζ)
    cgrtol = max(mintol, cgrtol) # add some feasible limit
```
- April 15th, use 2-norm in optimality check

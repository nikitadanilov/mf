Overview
--------

To everyone sharing von Neumann's passionate ha^W dislike for floating point
arithmetic and with a sincere apology to Proffs. W. Kahan and V. S. Ryabenkiy.

`rat.hpp` is a simple fixed point arithmetics implementation. Fixed point
rationals are a particular type of what is called _mixed fractions_ in
elementary schools.

`rat_t<E, F>` is a class of fixed point rationals, where the integer (aka
_entier_) part has type `E` and the fractional part has type `F`. An instance of
this class is a pair `(e, f)`, where `e` has type `E` and `f`---`F`. This pair
represents the rational number `e + f/N` where `N` is the (not actually
representable) value one greater than the maximal value representable by
`F`. For example

    rat_t<int64_t, uint64_t> epsilon(0, 1)
    
is `0 + 1/2^64 == 2^{-64} == 0.0000000000000000000542101086242752217003726400434970855712890625`.

The trick is to implement arithmetic operations for `rat_t<E, F>` without using
integer types wider than `E` and `F` (for one thing, they can already be the
widest types accessible on the architecture). This requires manual [long
division](https://www.cofault.com/2025/02/long-story-of-division.html) and
_half-width multiplication_.



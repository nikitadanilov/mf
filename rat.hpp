/* -*- C++ -*- */

/* Copyright 2025 Nikita Danilov <danilov@gmail.com> */

/* c++ -std=c++20 rat.cpp */

#pragma once

#include <stdexcept>
#include <string>
#include <type_traits>
#include <cstdint>
#include <limits>
#include <climits>
#include <cmath>

/*
 * Mixed fractions.
 *
 * Notation:
 *
 *     Z     - the set of integers
 *     Q     - the set of rationals
 *     x : A - x belongs to the set A
 *     S(n)  - the set { 0, 1, ..., n - 1 }
 *
 * For n : Z, n > 0, define mixed fractions MF(n) = Z * S(n), i.e., r : MF(n) iff r = (a, b), a : Z, b : S(n)
 *
 * Define v : MF(n) -> Q, v(a, b) = a + b/n.
 *
 * Proposition: v is an injection.
 *
 * Let v(a, b) = v(A, B).
 *
 *     a + b/n = A + B/n
 *
 *     n * (a - A) = B - b
 *
 *     As b and B vary from 0 to n - 1, the right hand side varies from -(n - 1)
 *     to n - 1, so is not a non-zero multiple of n.
 *
 *     Hence, the left hand side is not a non-zero multiple of n.
 *
 *     Hence n = 0, B = b, A = a.
 *
 * Define carry(x) = floor(x/n), mod(x) = x mod n.
 *
 *     Then by definition of mod, x = n * carry(x) + mod(x) (Note that n > 0.)
 *
 *     x/n = carry(x) + mod(x)/n = v(carry(x), mod(x))
 *
 * Define "normalisation" N : Z * Z -> MF(n), N(a, b) = (a + carry(b), mod(b)).
 *
 * Then for a, b : Z, a + b/n = a + carry(b) + mod(b)/n = v(N(a, b)).
 */
template <typename E, typename F> struct rat_t {
        using en_t = E; /* Entier (integral) part. */
        using fr_t = F; /* Fractional part. */
#define CHECK(field, msg) \
        static_assert((std::numeric_limits<en_t>::field) && \
                      (std::numeric_limits<fr_t>::field), "Entier and fractional types must " msg)
        static_assert(std::is_signed<en_t>{},        "Entier type must be signed.");
        static_assert(std::is_unsigned<fr_t>{},      "Fractional type must be unsigned.");
        CHECK(is_exact,                              "be exact.");
        CHECK(is_integer,                            "be integer.");
        CHECK(round_style == std::round_toward_zero, "round toward zero.");
        CHECK(radix == 2,                            "be binary.");
        en_t e;
        fr_t f;
        constexpr static double n        = (double)std::numeric_limits<F>::max() + 1;
        constexpr static int    bits     = std::numeric_limits<fr_t>::digits;
        constexpr static int    halfbits = bits / 2;
        static_assert(bits == halfbits * 2, "Fractional type needs even width.");
        constexpr static fr_t mask = ((fr_t)1 << halfbits) - 1;
        constexpr static fr_t lo(fr_t x) { return x & mask; };
        constexpr static fr_t hi(fr_t x) { return (x >> halfbits) & mask; };
        rat_t(en_t ee, fr_t ff) : e(ee), f(ff & (((fr_t)1 << bits) - 1)) {}
        rat_t(en_t ee) : rat_t(ee, 0) {}
        rat_t() : rat_t(0) {}
        static rat_t fromdouble(double x) {
                return rat_t(floor(x), (x - floor(x)) * n);
        }

        /* Returns (carry(a * b), mod(a * b)), without using wider types. */
        static rat_t mul(fr_t n0, fr_t n1) {
                rat_t result;
                /*
                 * Split n0 and n1 into low and high parts:
                 *
                 *   n0 = (h0 << halfbits) + l0
                 *   n1 = (h1 << halfbits) + l1
                 */
                fr_t l0 = lo(n0);
                fr_t h0 = hi(n0);
                fr_t l1 = lo(n1);
                fr_t h1 = hi(n1);

                /*
                 * n0 * n1 = ll
                 *         + (lh + hl) << halfbits
                 *         + hh        << bits
                 *
                 * where
                 */
                fr_t ll = l0 * l1;
                fr_t lh = l0 * h1;
                fr_t hl = h0 * l1;
                fr_t hh = h0 * h1;
                /*
                 * ll, lh, hl and hh all fit in fr_t.
                 *
                 * Example: bits = 32
                 *
                 * 6   6         5         4         3         2         1         0
                 * 43210987654321098765432109876543210987654321098765432109876543210
                 *                                 |---------- 1 << bits ----------|
                 *                                                 |------ l0 -----|
                 *                                                 |------ h0 -----|
                 *                                                 |------ l1 -----|
                 *                                                 |------ h0 -----|
                 *                                 |-------------- ll -------------|
                 *                 |-------- lh << halfbits -------|
                 *                 |-------- hl << halfbits -------|
                 * |---------- hh << bits ---------|
                 *
                 *                 |------ sum << halfbits --------|
                 * |---- carry ----| <.................................. (carry << 3 * halfbits)
                 *
                 * |--------------------------- result ----------------------------|
                 *                                 |----------- result.f ----------|
                 * |------ result.e << bits -------|
                 */

                /* First add lh + hl in fr_t, capturing any 1-bit overflow. */
                fr_t sum   = lh + hl;
                fr_t carry = sum < lh;

                /* Now add (ll >> halfbits) into sum, again track carry. */
                sum   += ll >> halfbits;
                carry += sum < (ll >> halfbits);

                /*
                 * The lower bits come from:
                 *
                 *  (ll & mask)  in the lowest halfbits
                 *  (sum & mask) in the next halfbits
                 */
                result.f = (ll & mask) | ((sum & mask) << halfbits);

                /*
                 * The high bits come from:
                 *
                 *   hh + (sum >> halfbits) + (carry << halfbits)
                 *
                 * This adds in the top half of sum plus whatever carry is left.
                 */
                result.e = hh + (sum >> halfbits) + (carry << halfbits);
                return result;
        }
        static fr_t topbits(fr_t x, int nr) {
                return x >> (bits - nr);
        }
        rat_t operator <<(int nr) const {
                return rat_t((e << nr) | topbits(f, nr), f << nr);
        }
        void set(int nr) {
                if (nr >= bits) {
                        e |= (fr_t)1 << (nr - bits);
                } else {
                        f |= (fr_t)1 << nr;
                }
        }
        /*
         * Calculate 1/v(a, b).
         *
         *     1/v(a, b) = 1/(a + b/n) = n/(n*a + b) = 0 + (n^2/(n*a + b))/n =
         *
         *     v(N(0, n^2/(n*a + b))
         *
         * Calculate n^2/(n*a + b) via long division.
         */
        rat_t inverse() const {
                /* Long division: (1 << (2 * bits)) / (*this) */
                rat_t q;
                rat_t r;
                for (int i = 2 * bits; i >= 0; i--) {
                        r = r << 1;
                        r.f |= (i == 2 * bits);
                        if (r >= *this) {
                                r = r - *this;
                                q.set(i);
                        }
                }
                return q;
        }
        /*
         * Calculate v(a, b) + v(A, B).
         *
         *     v(a, b) + v(A, B) = a + A + (b + B)/n = v(N(a + A, b + B)) =
         *
         *     v(a + A + carry(b + B), mod(b + B))
         */
        rat_t operator +(const rat_t &o) const {
                return rat_t(e + o.e + (f > (fr_t)(f + o.f)), (fr_t)f + o.f);
        }
        /*
         * Calculate -v(a, b).
         *
         *     -v(a, b) = -a - b/n = v(N(-a, -b)), or
         *
         *     if b = 0, then -v(a, 0) = -a = v(-a, 0)
         *
         *     if b > 0, then -v(a, b) = -a - b/n = -a - 1 + n/n - b/n = -a - 1 + (n - b)/n =
         *
         *                     v(-a - 1, n - b), as n - b : S(n)
         */
        rat_t operator -() const {
                return f == 0 ? rat_t(-e, 0) : rat_t(-e - 1, ~f);
        }
        rat_t operator -(const rat_t &o) const {
                return *this + -o;
        }

        static inline rat_t zero = rat_t();

        rat_t rabs() const {
                return *this > zero ? *this : -*this;
        }
        /*
         * Calculate v(a, b) * v(A, B).
         *
         *     v(a, b) * v(A, B) = (a + b/n) * (A + B/n) =
         *
         *     a*A + (a*B + b*A + b*B/n)/n = v(N(a*A, a*B + b*A + carry(b * B))) + mod(b * B)/n^2
         */
        rat_t operator *(const rat_t &o) const {
                return  rat_t(abs(e) * abs(o.e)).sgn(e * o.e) + /* a * A */
                        mul(abs(e), o.f).sgn(e) +               /* a * B / n */
                        mul(f, abs(o.e)).sgn(o.e) +             /* b * A / n */
                        rat_t(0, mul(f, o.f).e);                /* b * B / n^2 */
        }
        rat_t operator /(const rat_t &o) const {
                if (o.e == 0 && o.f == 0)
                        throw std::domain_error("Division by 0.");
                else if (o.e == 0 && o.f == 1)
                        throw std::overflow_error("Division result is too large.");
                else if (o.e == 1 && o.f == 1)
                        return *this;
                else
                        return *this * o.rabs().inverse().sgn(o.e);
        }
        double todouble() const {
                return (1.0 * f) / n + e;
        }
        explicit operator double() const {
                return todouble();
        }
        std::string tostring() const {
                return std::string("[") + std::to_string(e) + ", " + std::to_string(f) + "]";
        }
        operator std::string () const {
                return tostring();
        }
        std::strong_ordering operator <=>(const rat_t &o) const {
                std::strong_ordering ecmp = e <=> o.e;
                return ecmp != 0 ? ecmp : f <=> o.f;
        }
        rat_t sgn(en_t x) const {
                return x < 0 ? -(*this) : *this;
        }
};

template <typename E, typename F, std::uintmax_t M = std::numeric_limits<F>::max()>
std::ostream &operator <<(std::ostream &s, const rat_t<E, F> &r) {
        return s << r.todouble() << ": " << r.tostring();
}

/*
 *  Local variables:
 *  c-indentation-style: "K&R"
 *  c-basic-offset: 8
 *  tab-width: 8
 *  scroll-step: 1
 *  indent-tabs-mode: nil
 *  End:
 */

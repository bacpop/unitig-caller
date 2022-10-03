// Copyright (c) 2016, the SDSL Project Authors.  All rights reserved.
// Please see the AUTHORS file for details.  Use of this source code is governed
// by a BSD license that can be found in the LICENSE file.
/*!\file k2_treap_helper.hpp
 * \brief k2_treap_helper.hpp contains helper functions and definitions for a k^2-treap implementation.
 * \author Simon Gog
 */
#ifndef INCLUDED_SDSL_K2_TREAP_HELPER
#define INCLUDED_SDSL_K2_TREAP_HELPER

#include <algorithm>
#include <array>
#include <complex>
#include <iterator>
#include <queue>
#include <tuple>
#include <vector>

#include <sdsl/bits.hpp>
#include <sdsl/vectors.hpp>

//! Namespace for the succinct data structure library.
namespace sdsl
{

namespace k2_treap_ns
{

// Precomputed value for fast k^2 treap operations
template <uint8_t t_k>
struct precomp
{
    static struct impl
    {
        uint64_t exp[65];
        impl()
        {
            exp[0] = 1;
            for (uint8_t i = 1; i < 65; ++i) { exp[i] = t_k * exp[i - 1]; }
        }
    } data;

    static uint64_t exp(uint8_t l) { return data.exp[l]; }

    static uint64_t divexp(uint64_t x, uint8_t l) { return x / data.exp[l]; }

    static uint64_t modexp(uint64_t x, uint8_t l) { return x % data.exp[l]; }
};

template <>
struct precomp<2>
{
    static uint64_t exp(uint8_t l) { return 1ULL << l; }

    static uint64_t divexp(uint64_t x, uint8_t l) { return x >> l; }

    static uint64_t modexp(uint64_t x, uint8_t l) { return x & bits::lo_set[l]; }
};

template <>
struct precomp<4>
{
    static uint64_t exp(uint8_t l) { return 1ULL << (2 * l); }

    static uint64_t divexp(uint64_t x, uint8_t l) { return x >> (2 * l); }

    static uint64_t modexp(uint64_t x, uint8_t l) { return x & bits::lo_set[2 * l]; }
};

template <>
struct precomp<8>
{
    static uint64_t exp(uint8_t l) { return 1ULL << (3 * l); }

    static uint64_t divexp(uint64_t x, uint8_t l) { return x >> (3 * l); }

    static uint64_t modexp(uint64_t x, uint8_t l) { return x & bits::lo_set[3 * l]; }
};

template <>
struct precomp<16>
{
    static uint64_t exp(uint8_t l) { return 1ULL << (4 * l); }

    static uint64_t divexp(uint64_t x, uint8_t l) { return x >> (4 * l); }

    static uint64_t modexp(uint64_t x, uint8_t l) { return x & bits::lo_set[4 * l]; }
};

template <uint8_t t_k>
typename precomp<t_k>::impl precomp<t_k>::data;

typedef std::complex<uint64_t> t_p;
typedef t_p point_type;
typedef t_p range_type;

struct node_type
{
    uint8_t t;      // level; size of node 1<<t
    t_p p;          // lower left corner
    uint64_t idx;   // index in bp
    uint64_t max_v; // maximal value
    t_p max_p;      // maximal point

    node_type() = default;
    node_type(uint8_t _t, t_p _p, uint64_t _idx, uint64_t _max_v, t_p _max_p)
      : t(_t)
      , p(_p)
      , idx(_idx)
      , max_v(_max_v)
      , max_p(_max_p)
    {}
    node_type(node_type &&) = default;
    node_type(const node_type &) = default;
    node_type & operator=(node_type &&) = default;
    node_type & operator=(const node_type &) = default;

    bool operator<(const node_type & v) const
    {
        if (max_v != v.max_v) { return max_v < v.max_v; }
        if (real(max_p) != real(v.max_p)) { return real(max_p) > real(v.max_p); }
        return imag(max_p) > imag(v.max_p);
    }
};

} // namespace k2_treap_ns

} // namespace sdsl
#endif

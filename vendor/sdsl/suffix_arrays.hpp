// Copyright (c) 2016, the SDSL Project Authors.  All rights reserved.
// Please see the AUTHORS file for details.  Use of this source code is governed
// by a BSD license that can be found in the LICENSE file.
/*!\file suffix_arrays.hpp
 * \brief suffix_arrays.hpp contains generic classes for different suffix array classes.
 * \author Simon Gog
 */
#ifndef INCLUDED_SDSL_SUFFIX_ARRAYS
#define INCLUDED_SDSL_SUFFIX_ARRAYS

#include <sdsl/sdsl_concepts.hpp>

/** \defgroup csa Compressed Suffix Arrays (CSA) */

#include <sdsl/construct.hpp>
#include <sdsl/csa_bitcompressed.hpp>
#include <sdsl/csa_sada.hpp>
#include <sdsl/csa_wt.hpp>
#include <sdsl/suffix_array_algorithm.hpp>
#include <sdsl/wavelet_trees.hpp>

namespace sdsl
{

//! Typedef for convenient usage of std integer alphabet strategy
template <class t_wt = wt_int<>,
          uint32_t t_dens = 32,
          uint32_t t_inv_dens = 64,
          class t_sa_sample_strat = sa_order_sa_sampling<>,
          class t_isa_sample_strat = isa_sampling<>>
using csa_wt_int = csa_wt<t_wt, t_dens, t_inv_dens, t_sa_sample_strat, t_isa_sample_strat, int_alphabet<>>;

template <class t_enc_vec = enc_vector<>,                   // Vector type used to store the Psi-function
          uint32_t t_dens = 32,                             // Sample density for suffix array (SA) values
          uint32_t t_inv_dens = 64,                         // Sample density for inverse suffix array (ISA) values
          class t_sa_sample_strat = sa_order_sa_sampling<>, // Policy class for the SA sampling. Alternative
                                                            // text_order_sa_sampling.
          class t_isa_sample_strat = isa_sampling<>         // Policy class for the ISA sampling.
          >
using csa_sada_int = csa_sada<t_enc_vec, t_dens, t_inv_dens, t_sa_sample_strat, t_isa_sample_strat, int_alphabet<>>;
} // namespace sdsl

#endif

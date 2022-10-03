// Copyright (c) 2016, the SDSL Project Authors.  All rights reserved.
// Please see the AUTHORS file for details.  Use of this source code is governed
// by a BSD license that can be found in the LICENSE file.
/*!\file coder.hpp
 * \brief coder.hpp contains the coder namespace and includes the header files of sdsl::coder::fibonacci,
 * sdsl::coder::elias_delta, and sdsl::coder::run_length \author Simon Gog
 */
#ifndef SDSL_CODER
#define SDSL_CODER

#include <sdsl/coder_comma.hpp>
#include <sdsl/coder_elias_delta.hpp>
#include <sdsl/coder_elias_gamma.hpp>
#include <sdsl/coder_fibonacci.hpp>
#include <sdsl/int_vector.hpp>

namespace sdsl
{

//! Namespace for the different coder of the sdsl.
namespace coder
{

template <class Coder>
class run_length
{
  public:
    typedef uint64_t size_type;
    static void encode(uint64_t x, uint64_t *& z, uint8_t offset);
    static uint64_t encoding_length(const uint64_t * s, uint8_t s_offset, size_type bit_length);
};

template <class Coder>
typename run_length<Coder>::size_type run_length<Coder>::encoding_length(const uint64_t * s,
                                                                         uint8_t s_offset,
                                                                         size_type bit_length)
{
    assert(s_offset < 64);
    size_type i = 0;
    uint64_t w = (*s >> s_offset);
    uint8_t last_bit = w & 1;
    size_type result = 0;
    while (i < bit_length)
    {
        size_type len = 0;
        while (last_bit == (w & 1) and i < bit_length)
        {
            //			std::cout<<w<<" "<<i<<std::endl;
            ++len;
            ++i;
            ++s_offset;
            w >>= 1;
            if (s_offset == 64)
            {
                s_offset = 0;
                w = *(++s);
            }
        }
        //		std::cout<<"len="<<Coder::encoding_length(len)<<std::endl;
        last_bit = (w & 1);
        result += Coder::encoding_length(len);
    }
    return result;
}

} // end namespace coder

} // end namespace sdsl

#endif

#pragma once

#include "maj_bitvec.hpp"
#include "maj_popcount.hpp"

#include <cassert>
#include <vector>

namespace maj
{

struct approx_block3_stats
{
  int group_size{ 0 };
  int compressed_size{ 0 };
  int groups3{ 0 };
  int passthrough_bits{ 0 };
  int popcount_levels{ 0 };
};

namespace detail
{

template<typename Ntk>
signal_t<Ntk> majority3( Ntk& ntk, signal_t<Ntk> a, signal_t<Ntk> b, signal_t<Ntk> c )
{
  auto const ab = ntk.create_and( a, b );
  auto const ac = ntk.create_and( a, c );
  auto const bc = ntk.create_and( b, c );
  return make_or( ntk, make_or( ntk, ab, ac ), bc );
}

} // namespace detail

template<typename Ntk>
signal_t<Ntk> create_maj_approx_block3_popcount( Ntk& ntk, std::vector<signal_t<Ntk>> const& xs, approx_block3_stats* st = nullptr )
{
  assert( !xs.empty() );
  assert( ( xs.size() & 1u ) == 1u );

  approx_block3_stats local_stats;
  approx_block3_stats* stats = ( st == nullptr ) ? &local_stats : st;
  *stats = {};
  stats->group_size = 3;

  int const groups = static_cast<int>( xs.size() / 3u );
  int const rem = static_cast<int>( xs.size() % 3u );

  std::vector<signal_t<Ntk>> compressed;
  compressed.reserve( static_cast<size_t>( groups + rem ) );

  for ( int i = 0; i < groups; ++i )
  {
    int const b = 3 * i;
    compressed.emplace_back(
        detail::majority3( ntk, xs[static_cast<size_t>( b )], xs[static_cast<size_t>( b + 1 )], xs[static_cast<size_t>( b + 2 )] ) );
  }
  for ( int i = groups * 3; i < static_cast<int>( xs.size() ); ++i )
  {
    compressed.emplace_back( xs[static_cast<size_t>( i )] );
  }

  popcount_stats pop_stats{};
  auto const y = create_maj_popcount_exact( ntk, compressed, &pop_stats );

  stats->groups3 = groups;
  stats->passthrough_bits = rem;
  stats->compressed_size = static_cast<int>( compressed.size() );
  stats->popcount_levels = pop_stats.popcount_levels;
  return y;
}

template<typename Ntk>
signal_t<Ntk> create_maj_approx_block3( Ntk& ntk, std::vector<signal_t<Ntk>> const& xs, approx_block3_stats* st = nullptr )
{
  return create_maj_approx_block3_popcount( ntk, xs, st );
}

} // namespace maj

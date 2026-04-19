#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <utility>
#include <vector>

namespace maj
{

template<typename Ntk>
using signal_t = typename Ntk::signal;

inline int width_for_size( int size )
{
  assert( size >= 0 );
  int w = 0;
  int limit = 1;
  while ( limit < size + 1 )
  {
    limit <<= 1;
    ++w;
  }
  return w;
}

template<typename Ntk>
signal_t<Ntk> make_or( Ntk& ntk, signal_t<Ntk> a, signal_t<Ntk> b )
{
  return !ntk.create_and( !a, !b );
}

template<typename Ntk>
signal_t<Ntk> eq_bit( Ntk& ntk, signal_t<Ntk> a, signal_t<Ntk> b )
{
  return !ntk.create_xor( a, b );
}

template<typename Ntk>
signal_t<Ntk> mux_bit( Ntk& ntk, signal_t<Ntk> sel, signal_t<Ntk> t, signal_t<Ntk> f )
{
  auto const s_t = ntk.create_and( sel, t );
  auto const ns_f = ntk.create_and( !sel, f );
  return make_or( ntk, s_t, ns_f );
}

template<typename Ntk>
std::vector<signal_t<Ntk>> make_const_vec( Ntk& ntk, int width, uint64_t value )
{
  assert( width >= 0 );
  std::vector<signal_t<Ntk>> out;
  out.reserve( static_cast<size_t>( width ) );
  for ( int i = 0; i < width; ++i )
  {
    out.emplace_back( ntk.get_constant( ( ( value >> i ) & 1ULL ) != 0ULL ) );
  }
  return out;
}

template<typename Ntk>
std::vector<signal_t<Ntk>> zero_extend( Ntk& ntk, std::vector<signal_t<Ntk>> const& v, int new_w )
{
  assert( new_w >= static_cast<int>( v.size() ) );
  std::vector<signal_t<Ntk>> out = v;
  out.resize( static_cast<size_t>( new_w ), ntk.get_constant( false ) );
  return out;
}

template<typename Ntk>
signal_t<Ntk> is_zero( Ntk& ntk, std::vector<signal_t<Ntk>> const& v )
{
  auto any = ntk.get_constant( false );
  for ( auto const& bit : v )
  {
    any = make_or( ntk, any, bit );
  }
  return !any;
}

template<typename Ntk>
signal_t<Ntk> eq_vec( Ntk& ntk, std::vector<signal_t<Ntk>> const& a, std::vector<signal_t<Ntk>> const& b )
{
  assert( a.size() == b.size() );
  auto eq = ntk.get_constant( true );
  for ( size_t i = 0; i < a.size(); ++i )
  {
    eq = ntk.create_and( eq, eq_bit( ntk, a[i], b[i] ) );
  }
  return eq;
}

template<typename Ntk>
signal_t<Ntk> gt_unsigned( Ntk& ntk, std::vector<signal_t<Ntk>> const& a, std::vector<signal_t<Ntk>> const& b )
{
  assert( a.size() == b.size() );
  auto gt = ntk.get_constant( false );
  auto eq = ntk.get_constant( true );

  for ( int i = static_cast<int>( a.size() ) - 1; i >= 0; --i )
  {
    auto const ai_gt_bi = ntk.create_and( a[static_cast<size_t>( i )], !b[static_cast<size_t>( i )] );
    gt = make_or( ntk, gt, ntk.create_and( eq, ai_gt_bi ) );
    eq = ntk.create_and( eq, eq_bit( ntk, a[static_cast<size_t>( i )], b[static_cast<size_t>( i )] ) );
  }

  return gt;
}

template<typename Ntk>
std::vector<signal_t<Ntk>> add_unsigned( Ntk& ntk, std::vector<signal_t<Ntk>> const& a, std::vector<signal_t<Ntk>> const& b )
{
  assert( a.size() == b.size() );
  std::vector<signal_t<Ntk>> sum( a.size(), ntk.get_constant( false ) );
  auto carry = ntk.get_constant( false );

  for ( size_t i = 0; i < a.size(); ++i )
  {
    auto const axb = ntk.create_xor( a[i], b[i] );
    sum[i] = ntk.create_xor( axb, carry );

    auto const ab = ntk.create_and( a[i], b[i] );
    auto const carry_axb = ntk.create_and( carry, axb );
    carry = make_or( ntk, ab, carry_axb );
  }

  return sum;
}

template<typename Ntk>
std::vector<signal_t<Ntk>> sub_unsigned( Ntk& ntk, std::vector<signal_t<Ntk>> const& a, std::vector<signal_t<Ntk>> const& b )
{
  assert( a.size() == b.size() );
  std::vector<signal_t<Ntk>> diff( a.size(), ntk.get_constant( false ) );
  auto borrow = ntk.get_constant( false );

  for ( size_t i = 0; i < a.size(); ++i )
  {
    auto const axb = ntk.create_xor( a[i], b[i] );
    diff[i] = ntk.create_xor( axb, borrow );

    auto const not_a_and_b = ntk.create_and( !a[i], b[i] );
    auto const equal_ab = eq_bit( ntk, a[i], b[i] );
    auto const keep_borrow = ntk.create_and( equal_ab, borrow );
    borrow = make_or( ntk, not_a_and_b, keep_borrow );
  }

  return diff;
}

template<typename Ntk>
std::pair<std::vector<signal_t<Ntk>>, signal_t<Ntk>>
sub_unsigned_with_borrow( Ntk& ntk,
                          std::vector<signal_t<Ntk>> const& a,
                          std::vector<signal_t<Ntk>> const& b )
{
  assert( a.size() == b.size() );
  std::vector<signal_t<Ntk>> diff( a.size(), ntk.get_constant( false ) );
  auto borrow = ntk.get_constant( false );

  for ( size_t i = 0; i < a.size(); ++i )
  {
    auto const axb = ntk.create_xor( a[i], b[i] );
    diff[i] = ntk.create_xor( axb, borrow );

    auto const not_a_and_b = ntk.create_and( !a[i], b[i] );
    auto const equal_ab = eq_bit( ntk, a[i], b[i] );
    auto const keep_borrow = ntk.create_and( equal_ab, borrow );
    borrow = make_or( ntk, not_a_and_b, keep_borrow );
  }

  return { diff, borrow };
}

template<typename Ntk>
std::vector<signal_t<Ntk>> twos_complement( Ntk& ntk, std::vector<signal_t<Ntk>> const& v )
{
  std::vector<signal_t<Ntk>> inv;
  inv.reserve( v.size() );
  for ( auto const& bit : v )
  {
    inv.emplace_back( !bit );
  }

  auto const one = make_const_vec( ntk, static_cast<int>( v.size() ), 1u );
  return add_unsigned( ntk, inv, one );
}

template<typename Ntk>
std::vector<signal_t<Ntk>> mux_vec( Ntk& ntk, signal_t<Ntk> sel, std::vector<signal_t<Ntk>> const& t, std::vector<signal_t<Ntk>> const& f )
{
  assert( t.size() == f.size() );
  std::vector<signal_t<Ntk>> out;
  out.reserve( t.size() );
  for ( size_t i = 0; i < t.size(); ++i )
  {
    out.emplace_back( mux_bit( ntk, sel, t[i], f[i] ) );
  }
  return out;
}

template<typename Ntk>
int count_complemented_edges( Ntk const& ntk )
{
  int inv_count = 0;
  ntk.foreach_gate( [&]( auto const& n ) {
    ntk.foreach_fanin( n, [&]( auto const& fi ) {
      if ( ntk.is_complemented( fi ) )
      {
        ++inv_count;
      }
    } );
  } );
  ntk.foreach_po( [&]( auto const& f ) {
    if ( ntk.is_complemented( f ) )
    {
      ++inv_count;
    }
  } );
  return inv_count;
}

template<typename Ntk>
std::pair<int, int> count_and_xor_gates( Ntk const& ntk )
{
  int and_count = 0;
  int xor_count = 0;

  ntk.foreach_gate( [&]( auto const& n ) {
    if ( ntk.is_and( n ) )
    {
      ++and_count;
    }
    else if ( ntk.is_xor( n ) )
    {
      ++xor_count;
    }
  } );

  return { and_count, xor_count };
}

} // namespace maj

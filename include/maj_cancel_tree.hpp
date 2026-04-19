#pragma once

#include "maj_bitvec.hpp"

#include <algorithm>
#include <cassert>
#include <functional>
#include <vector>

namespace maj
{

struct cancel_tree_stats
{
  int max_k_width{ 0 };
  int num_merges{ 0 };
  int depth_merges{ 0 };
};

template<typename Ntk>
struct summary_t
{
  signal_t<Ntk> cand;
  std::vector<signal_t<Ntk>> k; // LSB-first
  int size{ 0 };
  int depth{ 0 };
};

namespace detail
{

template<typename Ntk>
summary_t<Ntk> make_leaf_summary( Ntk& ntk, signal_t<Ntk> x )
{
  summary_t<Ntk> out;
  out.cand = x;
  out.k = make_const_vec( ntk, width_for_size( 1 ), 1u );
  out.size = 1;
  out.depth = 0;
  return out;
}

template<typename Ntk>
summary_t<Ntk> combine_summaries( Ntk& ntk, summary_t<Ntk> const& a, summary_t<Ntk> const& b, cancel_tree_stats* st )
{
  assert( a.size > 0 );
  assert( b.size > 0 );

  int const out_size = a.size + b.size;
  int const w_out = width_for_size( out_size );

  if ( st != nullptr )
  {
    st->max_k_width = std::max( st->max_k_width, w_out );
  }

  auto const k_a_ext = zero_extend( ntk, a.k, w_out );
  auto const k_b_ext = zero_extend( ntk, b.k, w_out );

  auto const za = is_zero( ntk, a.k );
  auto const zb = is_zero( ntk, b.k );

  auto const same = eq_bit( ntk, a.cand, b.cand );
  auto const sum = add_unsigned( ntk, k_a_ext, k_b_ext );

  auto const gt = gt_unsigned( ntk, k_a_ext, k_b_ext );
  auto const eq = eq_vec( ntk, k_a_ext, k_b_ext );

  auto const diff_ab = sub_unsigned( ntk, k_a_ext, k_b_ext );
  auto const diff_ba = sub_unsigned( ntk, k_b_ext, k_a_ext );

  auto const zeros = make_const_vec( ntk, w_out, 0u );
  auto const k_diff_b = mux_vec( ntk, eq, zeros, diff_ba );
  auto const k_diff = mux_vec( ntk, gt, diff_ab, k_diff_b );

  auto const cand_diff = mux_bit( ntk, gt, a.cand, b.cand );

  auto const cand_core = mux_bit( ntk, same, a.cand, cand_diff );
  auto const k_core = mux_vec( ntk, same, sum, k_diff );

  auto const cand_no_za = mux_bit( ntk, zb, a.cand, cand_core );
  auto const k_no_za = mux_vec( ntk, zb, k_a_ext, k_core );

  auto const cand_raw = mux_bit( ntk, za, b.cand, cand_no_za );
  auto const k_raw = mux_vec( ntk, za, k_b_ext, k_no_za );

  auto const z_out = is_zero( ntk, k_raw );
  auto const cand_final = mux_bit( ntk, z_out, ntk.get_constant( false ), cand_raw );

  summary_t<Ntk> out;
  out.cand = cand_final;
  out.k = k_raw;
  out.size = out_size;
  out.depth = 1 + std::max( a.depth, b.depth );

  if ( st != nullptr )
  {
    ++st->num_merges;
    st->depth_merges = std::max( st->depth_merges, out.depth );
  }

  return out;
}

template<typename Ntk>
summary_t<Ntk> combine_summaries_v2( Ntk& ntk, summary_t<Ntk> const& a, summary_t<Ntk> const& b, cancel_tree_stats* st )
{
  assert( a.size > 0 );
  assert( b.size > 0 );

  int const out_size = a.size + b.size;
  int const w_out = width_for_size( out_size );

  if ( st != nullptr )
  {
    st->max_k_width = std::max( st->max_k_width, w_out );
  }

  auto const k_a_ext = zero_extend( ntk, a.k, w_out );
  auto const k_b_ext = zero_extend( ntk, b.k, w_out );

  auto const za = is_zero( ntk, a.k );
  auto const zb = is_zero( ntk, b.k );

  auto const same = eq_bit( ntk, a.cand, b.cand );
  auto const sum = add_unsigned( ntk, k_a_ext, k_b_ext );

  auto const [diff_ab, borrow] = sub_unsigned_with_borrow( ntk, k_a_ext, k_b_ext );
  auto const diff_ba_mag = twos_complement( ntk, diff_ab );
  auto const mag = mux_vec( ntk, borrow, diff_ba_mag, diff_ab );
  auto const cand_diff = mux_bit( ntk, borrow, b.cand, a.cand );

  auto const cand_core = mux_bit( ntk, same, a.cand, cand_diff );
  auto const k_core = mux_vec( ntk, same, sum, mag );

  auto const cand_no_za = mux_bit( ntk, zb, a.cand, cand_core );
  auto const k_no_za = mux_vec( ntk, zb, k_a_ext, k_core );

  auto const cand_raw = mux_bit( ntk, za, b.cand, cand_no_za );
  auto const k_raw = mux_vec( ntk, za, k_b_ext, k_no_za );

  auto const z_out = is_zero( ntk, k_raw );
  auto const cand_final = mux_bit( ntk, z_out, ntk.get_constant( false ), cand_raw );

  summary_t<Ntk> out;
  out.cand = cand_final;
  out.k = k_raw;
  out.size = out_size;
  out.depth = 1 + std::max( a.depth, b.depth );

  if ( st != nullptr )
  {
    ++st->num_merges;
    st->depth_merges = std::max( st->depth_merges, out.depth );
  }

  return out;
}

template<typename Ntk>
summary_t<Ntk> build_range( Ntk& ntk,
                            std::vector<signal_t<Ntk>> const& xs,
                            int l,
                            int r,
                            cancel_tree_stats* st )
{
  assert( l >= 0 );
  assert( l < r );
  assert( r <= static_cast<int>( xs.size() ) );

  if ( r - l == 1 )
  {
    return make_leaf_summary( ntk, xs[static_cast<size_t>( l )] );
  }

  int const mid = l + ( r - l ) / 2;
  auto const left = build_range( ntk, xs, l, mid, st );
  auto const right = build_range( ntk, xs, mid, r, st );
  return combine_summaries( ntk, left, right, st );
}

template<typename Ntk>
summary_t<Ntk> build_range_v2( Ntk& ntk,
                               std::vector<signal_t<Ntk>> const& xs,
                               int l,
                               int r,
                               cancel_tree_stats* st )
{
  assert( l >= 0 );
  assert( l < r );
  assert( r <= static_cast<int>( xs.size() ) );

  if ( r - l == 1 )
  {
    return make_leaf_summary( ntk, xs[static_cast<size_t>( l )] );
  }

  int const mid = l + ( r - l ) / 2;
  auto const left = build_range_v2( ntk, xs, l, mid, st );
  auto const right = build_range_v2( ntk, xs, mid, r, st );
  return combine_summaries_v2( ntk, left, right, st );
}

} // namespace detail

template<typename Ntk>
signal_t<Ntk> create_maj_cancel_tree( Ntk& ntk, std::vector<signal_t<Ntk>> const& xs, cancel_tree_stats* st = nullptr )
{
  assert( !xs.empty() );
  assert( ( xs.size() & 1u ) == 1u );

  cancel_tree_stats local_stats;
  cancel_tree_stats* stats = ( st == nullptr ) ? &local_stats : st;
  *stats = {};

  auto const root = detail::build_range( ntk, xs, 0, static_cast<int>( xs.size() ), stats );

  stats->max_k_width = std::max( stats->max_k_width, static_cast<int>( root.k.size() ) );
  stats->depth_merges = std::max( stats->depth_merges, root.depth );

  auto const k_nonzero = !is_zero( ntk, root.k );
  return mux_bit( ntk, k_nonzero, root.cand, ntk.get_constant( false ) );
}

template<typename Ntk>
signal_t<Ntk> create_boyermoore_tree( Ntk& ntk, std::vector<signal_t<Ntk>> const& xs, cancel_tree_stats* st = nullptr )
{
  return create_maj_cancel_tree( ntk, xs, st );
}

template<typename Ntk>
signal_t<Ntk> create_maj_cancel_tree_v2( Ntk& ntk, std::vector<signal_t<Ntk>> const& xs, cancel_tree_stats* st = nullptr )
{
  assert( !xs.empty() );
  assert( ( xs.size() & 1u ) == 1u );

  cancel_tree_stats local_stats;
  cancel_tree_stats* stats = ( st == nullptr ) ? &local_stats : st;
  *stats = {};

  auto const root = detail::build_range_v2( ntk, xs, 0, static_cast<int>( xs.size() ), stats );

  stats->max_k_width = std::max( stats->max_k_width, static_cast<int>( root.k.size() ) );
  stats->depth_merges = std::max( stats->depth_merges, root.depth );

  auto const k_nonzero = !is_zero( ntk, root.k );
  return mux_bit( ntk, k_nonzero, root.cand, ntk.get_constant( false ) );
}

template<typename Ntk>
signal_t<Ntk> create_boyermoore_tree_v2( Ntk& ntk, std::vector<signal_t<Ntk>> const& xs, cancel_tree_stats* st = nullptr )
{
  return create_maj_cancel_tree_v2( ntk, xs, st );
}

} // namespace maj

#pragma once

#include "maj_bitvec.hpp"

#include <algorithm>
#include <cassert>
#include <deque>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace maj
{

struct popcount_stats
{
  int popcount_levels{ 0 };
  int max_count_width{ 0 };
};

enum class csa_schedule_mode
{
  serial,
  wallace,
  dadda
};

inline csa_schedule_mode parse_csa_schedule_mode( std::string const& mode )
{
  if ( mode == "serial" )
  {
    return csa_schedule_mode::serial;
  }
  if ( mode == "wallace" )
  {
    return csa_schedule_mode::wallace;
  }
  if ( mode == "dadda" )
  {
    return csa_schedule_mode::dadda;
  }
  throw std::invalid_argument( "unsupported CSA schedule mode: " + mode );
}

inline std::string to_string( csa_schedule_mode mode )
{
  if ( mode == csa_schedule_mode::serial )
  {
    return "serial";
  }
  if ( mode == csa_schedule_mode::wallace )
  {
    return "wallace";
  }
  return "dadda";
}

struct popcount_strict_stats
{
  int scaffold_p{ 0 };
  int scaffold_inputs{ 0 };
  int scaffold_threshold{ 0 };
  int comparator_width{ 0 };
  int num_fixed_pairs{ 0 };
  int csa_fa_count{ 0 };
  int comparator_fa_count{ 0 };
  int total_fa_count{ 0 };
  int csa_levels{ 0 };
  int total_levels{ 0 };
  csa_schedule_mode schedule_mode{ csa_schedule_mode::dadda };
};

template<typename Ntk>
struct count_summary_t
{
  std::vector<signal_t<Ntk>> bits;
  int size{ 0 };
  int depth{ 0 };
};

namespace detail
{

inline std::vector<int> dadda_stage_targets( int height )
{
  if ( height <= 2 )
  {
    return {};
  }

  std::vector<int> seq;
  seq.push_back( 2 );
  while ( seq.back() < height )
  {
    seq.push_back( ( 3 * seq.back() ) / 2 );
  }
  seq.pop_back();
  std::reverse( seq.begin(), seq.end() );
  return seq;
}

template<typename Ntk>
struct leveled_signal_t
{
  signal_t<Ntk> sig;
  int level{ 0 };
};

template<typename Ntk>
std::pair<leveled_signal_t<Ntk>, leveled_signal_t<Ntk>>
full_adder_1bit( Ntk& ntk,
                 leveled_signal_t<Ntk> const& a,
                 leveled_signal_t<Ntk> const& b,
                 leveled_signal_t<Ntk> const& cin )
{
  auto const axb = ntk.create_xor( a.sig, b.sig );
  auto const sum = ntk.create_xor( axb, cin.sig );

  auto const ab = ntk.create_and( a.sig, b.sig );
  auto const ac = ntk.create_and( a.sig, cin.sig );
  auto const bc = ntk.create_and( b.sig, cin.sig );
  auto const carry = make_or( ntk, make_or( ntk, ab, ac ), bc );

  int const level = 1 + std::max( a.level, std::max( b.level, cin.level ) );
  return { { sum, level }, { carry, level } };
}

template<typename Ntk>
struct column_reduce_result_t
{
  bool has_residual{ false };
  leveled_signal_t<Ntk> residual;
  std::vector<leveled_signal_t<Ntk>> carries;
};

template<typename Ntk>
column_reduce_result_t<Ntk> reduce_column_bits( Ntk& ntk,
                                                std::vector<leveled_signal_t<Ntk>> const& bits_in,
                                                csa_schedule_mode mode,
                                                int* fa_count,
                                                int* csa_levels )
{
  column_reduce_result_t<Ntk> out;
  if ( bits_in.empty() )
  {
    return out;
  }

  auto const zero = leveled_signal_t<Ntk>{ ntk.get_constant( false ), 0 };
  auto emit_fa = [&]( leveled_signal_t<Ntk> const& a,
                      leveled_signal_t<Ntk> const& b,
                      leveled_signal_t<Ntk> const& cin ) -> std::pair<leveled_signal_t<Ntk>, leveled_signal_t<Ntk>> {
    auto const [sum, carry] = full_adder_1bit( ntk, a, b, cin );
    if ( fa_count != nullptr )
    {
      ++( *fa_count );
    }
    if ( csa_levels != nullptr )
    {
      *csa_levels = std::max( *csa_levels, sum.level );
      *csa_levels = std::max( *csa_levels, carry.level );
    }
    return { sum, carry };
  };

  if ( mode == csa_schedule_mode::serial )
  {
    auto acc = bits_in.front();
    size_t idx = 1;

    while ( idx + 1 < bits_in.size() )
    {
      auto const [sum, carry] = emit_fa( acc, bits_in[idx], bits_in[idx + 1] );
      out.carries.push_back( carry );
      acc = sum;
      idx += 2;
    }

    if ( idx < bits_in.size() )
    {
      auto const [sum, carry] = emit_fa( acc, bits_in[idx], zero );
      out.carries.push_back( carry );
      acc = sum;
    }

    out.has_residual = true;
    out.residual = acc;
    return out;
  }

  if ( mode == csa_schedule_mode::wallace )
  {
    std::vector<leveled_signal_t<Ntk>> current = bits_in;

    while ( current.size() > 2 )
    {
      std::vector<leveled_signal_t<Ntk>> next;
      size_t i = 0;
      while ( i + 2 < current.size() )
      {
        auto const [sum, carry] = emit_fa( current[i], current[i + 1], current[i + 2] );
        out.carries.push_back( carry );
        next.push_back( sum );
        i += 3;
      }
      while ( i < current.size() )
      {
        next.push_back( current[i] );
        ++i;
      }
      current = next;
    }

    if ( current.size() == 2 )
    {
      auto const [sum, carry] = emit_fa( current[0], current[1], zero );
      out.carries.push_back( carry );
      current = { sum };
    }

    out.has_residual = true;
    out.residual = current[0];
    return out;
  }

  // Dadda mode.
  std::vector<leveled_signal_t<Ntk>> current = bits_in;
  auto const targets = dadda_stage_targets( static_cast<int>( current.size() ) );

  for ( int const target : targets )
  {
    if ( static_cast<int>( current.size() ) <= target )
    {
      continue;
    }

    std::deque<leveled_signal_t<Ntk>> work( current.begin(), current.end() );
    std::vector<leveled_signal_t<Ntk>> next;

    while ( static_cast<int>( work.size() ) > target && work.size() >= 3 )
    {
      auto const a = work.front();
      work.pop_front();
      auto const b = work.front();
      work.pop_front();
      auto const c = work.front();
      work.pop_front();
      auto const [sum, carry] = emit_fa( a, b, c );
      out.carries.push_back( carry );
      next.push_back( sum );
    }

    while ( !work.empty() )
    {
      next.push_back( work.front() );
      work.pop_front();
    }
    current = next;
  }

  while ( current.size() > 2 )
  {
    std::vector<leveled_signal_t<Ntk>> next;
    size_t i = 0;
    while ( i + 2 < current.size() )
    {
      auto const [sum, carry] = emit_fa( current[i], current[i + 1], current[i + 2] );
      out.carries.push_back( carry );
      next.push_back( sum );
      i += 3;
    }
    while ( i < current.size() )
    {
      next.push_back( current[i] );
      ++i;
    }
    current = next;
  }

  if ( current.size() == 2 )
  {
    auto const [sum, carry] = emit_fa( current[0], current[1], zero );
    out.carries.push_back( carry );
    current = { sum };
  }

  out.has_residual = true;
  out.residual = current[0];
  return out;
}

template<typename Ntk>
std::map<int, leveled_signal_t<Ntk>> csa_macro_schedule_all_columns( Ntk& ntk,
                                                                     std::vector<signal_t<Ntk>> const& raw_inputs,
                                                                     csa_schedule_mode mode,
                                                                     int* fa_count,
                                                                     int* csa_levels )
{
  std::deque<leveled_signal_t<Ntk>> raw;
  for ( auto const& x : raw_inputs )
  {
    raw.push_back( { x, 0 } );
  }

  std::deque<leveled_signal_t<Ntk>> col0_sums;
  std::map<int, std::vector<leveled_signal_t<Ntk>>> col_bits;
  col_bits[1] = {};

  auto emit_fa = [&]( leveled_signal_t<Ntk> const& a,
                      leveled_signal_t<Ntk> const& b,
                      leveled_signal_t<Ntk> const& cin ) -> std::pair<leveled_signal_t<Ntk>, leveled_signal_t<Ntk>> {
    auto const [sum, carry] = full_adder_1bit( ntk, a, b, cin );
    if ( fa_count != nullptr )
    {
      ++( *fa_count );
    }
    if ( csa_levels != nullptr )
    {
      *csa_levels = std::max( *csa_levels, sum.level );
      *csa_levels = std::max( *csa_levels, carry.level );
    }
    return { sum, carry };
  };

  while ( raw.size() >= 3 )
  {
    auto const a = raw.front();
    raw.pop_front();
    auto const b = raw.front();
    raw.pop_front();
    auto const c = raw.front();
    raw.pop_front();
    auto const [sum, carry] = emit_fa( a, b, c );
    col0_sums.push_back( sum );
    col_bits[1].push_back( carry );
  }

  std::map<int, leveled_signal_t<Ntk>> residual_by_col;
  std::vector<leveled_signal_t<Ntk>> col0_queue;
  while ( !col0_sums.empty() )
  {
    col0_queue.push_back( col0_sums.front() );
    col0_sums.pop_front();
  }
  while ( !raw.empty() )
  {
    col0_queue.push_back( raw.front() );
    raw.pop_front();
  }

  auto const reduced0 = reduce_column_bits( ntk, col0_queue, mode, fa_count, csa_levels );
  if ( reduced0.has_residual )
  {
    residual_by_col[0] = reduced0.residual;
  }
  col_bits[1].insert( col_bits[1].end(), reduced0.carries.begin(), reduced0.carries.end() );

  int current_col = 1;
  while ( true )
  {
    int next_col = -1;
    for ( auto const& [col, bits] : col_bits )
    {
      if ( col < current_col )
      {
        continue;
      }
      if ( !bits.empty() )
      {
        next_col = col;
        break;
      }
    }
    if ( next_col < 0 )
    {
      break;
    }

    auto bits = std::move( col_bits[next_col] );
    col_bits[next_col].clear();
    auto const reduced = reduce_column_bits( ntk, bits, mode, fa_count, csa_levels );
    if ( reduced.has_residual )
    {
      residual_by_col[next_col] = reduced.residual;
    }
    if ( !reduced.carries.empty() )
    {
      auto& dst = col_bits[next_col + 1];
      dst.insert( dst.end(), reduced.carries.begin(), reduced.carries.end() );
    }
    current_col = next_col + 1;
  }

  return residual_by_col;
}

template<typename Ntk>
count_summary_t<Ntk> make_leaf_count( Ntk& ntk, signal_t<Ntk> x )
{
  count_summary_t<Ntk> out;
  out.bits = { x };
  out.size = 1;
  out.depth = 0;
  (void)ntk;
  return out;
}

template<typename Ntk>
count_summary_t<Ntk> combine_counts( Ntk& ntk, count_summary_t<Ntk> const& a, count_summary_t<Ntk> const& b )
{
  assert( a.size > 0 );
  assert( b.size > 0 );

  int const out_size = a.size + b.size;
  int const w_out = width_for_size( out_size );

  auto const a_ext = zero_extend( ntk, a.bits, w_out );
  auto const b_ext = zero_extend( ntk, b.bits, w_out );

  count_summary_t<Ntk> out;
  out.bits = add_unsigned( ntk, a_ext, b_ext );
  out.size = out_size;
  out.depth = 1 + std::max( a.depth, b.depth );
  return out;
}

template<typename Ntk>
count_summary_t<Ntk> build_count_range( Ntk& ntk,
                                       std::vector<signal_t<Ntk>> const& xs,
                                       int l,
                                       int r )
{
  assert( l >= 0 );
  assert( l < r );
  assert( r <= static_cast<int>( xs.size() ) );

  if ( r - l == 1 )
  {
    return make_leaf_count( ntk, xs[static_cast<size_t>( l )] );
  }

  int const mid = l + ( r - l ) / 2;
  auto const left = build_count_range( ntk, xs, l, mid );
  auto const right = build_count_range( ntk, xs, mid, r );
  return combine_counts( ntk, left, right );
}

} // namespace detail

template<typename Ntk>
signal_t<Ntk> create_maj_popcount_exact( Ntk& ntk, std::vector<signal_t<Ntk>> const& xs, popcount_stats* st = nullptr )
{
  assert( !xs.empty() );
  assert( ( xs.size() & 1u ) == 1u );

  auto const root = detail::build_count_range( ntk, xs, 0, static_cast<int>( xs.size() ) );
  auto const threshold = make_const_vec( ntk, static_cast<int>( root.bits.size() ), static_cast<uint64_t>( xs.size() / 2u ) );

  if ( st != nullptr )
  {
    st->popcount_levels = root.depth;
    st->max_count_width = static_cast<int>( root.bits.size() );
  }

  return gt_unsigned( ntk, root.bits, threshold );
}

template<typename Ntk>
signal_t<Ntk> create_maj_popcount_strict( Ntk& ntk,
                                          std::vector<signal_t<Ntk>> const& xs,
                                          popcount_stats* st = nullptr,
                                          popcount_strict_stats* strict_st = nullptr,
                                          csa_schedule_mode schedule = csa_schedule_mode::dadda )
{
  assert( !xs.empty() );
  assert( ( xs.size() & 1u ) == 1u );

  int const n = static_cast<int>( xs.size() );
  int p = 0;
  int pow2 = 1;
  while ( pow2 < n + 1 )
  {
    pow2 <<= 1;
    ++p;
  }

  int const N = pow2 - 1;
  int const m = p;
  int const th_N = ( N + 1 ) / 2;
  int const num_fix = ( N - n ) / 2;
  assert( num_fix >= 0 );

  std::vector<signal_t<Ntk>> hw_inputs;
  hw_inputs.reserve( static_cast<size_t>( N ) );
  for ( auto const& x : xs )
  {
    hw_inputs.push_back( x );
  }
  for ( int i = 0; i < num_fix; ++i )
  {
    hw_inputs.push_back( ntk.get_constant( true ) );
  }
  for ( int i = 0; i < num_fix; ++i )
  {
    hw_inputs.push_back( ntk.get_constant( false ) );
  }
  assert( static_cast<int>( hw_inputs.size() ) == N );

  int csa_fa_count = 0;
  int csa_levels = 0;
  auto const residual_by_col = detail::csa_macro_schedule_all_columns( ntk, hw_inputs, schedule, &csa_fa_count, &csa_levels );

  std::vector<detail::leveled_signal_t<Ntk>> hw_bits;
  hw_bits.reserve( static_cast<size_t>( m ) );
  for ( int i = 0; i < m; ++i )
  {
    auto it = residual_by_col.find( i );
    if ( it != residual_by_col.end() )
    {
      hw_bits.push_back( it->second );
    }
    else
    {
      hw_bits.push_back( { ntk.get_constant( false ), 0 } );
    }
  }

  auto carry = detail::leveled_signal_t<Ntk>{ ntk.get_constant( true ), 0 };
  for ( int i = 0; i < m; ++i )
  {
    bool const bit_i = ( ( ( th_N - 1 ) >> i ) & 1 ) != 0;
    auto const t_i = detail::leveled_signal_t<Ntk>{ ntk.get_constant( bit_i ), 0 };
    auto const [sum, next_carry] = detail::full_adder_1bit( ntk, hw_bits[static_cast<size_t>( i )], t_i, carry );
    (void)sum;
    carry = next_carry;
  }

  int const comparator_fa_count = m;
  int const total_fa_count = csa_fa_count + comparator_fa_count;
  int const total_levels = csa_levels + m;

  if ( st != nullptr )
  {
    st->popcount_levels = total_levels;
    st->max_count_width = m;
  }
  if ( strict_st != nullptr )
  {
    strict_st->scaffold_p = p;
    strict_st->scaffold_inputs = N;
    strict_st->scaffold_threshold = th_N;
    strict_st->comparator_width = m;
    strict_st->num_fixed_pairs = num_fix;
    strict_st->csa_fa_count = csa_fa_count;
    strict_st->comparator_fa_count = comparator_fa_count;
    strict_st->total_fa_count = total_fa_count;
    strict_st->csa_levels = csa_levels;
    strict_st->total_levels = total_levels;
    strict_st->schedule_mode = schedule;
  }

  return carry.sig;
}

template<typename Ntk>
signal_t<Ntk> create_maj_baseline_strict( Ntk& ntk,
                                          std::vector<signal_t<Ntk>> const& xs,
                                          popcount_stats* st = nullptr,
                                          popcount_strict_stats* strict_st = nullptr,
                                          csa_schedule_mode schedule = csa_schedule_mode::dadda )
{
  return create_maj_popcount_strict( ntk, xs, st, strict_st, schedule );
}

} // namespace maj

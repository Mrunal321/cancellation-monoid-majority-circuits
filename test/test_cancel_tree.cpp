#include "maj_approx.hpp"
#include "maj_cancel_tree.hpp"
#include "maj_popcount.hpp"

#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/io/write_dot.hpp>
#include <mockturtle/networks/xag.hpp>

#include <algorithm>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace
{

struct soft_summary
{
  bool cand{ false };
  int k{ 0 };
  int size{ 0 };
  int depth{ 0 };
  int l{ 0 };
  int r{ 0 };
};

bool maj_ref( uint64_t mask, int n )
{
  int ones = 0;
  for ( int i = 0; i < n; ++i )
  {
    ones += ( ( mask >> i ) & 1ULL ) ? 1 : 0;
  }
  return ones > ( n / 2 );
}

bool maj_ref_random( std::vector<bool> const& bits )
{
  int ones = 0;
  for ( auto bit : bits )
  {
    ones += bit ? 1 : 0;
  }
  return ones > ( static_cast<int>( bits.size() ) / 2 );
}

mockturtle::xag_network build_cancel_tree_ntk( int n, maj::cancel_tree_stats* st = nullptr )
{
  mockturtle::xag_network ntk;
  std::vector<mockturtle::xag_network::signal> xs;
  xs.reserve( static_cast<size_t>( n ) );

  for ( int i = 0; i < n; ++i )
  {
    (void)i;
    xs.emplace_back( ntk.create_pi() );
  }

  auto const y = maj::create_maj_cancel_tree( ntk, xs, st );
  ntk.create_po( y );
  return mockturtle::cleanup_dangling( ntk );
}

mockturtle::xag_network build_cancel_tree_v2_ntk( int n, maj::cancel_tree_stats* st = nullptr )
{
  mockturtle::xag_network ntk;
  std::vector<mockturtle::xag_network::signal> xs;
  xs.reserve( static_cast<size_t>( n ) );

  for ( int i = 0; i < n; ++i )
  {
    (void)i;
    xs.emplace_back( ntk.create_pi() );
  }

  auto const y = maj::create_maj_cancel_tree_v2( ntk, xs, st );
  ntk.create_po( y );
  return mockturtle::cleanup_dangling( ntk );
}

mockturtle::xag_network build_popcount_ntk( int n )
{
  mockturtle::xag_network ntk;
  std::vector<mockturtle::xag_network::signal> xs;
  xs.reserve( static_cast<size_t>( n ) );

  for ( int i = 0; i < n; ++i )
  {
    (void)i;
    xs.emplace_back( ntk.create_pi() );
  }

  auto const y = maj::create_maj_popcount_exact( ntk, xs, nullptr );
  ntk.create_po( y );
  return mockturtle::cleanup_dangling( ntk );
}

mockturtle::xag_network build_popcount_strict_ntk( int n, maj::popcount_strict_stats* st = nullptr )
{
  mockturtle::xag_network ntk;
  std::vector<mockturtle::xag_network::signal> xs;
  xs.reserve( static_cast<size_t>( n ) );

  for ( int i = 0; i < n; ++i )
  {
    (void)i;
    xs.emplace_back( ntk.create_pi() );
  }

  maj::popcount_stats pop_st{};
  auto const y = maj::create_maj_popcount_strict( ntk, xs, &pop_st, st, maj::csa_schedule_mode::dadda );
  ntk.create_po( y );
  return mockturtle::cleanup_dangling( ntk );
}

mockturtle::xag_network build_approx_block3_ntk( int n, maj::approx_block3_stats* st = nullptr )
{
  mockturtle::xag_network ntk;
  std::vector<mockturtle::xag_network::signal> xs;
  xs.reserve( static_cast<size_t>( n ) );

  for ( int i = 0; i < n; ++i )
  {
    (void)i;
    xs.emplace_back( ntk.create_pi() );
  }

  auto const y = maj::create_maj_approx_block3_popcount( ntk, xs, st );
  ntk.create_po( y );
  return mockturtle::cleanup_dangling( ntk );
}

bool eval_ntk( mockturtle::xag_network const& ntk, std::vector<bool> const& assignment )
{
  mockturtle::default_simulator<bool> sim( assignment );
  auto const out = mockturtle::simulate<bool>( ntk, sim );
  return out.at( 0u );
}

bool approx_block3_ref( std::vector<bool> const& bits )
{
  int const n = static_cast<int>( bits.size() );
  int const groups = n / 3;
  int const rem = n % 3;

  std::vector<bool> compressed;
  compressed.reserve( static_cast<size_t>( groups + rem ) );

  for ( int i = 0; i < groups; ++i )
  {
    int const b = 3 * i;
    int const ones = static_cast<int>( bits[static_cast<size_t>( b )] )
                     + static_cast<int>( bits[static_cast<size_t>( b + 1 )] )
                     + static_cast<int>( bits[static_cast<size_t>( b + 2 )] );
    compressed.push_back( ones >= 2 );
  }
  for ( int i = groups * 3; i < n; ++i )
  {
    compressed.push_back( bits[static_cast<size_t>( i )] );
  }

  return maj_ref_random( compressed );
}

soft_summary soft_combine( soft_summary const& a, soft_summary const& b )
{
  assert( a.size > 0 );
  assert( b.size > 0 );
  assert( a.k >= 0 && a.k <= a.size );
  assert( b.k >= 0 && b.k <= b.size );

  soft_summary out;
  out.size = a.size + b.size;
  out.depth = 1 + std::max( a.depth, b.depth );
  out.l = a.l;
  out.r = b.r;

  if ( a.k == 0 )
  {
    out.cand = b.cand;
    out.k = b.k;
  }
  else if ( b.k == 0 )
  {
    out.cand = a.cand;
    out.k = a.k;
  }
  else if ( a.cand == b.cand )
  {
    out.cand = a.cand;
    out.k = a.k + b.k;
  }
  else if ( a.k > b.k )
  {
    out.cand = a.cand;
    out.k = a.k - b.k;
  }
  else if ( b.k > a.k )
  {
    out.cand = b.cand;
    out.k = b.k - a.k;
  }
  else
  {
    out.cand = false;
    out.k = 0;
  }

  if ( out.k == 0 )
  {
    out.cand = false;
  }

  assert( out.k >= 0 );
  assert( out.k <= out.size );

  int const width = maj::width_for_size( out.size );
  if ( width > 0 )
  {
    int const max_representable = ( 1 << width ) - 1;
    assert( out.k <= max_representable );
  }

  return out;
}

soft_summary build_soft_range( std::vector<bool> const& bits,
                               int l,
                               int r,
                               std::vector<std::vector<soft_summary>>& by_depth )
{
  if ( r - l == 1 )
  {
    soft_summary leaf;
    leaf.cand = bits.at( static_cast<size_t>( l ) );
    leaf.k = 1;
    leaf.size = 1;
    leaf.depth = 0;
    leaf.l = l;
    leaf.r = r;
    return leaf;
  }

  int const mid = l + ( r - l ) / 2;
  auto const left = build_soft_range( bits, l, mid, by_depth );
  auto const right = build_soft_range( bits, mid, r, by_depth );
  auto const out = soft_combine( left, right );

  if ( static_cast<int>( by_depth.size() ) <= out.depth )
  {
    by_depth.resize( static_cast<size_t>( out.depth + 1 ) );
  }
  by_depth[static_cast<size_t>( out.depth )].push_back( out );

  return out;
}

void dump_debug_summaries( int n, std::string const& path, int num_patterns )
{
  fs::create_directories( fs::path( path ).parent_path() );
  std::ofstream os( path, std::ios::out | std::ios::trunc );
  if ( !os )
  {
    throw std::runtime_error( "cannot open debug dump path: " + path );
  }

  uint64_t const max_patterns = ( n <= 20 ) ? ( 1ULL << n ) : static_cast<uint64_t>( num_patterns );
  uint64_t const patterns = std::min<uint64_t>( static_cast<uint64_t>( num_patterns ), max_patterns );

  for ( uint64_t mask = 0; mask < patterns; ++mask )
  {
    std::vector<bool> bits( static_cast<size_t>( n ), false );
    for ( int i = 0; i < n; ++i )
    {
      bits[static_cast<size_t>( i )] = ( ( mask >> i ) & 1ULL ) != 0ULL;
    }

    std::vector<std::vector<soft_summary>> by_depth;
    auto const root = build_soft_range( bits, 0, n, by_depth );

    os << "pattern=" << mask << " bits=";
    for ( int i = n - 1; i >= 0; --i )
    {
      os << ( bits[static_cast<size_t>( i )] ? '1' : '0' );
    }
    os << "\n";

    for ( size_t d = 1; d < by_depth.size(); ++d )
    {
      os << "  depth=" << d << "\n";
      for ( auto const& s : by_depth[d] )
      {
        os << "    range=[" << s.l << "," << s.r << ") cand=" << static_cast<int>( s.cand )
           << " k=" << s.k << " size=" << s.size << "\n";
      }
    }

    os << "  root: cand=" << static_cast<int>( root.cand ) << " k=" << root.k << " size=" << root.size << "\n\n";
  }
}

bool run_exhaustive_tests()
{
  std::vector<int> const ns = { 3, 5, 7, 9, 11 };

  for ( int n : ns )
  {
    maj::cancel_tree_stats stats_v1;
    maj::cancel_tree_stats stats_v2;
    maj::approx_block3_stats stats_approx;
    maj::popcount_strict_stats stats_pop_strict;
    auto const ntk_v1 = build_cancel_tree_ntk( n, &stats_v1 );
    auto const ntk_v2 = build_cancel_tree_v2_ntk( n, &stats_v2 );
    auto const ntk_approx = build_approx_block3_ntk( n, &stats_approx );
    auto const ntk_pop_strict = build_popcount_strict_ntk( n, &stats_pop_strict );

    uint64_t const limit = 1ULL << n;
    uint64_t approx_hits_exact = 0ULL;
    for ( uint64_t mask = 0; mask < limit; ++mask )
    {
      std::vector<bool> assignment( static_cast<size_t>( n ), false );
      for ( int i = 0; i < n; ++i )
      {
        assignment[static_cast<size_t>( i )] = ( ( mask >> i ) & 1ULL ) != 0ULL;
      }

      bool const got_v1 = eval_ntk( ntk_v1, assignment );
      bool const got_v2 = eval_ntk( ntk_v2, assignment );
      bool const got_approx = eval_ntk( ntk_approx, assignment );
      bool const got_pop_strict = eval_ntk( ntk_pop_strict, assignment );
      bool const exp = maj_ref( mask, n );
      bool const exp_approx = approx_block3_ref( assignment );
      if ( got_v1 != exp )
      {
        std::cerr << "[FAIL] exhaustive mismatch method=cancel_tree n=" << n << " mask=" << mask
                  << " got=" << got_v1 << " exp=" << exp << "\n";
        return false;
      }
      if ( got_v2 != exp )
      {
        std::cerr << "[FAIL] exhaustive mismatch method=cancel_tree_v2 n=" << n << " mask=" << mask
                  << " got=" << got_v2 << " exp=" << exp << "\n";
        return false;
      }
      if ( got_v1 != got_v2 )
      {
        std::cerr << "[FAIL] exhaustive mismatch method=v1_vs_v2 n=" << n << " mask=" << mask
                  << " v1=" << got_v1 << " v2=" << got_v2 << "\n";
        return false;
      }
      if ( got_approx != exp_approx )
      {
        std::cerr << "[FAIL] exhaustive mismatch method=approx_block3_ref n=" << n << " mask=" << mask
                  << " got=" << got_approx << " exp=" << exp_approx << "\n";
        return false;
      }
      if ( got_pop_strict != exp )
      {
        std::cerr << "[FAIL] exhaustive mismatch method=popcount_strict n=" << n << " mask=" << mask
                  << " got=" << got_pop_strict << " exp=" << exp << "\n";
        return false;
      }
      if ( got_approx == exp )
      {
        ++approx_hits_exact;
      }
    }

    std::cout << "[PASS] exhaustive method=cancel_tree n=" << n << " patterns=" << limit
              << " merge_depth=" << stats_v1.depth_merges << " max_k_width=" << stats_v1.max_k_width << "\n";
    std::cout << "[PASS] exhaustive method=cancel_tree_v2 n=" << n << " patterns=" << limit
              << " merge_depth=" << stats_v2.depth_merges << " max_k_width=" << stats_v2.max_k_width << "\n";
    std::cout << "[PASS] exhaustive method=approx_block3 n=" << n << " patterns=" << limit
              << " compressed_n=" << stats_approx.compressed_size
              << " groups3=" << stats_approx.groups3
              << " exact_match_rate=" << static_cast<double>( approx_hits_exact ) / static_cast<double>( limit ) << "\n";
    std::cout << "[PASS] exhaustive method=popcount_strict n=" << n << " patterns=" << limit
              << " scaffold_N=" << stats_pop_strict.scaffold_inputs
              << " total_fa=" << stats_pop_strict.total_fa_count << "\n";
  }

  return true;
}

bool run_random_tests()
{
  std::vector<int> const ns = { 31, 63, 127 };
  std::mt19937_64 rng( 0xC0FFEEULL );

  for ( int n : ns )
  {
    auto const cancel_ntk = build_cancel_tree_ntk( n );
    auto const cancel_v2_ntk = build_cancel_tree_v2_ntk( n );
    auto const pop_ntk = build_popcount_ntk( n );
    auto const pop_strict_ntk = build_popcount_strict_ntk( n );
    auto const approx_ntk = build_approx_block3_ntk( n );

    int approx_hits_exact = 0;

    for ( int t = 0; t < 10000; ++t )
    {
      std::vector<bool> assignment( static_cast<size_t>( n ), false );
      for ( int i = 0; i < n; ++i )
      {
        assignment[static_cast<size_t>( i )] = ( rng() & 1ULL ) != 0ULL;
      }

      bool const exp = maj_ref_random( assignment );
      bool const got_cancel = eval_ntk( cancel_ntk, assignment );
      bool const got_cancel_v2 = eval_ntk( cancel_v2_ntk, assignment );
      bool const got_pop = eval_ntk( pop_ntk, assignment );
      bool const got_pop_strict = eval_ntk( pop_strict_ntk, assignment );
      bool const got_approx = eval_ntk( approx_ntk, assignment );
      bool const exp_approx = approx_block3_ref( assignment );

      if ( got_cancel != exp )
      {
        std::cerr << "[FAIL] random mismatch (golden) n=" << n << " trial=" << t << " got=" << got_cancel << " exp=" << exp << "\n";
        return false;
      }
      if ( got_cancel != got_pop )
      {
        std::cerr << "[FAIL] random mismatch (popcount equivalence) n=" << n << " trial=" << t
                  << " cancel=" << got_cancel << " pop=" << got_pop << "\n";
        return false;
      }
      if ( got_cancel_v2 != exp )
      {
        std::cerr << "[FAIL] random mismatch (golden) method=cancel_tree_v2 n=" << n
                  << " trial=" << t << " got=" << got_cancel_v2 << " exp=" << exp << "\n";
        return false;
      }
      if ( got_cancel_v2 != got_pop )
      {
        std::cerr << "[FAIL] random mismatch (popcount equivalence) method=cancel_tree_v2 n=" << n
                  << " trial=" << t << " cancel_v2=" << got_cancel_v2 << " pop=" << got_pop << "\n";
        return false;
      }
      if ( got_cancel_v2 != got_cancel )
      {
        std::cerr << "[FAIL] random mismatch (v1_vs_v2) n=" << n
                  << " trial=" << t << " cancel=" << got_cancel << " cancel_v2=" << got_cancel_v2 << "\n";
        return false;
      }
      if ( got_pop_strict != got_pop )
      {
        std::cerr << "[FAIL] random mismatch (strict_vs_popcount) n=" << n
                  << " trial=" << t << " pop_strict=" << got_pop_strict << " pop=" << got_pop << "\n";
        return false;
      }
      if ( got_approx != exp_approx )
      {
        std::cerr << "[FAIL] random mismatch (approx software model) n=" << n
                  << " trial=" << t << " approx_ntk=" << got_approx << " approx_sw=" << exp_approx << "\n";
        return false;
      }
      if ( got_approx == exp )
      {
        ++approx_hits_exact;
      }
    }

    std::cout << "[PASS] random method=cancel_tree n=" << n << " vectors=10000\n";
    std::cout << "[PASS] random method=cancel_tree_v2 n=" << n << " vectors=10000\n";
    std::cout << "[PASS] random method=popcount_strict n=" << n << " vectors=10000\n";
    std::cout << "[PASS] random method=approx_block3 n=" << n << " vectors=10000"
              << " exact_match_rate=" << ( static_cast<double>( approx_hits_exact ) / 10000.0 ) << "\n";
  }

  return true;
}

bool emit_debug_artifacts()
{
  fs::create_directories( "artifacts" );

  dump_debug_summaries( 7, "artifacts/debug_summaries_n7.txt", 16 );
  dump_debug_summaries( 11, "artifacts/debug_summaries_n11.txt", 16 );

  auto const ntk7 = build_cancel_tree_ntk( 7 );
  mockturtle::write_dot( ntk7, "artifacts/maj_7_cancel_tree.dot" );

  std::cout << "[PASS] debug artifacts written under artifacts/\n";
  return true;
}

} // namespace

int main()
{
  try
  {
    if ( !run_exhaustive_tests() )
    {
      return 1;
    }
    if ( !run_random_tests() )
    {
      return 1;
    }
    if ( !emit_debug_artifacts() )
    {
      return 1;
    }

    std::cout << "[PASS] all cancel_tree checks succeeded\n";
    return 0;
  }
  catch ( std::exception const& e )
  {
    std::cerr << "[FAIL] exception: " << e.what() << "\n";
    return 1;
  }
}

#include "maj_approx.hpp"
#include "maj_cancel_tree.hpp"
#include "maj_popcount.hpp"

#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/io/write_dot.hpp>
#include <mockturtle/io/write_verilog.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/views/depth_view.hpp>

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;

namespace
{

struct cli_options
{
  std::string maj_impl{ "cancel_tree" };
  std::string logic_repr{ "xag" };
  std::string strict_schedule{ "dadda" };
  int n{ 0 };
  std::string out_path;
  std::string format;
  std::string module_name;
  std::string stats_json;
  std::string dot_path;
  bool show_help{ false };
  bool show_version{ false };
};

void print_usage()
{
  std::cout
      << "Usage: majgen --maj_impl <cancel_tree|cancel_tree_v2|boyermoore_tree|boyermoore_tree_v2|approx_block3|approx_block3_popcount|popcount|popcount_strict|baseline_strict> --n <odd> --out <file> [options]\n"
      << "Options:\n"
      << "  --logic_repr <xag|aig|mig>         Internal logic representation for output/stats (default: xag)\n"
      << "  --strict_schedule <serial|wallace|dadda>\n"
      << "                                     CSA schedule for popcount_strict/baseline_strict (default: dadda)\n"
      << "  --format <verilog|blif|aig|dot>   Output format (default inferred from extension, else verilog)\n"
      << "  --module <name>                    Verilog module name\n"
      << "  --stats_json <file>                Write generator/network stats JSON\n"
      << "  --dot <file>                       Also write DOT graph for the selected representation\n"
      << "  --version                           Print tool version\n"
      << "  --help                              Show this message\n";
}

std::string infer_format( std::string const& path )
{
  auto const ext = fs::path( path ).extension().string();
  if ( ext == ".v" )
  {
    return "verilog";
  }
  if ( ext == ".blif" )
  {
    return "blif";
  }
  if ( ext == ".aig" || ext == ".aiger" )
  {
    return "aig";
  }
  if ( ext == ".dot" )
  {
    return "dot";
  }
  return "";
}

cli_options parse_args( int argc, char** argv )
{
  cli_options opt;

  auto require_value = [&]( int& i ) -> std::string {
    if ( i + 1 >= argc )
    {
      throw std::runtime_error( std::string( "missing value for " ) + argv[i] );
    }
    ++i;
    return argv[i];
  };

  for ( int i = 1; i < argc; ++i )
  {
    std::string const arg = argv[i];
    if ( arg == "--help" || arg == "-h" )
    {
      opt.show_help = true;
      continue;
    }
    if ( arg == "--version" )
    {
      opt.show_version = true;
      continue;
    }
    if ( arg == "--maj_impl" )
    {
      opt.maj_impl = require_value( i );
      continue;
    }
    if ( arg == "--n" )
    {
      opt.n = std::stoi( require_value( i ) );
      continue;
    }
    if ( arg == "--strict_schedule" )
    {
      opt.strict_schedule = require_value( i );
      continue;
    }
    if ( arg == "--logic_repr" )
    {
      opt.logic_repr = require_value( i );
      continue;
    }
    if ( arg == "--out" )
    {
      opt.out_path = require_value( i );
      continue;
    }
    if ( arg == "--format" )
    {
      opt.format = require_value( i );
      continue;
    }
    if ( arg == "--module" )
    {
      opt.module_name = require_value( i );
      continue;
    }
    if ( arg == "--stats_json" )
    {
      opt.stats_json = require_value( i );
      continue;
    }
    if ( arg == "--dot" )
    {
      opt.dot_path = require_value( i );
      continue;
    }

    throw std::runtime_error( "unknown argument: " + arg );
  }

  return opt;
}

bool is_valid_method( std::string const& method )
{
  return method == "cancel_tree" || method == "cancel_tree_v2"
         || method == "boyermoore_tree" || method == "boyermoore_tree_v2"
         || method == "approx_block3" || method == "approx_block3_popcount"
         || method == "popcount"
         || method == "popcount_strict" || method == "baseline_strict";
}

bool is_valid_logic_repr( std::string const& logic_repr )
{
  return logic_repr == "xag" || logic_repr == "aig" || logic_repr == "mig";
}

void ensure_parent_dir( std::string const& path )
{
  fs::path const p{ path };
  if ( p.has_parent_path() )
  {
    fs::create_directories( p.parent_path() );
  }
}

void write_stats_json( std::string const& path,
                       std::unordered_map<std::string, std::string> const& fields,
                       std::unordered_map<std::string, int> const& int_fields )
{
  ensure_parent_dir( path );
  std::ofstream os( path, std::ios::out | std::ios::trunc );
  if ( !os )
  {
    throw std::runtime_error( "failed to open stats JSON path: " + path );
  }

  os << "{\n";
  bool first = true;

  auto emit_sep = [&]() {
    if ( !first )
    {
      os << ",\n";
    }
    first = false;
  };

  for ( auto const& [k, v] : fields )
  {
    emit_sep();
    os << "  \"" << k << "\": \"" << v << "\"";
  }
  for ( auto const& [k, v] : int_fields )
  {
    emit_sep();
    os << "  \"" << k << "\": " << v;
  }
  os << "\n}\n";
}

template<typename Ntk>
int count_and_gates( Ntk const& ntk )
{
  int and_count = 0;
  ntk.foreach_gate( [&]( auto const& n ) {
    if ( ntk.is_and( n ) )
    {
      ++and_count;
    }
  } );
  return and_count;
}

template<typename Ntk>
int count_maj_gates( Ntk const& ntk )
{
  int maj_count = 0;
  ntk.foreach_gate( [&]( auto const& n ) {
    if ( ntk.is_maj( n ) )
    {
      ++maj_count;
    }
  } );
  return maj_count;
}

} // namespace

int main( int argc, char** argv )
{
  try
  {
    auto const opt = parse_args( argc, argv );

    if ( opt.show_help )
    {
      print_usage();
      return 0;
    }

    if ( opt.show_version )
    {
      std::cout << "majgen 0.1.0\n";
      return 0;
    }

    if ( !is_valid_method( opt.maj_impl ) )
    {
      throw std::runtime_error(
          "--maj_impl must be one of cancel_tree, cancel_tree_v2, boyermoore_tree, boyermoore_tree_v2, approx_block3, approx_block3_popcount, popcount, popcount_strict, baseline_strict" );
    }
    if ( !is_valid_logic_repr( opt.logic_repr ) )
    {
      throw std::runtime_error( "--logic_repr must be one of xag, aig, mig" );
    }

    auto const strict_schedule = maj::parse_csa_schedule_mode( opt.strict_schedule );
    if ( opt.n <= 0 || ( opt.n & 1 ) == 0 )
    {
      throw std::runtime_error( "--n must be a positive odd integer" );
    }
    if ( opt.out_path.empty() )
    {
      throw std::runtime_error( "--out is required" );
    }

    std::string format = opt.format;
    if ( format.empty() )
    {
      format = infer_format( opt.out_path );
      if ( format.empty() )
      {
        format = "verilog";
      }
    }

    if ( format != "verilog" && format != "blif" && format != "aig" && format != "dot" )
    {
      throw std::runtime_error( "unsupported --format: " + format );
    }

    std::string module_name = opt.module_name;
    if ( module_name.empty() )
    {
      module_name = "majority_" + opt.maj_impl + "_n" + std::to_string( opt.n );
    }

    mockturtle::xag_network ntk;
    std::vector<mockturtle::xag_network::signal> xs;
    xs.reserve( static_cast<size_t>( opt.n ) );

    for ( int i = 0; i < opt.n; ++i )
    {
      (void)i;
      xs.emplace_back( ntk.create_pi() );
    }

    maj::cancel_tree_stats cancel_stats{};
    maj::popcount_stats pop_stats{};
    maj::popcount_strict_stats strict_stats{};
    maj::approx_block3_stats approx_stats{};

    mockturtle::xag_network::signal y = ntk.get_constant( false );
    if ( opt.maj_impl == "cancel_tree" || opt.maj_impl == "boyermoore_tree" )
    {
      y = maj::create_maj_cancel_tree( ntk, xs, &cancel_stats );
    }
    else if ( opt.maj_impl == "cancel_tree_v2" || opt.maj_impl == "boyermoore_tree_v2" )
    {
      y = maj::create_maj_cancel_tree_v2( ntk, xs, &cancel_stats );
    }
    else if ( opt.maj_impl == "approx_block3" || opt.maj_impl == "approx_block3_popcount" )
    {
      y = maj::create_maj_approx_block3_popcount( ntk, xs, &approx_stats );
      pop_stats.popcount_levels = approx_stats.popcount_levels;
    }
    else if ( opt.maj_impl == "popcount_strict" || opt.maj_impl == "baseline_strict" )
    {
      y = maj::create_maj_popcount_strict( ntk, xs, &pop_stats, &strict_stats, strict_schedule );
    }
    else
    {
      y = maj::create_maj_popcount_exact( ntk, xs, &pop_stats );
    }

    ntk.create_po( y );

    auto const clean_xag = mockturtle::cleanup_dangling( ntk );
    auto const [and_count, xor_count] = maj::count_and_xor_gates( clean_xag );
    mockturtle::depth_view depth_xag{ clean_xag };
    auto const inv_count_xag = maj::count_complemented_edges( clean_xag );
    auto const logic_depth_xag = static_cast<int>( depth_xag.depth() );

    auto const clean_aig = mockturtle::cleanup_dangling<mockturtle::xag_network, mockturtle::aig_network>( clean_xag );
    mockturtle::depth_view depth_aig{ clean_aig };
    auto const aig_and_count = count_and_gates( clean_aig );
    auto const inv_count_aig = maj::count_complemented_edges( clean_aig );
    auto const logic_depth_aig = static_cast<int>( depth_aig.depth() );

    auto const clean_mig = mockturtle::cleanup_dangling<mockturtle::xag_network, mockturtle::mig_network>( clean_xag );
    mockturtle::depth_view depth_mig{ clean_mig };
    auto const mig_maj_count = count_maj_gates( clean_mig );
    auto const inv_count_mig = maj::count_complemented_edges( clean_mig );
    auto const logic_depth_mig = static_cast<int>( depth_mig.depth() );

    int repr_node_count = static_cast<int>( clean_xag.num_gates() );
    int repr_inv_count = inv_count_xag;
    int repr_depth = logic_depth_xag;

    if ( opt.logic_repr == "aig" )
    {
      repr_node_count = static_cast<int>( clean_aig.num_gates() );
      repr_inv_count = inv_count_aig;
      repr_depth = logic_depth_aig;
    }
    else if ( opt.logic_repr == "mig" )
    {
      repr_node_count = static_cast<int>( clean_mig.num_gates() );
      repr_inv_count = inv_count_mig;
      repr_depth = logic_depth_mig;
    }

    ensure_parent_dir( opt.out_path );

    if ( format == "verilog" )
    {
      mockturtle::write_verilog_params ps;
      ps.module_name = module_name;
      if ( opt.logic_repr == "xag" )
      {
        mockturtle::write_verilog( clean_xag, opt.out_path, ps );
      }
      else if ( opt.logic_repr == "aig" )
      {
        mockturtle::write_verilog( clean_aig, opt.out_path, ps );
      }
      else
      {
        mockturtle::write_verilog( clean_mig, opt.out_path, ps );
      }
    }
    else if ( format == "blif" )
    {
      if ( opt.logic_repr == "xag" )
      {
        mockturtle::write_blif( clean_xag, opt.out_path );
      }
      else if ( opt.logic_repr == "aig" )
      {
        mockturtle::write_blif( clean_aig, opt.out_path );
      }
      else
      {
        mockturtle::write_blif( clean_mig, opt.out_path );
      }
    }
    else if ( format == "aig" )
    {
      if ( opt.logic_repr == "mig" )
      {
        auto const mig_as_aig = mockturtle::cleanup_dangling<mockturtle::mig_network, mockturtle::aig_network>( clean_mig );
        mockturtle::write_aiger( mig_as_aig, opt.out_path );
      }
      else
      {
        mockturtle::write_aiger( clean_aig, opt.out_path );
      }
    }
    else
    {
      if ( opt.logic_repr == "xag" )
      {
        mockturtle::write_dot( clean_xag, opt.out_path );
      }
      else if ( opt.logic_repr == "aig" )
      {
        mockturtle::write_dot( clean_aig, opt.out_path );
      }
      else
      {
        mockturtle::write_dot( clean_mig, opt.out_path );
      }
    }

    if ( !opt.dot_path.empty() )
    {
      ensure_parent_dir( opt.dot_path );
      if ( opt.logic_repr == "xag" )
      {
        mockturtle::write_dot( clean_xag, opt.dot_path );
      }
      else if ( opt.logic_repr == "aig" )
      {
        mockturtle::write_dot( clean_aig, opt.dot_path );
      }
      else
      {
        mockturtle::write_dot( clean_mig, opt.dot_path );
      }
    }

    bool const is_strict = ( opt.maj_impl == "popcount_strict" || opt.maj_impl == "baseline_strict" );
    std::string const strict_mode_str = is_strict ? maj::to_string( strict_stats.schedule_mode ) : "";

    if ( !opt.stats_json.empty() )
    {
      write_stats_json(
          opt.stats_json,
          {
              { "method", opt.maj_impl },
              { "logic_repr", opt.logic_repr },
              { "strict_schedule_mode", strict_mode_str },
              { "module", module_name },
              { "file", fs::absolute( opt.out_path ).string() },
          },
          {
              { "n", opt.n },
              { "node_count", repr_node_count },
              { "inv_count", repr_inv_count },
              { "xag_and_count", and_count },
              { "xag_xor_count", xor_count },
              { "aig_and_count", aig_and_count },
              { "mig_maj_count", mig_maj_count },
              { "ntk_depth_logic", repr_depth },
              { "strict_scaffold_p", strict_stats.scaffold_p },
              { "strict_scaffold_inputs", strict_stats.scaffold_inputs },
              { "strict_scaffold_threshold", strict_stats.scaffold_threshold },
              { "strict_comparator_width", strict_stats.comparator_width },
              { "strict_num_fixed_pairs", strict_stats.num_fixed_pairs },
              { "strict_csa_fa_count", strict_stats.csa_fa_count },
              { "strict_comparator_fa_count", strict_stats.comparator_fa_count },
              { "strict_total_fa_count", strict_stats.total_fa_count },
              { "strict_csa_levels", strict_stats.csa_levels },
              { "strict_total_levels", strict_stats.total_levels },
              { "xag_node_count", static_cast<int>( clean_xag.num_gates() ) },
              { "aig_node_count", static_cast<int>( clean_aig.num_gates() ) },
              { "mig_node_count", static_cast<int>( clean_mig.num_gates() ) },
              { "xag_inv_count", inv_count_xag },
              { "aig_inv_count", inv_count_aig },
              { "mig_inv_count", inv_count_mig },
              { "xag_depth_logic", logic_depth_xag },
              { "aig_depth_logic", logic_depth_aig },
              { "mig_depth_logic", logic_depth_mig },
              { "cancel_merge_depth", cancel_stats.depth_merges },
              { "cancel_max_k_width", cancel_stats.max_k_width },
              { "cancel_num_merges", cancel_stats.num_merges },
              { "popcount_levels", pop_stats.popcount_levels },
              { "approx_group_size", approx_stats.group_size },
              { "approx_compressed_n", approx_stats.compressed_size },
              { "approx_groups3", approx_stats.groups3 },
              { "approx_passthrough_bits", approx_stats.passthrough_bits },
          } );
    }

    return 0;
  }
  catch ( std::exception const& e )
  {
    std::cerr << "[majgen] error: " << e.what() << "\n";
    return 1;
  }
}

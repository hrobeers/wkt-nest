#include <cstdlib>
#include <string>
#include <iostream>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "../version_autogen.hpp"

#include <boost/geometry/io/svg/svg_mapper.hpp>
#include "wkt-nest/wktio.hpp"
#include "wkt-nest/bbpack.hpp"

using namespace wktnest;

inline void showHelp(const po::options_description &cmdline_options, std::ostream &stream = std::cout)
{
    stream << cmdline_options << "\n";
}

int main(int argc, char *argv[])
{
  try
  {
    //
    // WKT-nest options
    //
    po::options_description nesting_params("Nesting options");
    nesting_params.add_options()
      ("sort,s", po::value<std::string>(),
       "Sorting strategy before packing. Defaults to 'none'.\n"
       "['none', 'height']");


    //
    // Positional options
    //
    po::positional_options_description positional_options;
    positional_options.add("input-file", -1);

    po::options_description positional_params;
    positional_params.add_options()
      ("input-file", po::value<std::string>(), "Optional input file (can also be provided through stdin)");


    //
    // Generic options
    //
    po::options_description generic_params("Generic options");
    generic_params.add_options()
      ("version,v", "Print version string")
      ("help,h", "Print help");


    //
    // Process the actual command line arguments given by the user
    //
    po::options_description cmdline_options;
    cmdline_options.add(positional_params).add(nesting_params).add(generic_params);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(positional_options).run(), vm);
    po::notify(vm);

    // Process the generic options
    bool do_exit = false;
    if (vm.count("version")) {
      std::cout << MAJOR_VERSION << '.'
                << MINOR_VERSION << '.'
                << REVISION << '.'
                << BUILD_NUMBER << '\n'
                << COMMIT_HASH << '\n';
      do_exit = true;
    }
    if (vm.count("help")) {
      showHelp(cmdline_options);
      do_exit = true;
    }
    if (do_exit) exit(EXIT_SUCCESS);


    // Read the cmd arguments to options struct
    nesting_opts opts;
    opts.sorting = vm.count("sort") &&
      vm["sort"].as<std::string>()=="height"? SORTING::HEIGHT : SORTING::NONE;

    //
    // Run the application
    //

    if (true)
    {
      std::optional<box_t> b = read_box(std::cin);
      if (!b) {
        std::cerr << "No bin box provided" << std::endl;
        exit(EXIT_FAILURE);
      }
      std::vector<polygon_t> ps = read_polygons(std::cin);

      auto state = bbpack::init(*b, opts);
      std::vector<matrix_t> fit = bbpack::fit(state, ps);

      // Declare a stream and an SVG mapper
      boost::geometry::svg_mapper<point_t> mapper(std::cout, 400, 400);

      // Add geometries such that all these geometries fit on the map
      mapper.add(*b);
      for (auto item : state.items)
        mapper.add(*item.bbox());

      mapper.map(*b, "fill-opacity:0.5;fill:rgb(153,204,0);stroke:rgb(153,204,0);stroke-width:2", 5);
      //for (auto p : ps)
      for (auto item : state.items) {
        // Draw the geometries on the SVG map, using a specific SVG style
        mapper.map(*item.bbox(), "fill-opacity:0.3;fill:rgb(51,51,153);stroke:rgb(51,51,153);stroke-width:2");
        mapper.map(*item.polygon(), "fill-opacity:0.3;fill:rgb(212,0,0);stroke:rgb(212,0,0);stroke-width:2");
      }

      // Destructor of map will be called - adding </svg>
      // Destructor of stream will be called, closing the file
    }
    else
    {
      showHelp(cmdline_options, std::cerr);
    }

    exit(EXIT_SUCCESS);
  }
  catch (std::exception &ex)
  {
    std::cerr << ex.what() << std::endl;
    exit(EXIT_FAILURE);
  }
}

#include <vector>
#include <map>
#include <list>
#include <iostream>
#include <boost/program_options.hpp>


namespace po = boost::program_options;


po::variables_map parse_program_options(const int argc, const char *const argv[]) {
    po::variables_map args;

    po::options_description desc("Allowed options");
    desc.add_options()
            ("help,help", "Show brief usage message")
            ("atol,a", po::value<double>()->default_value(1e-8), "atol")
            ("rtol,r", po::value<double>()->default_value(1e-8), "rtol")
            ("tol,t", po::value<double>()->default_value(1e-8), "tol")
            ("iter_num,z", po::value<double>()->default_value(3.14f), "iteration number")
            ("seed,h", po::value<int>()->default_value(-1), "iteration number")
            ("start_time,i", po::value<double>()->default_value(0.0f), "start_time")
            ("f,f", "print_flux")
            ("j,j", "enable_species_stats")
            ("k,k", "remove_zero")
            ("n,n", "print_save_net")
            ("s,s", "save file")
            ("s,s", "continuation")
            ("e,e", "print_end_net")
            ("c,c", "check_steady_state")
            ("u,u", "gillespie_update_interval")
            ("v,v", "verbose")
            ("?,?", "verbose")
            ("b,b", "use solver")
            ("d,d", "use dense solver")
            ("output_prefix,o", po::value<std::string>()->default_value(""), "outprefix")
            ("group_file,g", po::value<std::string>()->default_value(""), "group_file")
            ("max_steps,M", po::value<std::string>()->default_value(""), "maxSteps");
    ("step_int,I", po::value<std::string>()->default_value(""), "stepInterval");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    return vm;
}

char *usage = (char *) "run_network  [-bcdefkmsvx] [-a atol] [-g groupfile] [-h seed] [-i start_time] [-o outprefix] [-r rtol] [-t tol] [-z iteration number]";

int main2(int argc, char *argv[]) {

    po::variables_map options = parse_program_options(argc, argv);
    const auto atol = options["atol"].as<double>();
    const auto rtol = options["rtol"].as<double>();
    const auto seed = options["h"].as<int>();
    const auto t_start = options["start_time"].as<double>();
    const auto outtime = options["iter_num"].as<double>();
    const std::string max_steps = options["max_steps"].as<std::string>();
    const std::string step_int = options["step_int"].as<std::string>();
    const std::string outpre = options["outprefix"].as<std::string>();

    double maxSteps = INFINITY;

    if (max_steps == "INFINITY") { maxSteps = INFINITY; }
    else if (max_steps == "LONG_MAX") { maxSteps = LONG_MAX; }
    else if (max_steps == "INT_MAX") { maxSteps = INT_MAX; }


    double stepInterval = INFINITY;

    if (step_int == "INFINITY") { stepInterval = INFINITY; }
    else if (step_int == "LONG_MAX") { stepInterval = LONG_MAX; }
    else if (step_int == "INT_MAX") { stepInterval = INT_MAX; }
    else stepInterval = std::stod(step_int);

    enum {
        SSA, CVODE, EULER, RKCS, PLA, HAS
    };
    enum {
        DENSE, GMRES, DENSE_J, GMRES_J
    };

    int propagator = CVODE;
    int SOLVER = DENSE;

    const int check_steady_state = (options.count("c")) ? 1 : 0;
    const int print_end_net = (options.count("e")) ? 1 : 0;
    const int print_flux = (options.count("f")) ? 1 : 0;
    const int enable_species_stats = (options.count("j")) ? 1 : 0;
    const int remove_zero = (options.count("k")) ? 1 : 0;
    const int print_save_net = (options.count("n")) ? 1 : 0;
    const int save_file = (options.count("s")) ? 1 : 0;
    const int continuation = (options.count("x")) ? 1 : 0;
    const int verbose = (options.count("v")) ? 1 : 0;
    const int gillespie_update_interval = (options.count("x")) ? 1 : 0;
    propagator = (options.count("m")) ? SSA : CVODE;


    std::string group_input_file_name = nullptr;
    if (options.count("g")) { group_input_file_name = options["g"].as<std::string>(); }

    // set solver options
    if (options.count("d")) { SOLVER = GMRES; }
    if (options.count("b")) {
        if (SOLVER == DENSE) { SOLVER = DENSE_J; }
        else SOLVER = GMRES_J;
    }

    std::cout << "rtol " << rtol << std::endl;
    std::cout << "atol " << atol << std::endl;
    std::cout << "maxSteps " << maxSteps << std::endl;

    fflush(stdout);
    // Variables //
    int i;
    char *netfile_name, *network_name;

    char *save_file_name;
    FILE *netfile, *conc_file, *group_file/*, *func_file*/, *out, *flux_file, *species_stats_file;
    int net_line_number, group_line_number, n_read;

    int n, n_sample;
    double t, dt;
    double sample_time, *sample_times = nullptr/*, *st, t1*/;
    char c, buf[1000];
    int argleft, iarg = 1;


    double *conc, *conc_last, *derivs;
    // Allowed propagator types

    double scalelevel = 0.0;
    bool pScaleChecker = false;


    std::string pla_config; // No default
    bool print_cdat = true;
    bool print_func = false;
    bool additional_pla_output = false; // Print PLA-specific data (e.g., rxn classifications)
    bool print_on_stop = true; // Print to file if stopping condition met?

    std::string stop_string = "0";

    /* Process input options */
    while (argv[iarg][0] == '-') {
        c = argv[iarg++][1];
        switch (c) {
            case 'p':
                if (strcmp(argv[iarg], "ssa") == 0 || strcmp(argv[iarg], "has") == 0) propagator = SSA;
                else if (strcmp(argv[iarg], "cvode") == 0) propagator = CVODE;
                else if (strcmp(argv[iarg], "euler") == 0) propagator = EULER;
                else if (strcmp(argv[iarg], "rkcs") == 0) propagator = RKCS;
                else if (strcmp(argv[iarg], "pla") == 0) {
                    propagator = PLA;
                    if (argv[iarg + 1][0] != '-') pla_config = argv[++iarg];
                    else {
                        std::cout << "ERROR: To use the pla you must specify a simulation configuration."
                                     " Please try again." << std::endl;
                        exit(1);
                    }
                } else {
                    fprintf(stderr, "ERROR: Unrecognized propagator type %s.\n", argv[iarg]);
                    exit(1);
                }
                iarg++;
                break;

            case '-': // Process long options
//			cout << argv[iarg-1] << " ";
                std::string long_opt(argv[iarg - 1]);
                long_opt = long_opt.substr(2); // remove '--'
                //
                // Print to .cdat
                if (long_opt == "cdat") {
                    if (atoi(argv[iarg]) <= 0) {
                        print_cdat = false;
                        std::cout << "Suppressing concentrations (.cdat) output" << std::endl;
                    }
                }
                    // Print to .fdat
                else if (long_opt == "fdat") {
                    if (atoi(argv[iarg]) > 0) {
                        print_func = true;
                        std::cout << "Activating functions output (to .gdat)" << std::endl;
                    }
                }
                    // Print additional PLA data (e.g., rxn classifications)
                else if (long_opt == "pla_output") {
                    if (atoi(argv[iarg]) > 0) {
                        additional_pla_output = true;
                    }
                } else if (long_opt == "scalelevel") {
                    scalelevel = rint(atof(argv[iarg]));
                    if (scalelevel <= 1.0) {
                        std::cout << "Scaling target is too small (<= 1), using SSA without any scaling" << std::endl;
                        propagator = SSA;
                    } else {
                        propagator = HAS;
                        std::cout << "Using scaling method to accelerate simulation" << std::endl;
                    }
                } else if (long_opt == "check_product_scale") {
                    if (atoi(argv[iarg]) != 0) {
                        std::cout
                                << "The heterogeneous adaptive scaling method is also checking the scale of products (right hand side)"
                                << std::endl;
                        pScaleChecker = true;
                    }
                } else if (long_opt == "stop_cond") {
                    stop_string = (std::string) argv[iarg++];
                    std::cout << "Stopping condition specified: " << stop_string;
                    if (atoi(argv[iarg]) <= 0) {
                        print_on_stop = false;
                        std::cout << " (print-on-stop disabled)";
                    }
                    std::cout << std::endl;
//				cout << stop_string << endl;
                }
                    //...
                else {
//				cout << endl;
                    std::cout << "Sorry, don't recognize your long option "
                              << argv[iarg - 1] << ". Please try again." << std::endl;
                }
                iarg++;
                //
                break;
        }
    }

    /* Check input options for consistency */

    /* Check for correct number of input args */
    argleft = argc - iarg;

    /* Get net file name */
    netfile_name = strdup(argv[iarg++]);

    /* Process sample times */
    if ((argleft = argc - iarg) == 2) {
        /* input is sample_time n_sample */
        sample_time = atof(argv[iarg++]);
        n_sample = (int) atof(argv[iarg++]); // Read as float and cast to int to allow for exponential format
    } else {
        /* input is t1 t2 ... tn */
        n_sample = argleft;
        std::vector<double> st;
        std::vector<bool> keep;
        st.push_back(t_start);
        keep.push_back(true);

        // Collect all sample times
        for (int j = 0; j < n_sample; j++) {
            st.push_back(atof(argv[iarg++]));
            keep.push_back(true);
        }
        if (t_start > st[st.size() - 1]) { // BNG appends t_end to the sample_times array
            std::cout << "WARNING: t_start > t_end. Setting t_end = t_start, simulation will not run." << std::endl;
            st[st.size() - 1] = t_start;
        }
        double t_end = st[st.size() - 1];

        // Flag sample times <= t_start and >= t_end for removal
        for (unsigned int j = 1; j < st.size() - 1; j++) {
            if (st[j] <= t_start || st[j] >= t_end) {
//				cout << ": ERASE";
                keep.at(j) = false;
                n_sample--;
            }
//			cout << endl;
        }

        // Fill up sample_times array

        std::vector<double> sample_times;
        int k = 0;
        for (unsigned int j = 0; j < st.size(); j++) {
            if (keep.at(j)) {
                sample_times.push_back(st[j]);
                k++;
            }
        }
        // Error check
        if (k != n_sample + 1) {
            std::cout << "Oops, something went wrong while processing sample_times." << std::endl;
            exit(1);
        }

        // Make sure there are at least 2 elements (t_start and t_end)
        if (n_sample < 1) {
            fprintf(stderr, "ERROR: There must be at least one sample time (t_end).\n");
            exit(1);
        }

        // Check that final array is in ascending order with no negative elements
        for (i = 0; i <= n_sample; ++i) {
            if (sample_times[i] < 0.0) {
                fprintf(stderr, "ERROR: Negative sample times are not allowed.\n");
                exit(1);
            }
            if (i == 0) continue;
//			if (sample_times[i] <= sample_times[i-1]) {
            if (sample_times[i] < sample_times[i - 1]) { // Handle case where n_sample=2 and t_start=t_end
                fprintf(stderr, "ERROR: Sample times must be in ascending order.\n");
                exit(1);
            }
        }
    }
}

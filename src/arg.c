/*==============================================================================
arg.c : parse command line arguments
(C) 2008-2015 Jens Kleinjung and Alessandro Pandini
Read the COPYING file for license information.
==============================================================================*/

#include "arg.h"
#include "config.h"

/*____________________________________________________________________________*/
/** print version */
static void print_version()
{
	fprintf(stdout, "\n%s %s\n", PACKAGE, VERSION);
}

/*____________________________________________________________________________*/
/** print header */
static void print_header()
{
    fprintf(stdout, "\nOPTICS: single-linkage point ordering\n");
}

/*____________________________________________________________________________*/
/** print license */
static void print_license()
{
    fprintf(stderr, "\nCopyright (C) 2008-2015 by Jens Kleinjung and Alessandro Pandini\n"
					"You are welcome to redistribute it under certain conditions.\n"
                    "Read the COPYING file for license information.\n\n");
}

/*____________________________________________________________________________*/
/** print citation */
static void print_citation()
{
	printf("\n- Pandini et al., BMC Bioinformatics 11:97, 2010.\n");
	printf("- Ankerst et al., Proc.ACM SIGMOD'99 Int. Conf. on Management of Data, Philadelphia PA, 1999.\n");
	printf("- Daszykowski et al., J. Chem. Inf. Comput. Sci. 42:500-507, 2002.\n\n");
}

/*____________________________________________________________________________*/
/** set defaults */
static void set_defaults(Arg *arg)
{
    arg->dataInFileName = "input.dat";
    arg->dataOutFileName = "output.dat";
    arg->clusterOutFileName = "cluster.dat";
    arg->centerOutFileName = "center.dat";
    arg->uniqueOutFileName = "unique.dat";
    arg->eps = FLT_MAX;
    arg->minPts = 1;
    arg->w = 1;
    arg->outPathName = ".";
}

/*____________________________________________________________________________*/
/** parse command line long_options */
int parse_args(int argc, char **argv, Arg *arg)
{
    extern int silent;

    int c;
    const char usage[] = "\noptics [--datafile ...] [OPTIONS ...]\n\
      --datafile <filename>                (mode: mandatory, type: char,  default: input.dat)\n\
      --outpath <path>                     (mode: optional,  type: char,  default: .)\n\
      --outputfile <filename>              (mode: optional,  type: char,  default: output.dat)\n\
      --clusterfile <filename>             (mode: optional,  type: char,  default: cluster.dat)\n\
      --centerfile <filename>              (mode: optional,  type: char,  default: center.dat)\n\
      --uniquefile <filename>              (mode: optional,  type: char,  default: unique.dat)\n\
      --eps <epsilon cutoff>               (mode: optional,  type: float, default: FLT_MAX)\n\
      --minpts <min number of neighbours>  (mode: optional,  type: int,   default: 1)\n\
      --w <window size for string version> (mode: optional,  type: int,   default: 1)\n\
      --silent                             (mode: optional,  no argument, default: off)\n\
      --cite                               (mode: optional , type: no_arg, default: off)\n\
      --version                            (mode: optional , type: no_arg, default: off)\n\
      --help\n\n\
    optics_ang : input data in angle coordinates\n\
    optics_str : input data in string format\n\
    optics_vec : input data in vector (of arbitrary length) coordinates\n\
    optics_xyz : input data in Euclidean xyz coordinates\n\
    Use the configure options described in 'README Install/Uninstall' to compile a specific verion.\n\n";

    set_defaults(arg);

    if (argc < 2) {
		print_header();
        fprintf(stderr, "%s", usage);
		print_license();
        exit(1);
    }

    /** long option definition */
    static struct option long_options[] =
    {
        {"datafile", required_argument, 0, 1},
        {"outpath", required_argument, 0, 10},
        {"outputfile", required_argument, 0, 2},
        {"clusterfile", required_argument, 0, 3},
        {"centerfile", required_argument, 0, 4},
        {"uniquefile", required_argument, 0, 5},
        {"eps", required_argument, 0, 6},
        {"minpts", required_argument, 0, 7},
        {"w", required_argument, 0, 8},
        {"silent", no_argument, 0, 9},
        {"cite", no_argument, 0, 22},
        {"version", no_argument, 0, 23},
        {"help", no_argument, 0, 24},
        {0, 0, 0, 0}
    };

    /** assign parameters to long options */
    while ((c = getopt_long(argc, argv, "1:2:3:4:5:6:7:8:9 10: 22 23 24", long_options, NULL)) != -1)
    {
        switch(c)
        {
            case 1:
                arg->dataInFileName = optarg;
                break;
            case 2:
                arg->dataOutFileName = optarg;
                break;
            case 3:
                arg->clusterOutFileName = optarg;
                break;
            case 4:
                arg->centerOutFileName = optarg;
                break;
            case 5:
                arg->uniqueOutFileName = optarg;
                break;
            case 6:
                arg->eps = atof(optarg);
                break;
            case 7:
                arg->minPts = atoi(optarg);
                break;
            case 8:
                arg->w = atoi(optarg);
                break;
            case 9:
                ++ silent;
                break;
            case 10:
                arg->outPathName = optarg;
                break;
            case 22:
                print_citation();
                exit(0);
            case 23:
				print_version();
				print_license();
                exit(0);
            case 24:
                fprintf(stderr, "%s", usage);
				print_license();
                exit(0);
            default:
                fprintf(stderr, "%s", usage);
				print_license();
                exit(1);
        }
    }

    if (! silent) {
		print_header();
		fflush(stdout);
		print_license();
		fflush(stdout);
	}

    return 0;
}


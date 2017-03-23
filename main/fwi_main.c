//TODO put copyright

#include "fwi/fwi_core.h"

int main(int argc, char* argv[])
{
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <params_file> <frequency_file>\n", argv[0]);
        exit(0);
    }

    double tstart, tend;
    tstart = dtime();

    int err = execute_simulation(argc, argv);

    tend = dtime() - tstart;

    fprintf(stderr, "FWI Program finished in %lf seconds\n", tend);

    return err;
}

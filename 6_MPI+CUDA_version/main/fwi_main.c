//TODO put copyright

#include "fwi_core.h"

int main(int argc, char* argv[])
{
    double tstart, tend;
    tstart = dtime();

    int err = execute_simulation(argc, argv);

    tend = dtime() - tstart;

    fprintf(stderr, "FWI Program finished in %lf seconds\n", tend);

    return err;
}

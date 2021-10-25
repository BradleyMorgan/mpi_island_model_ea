//
//  main.cpp
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 2/26/20.
//  Copyright © 2020 Bradley Morgan. All rights reserved.
//

// EA -> Objective -> Population -> Genome

#include <mpi.h>
#include <algorithm>
#include <sys/time.h>
#include <ctime>
#include <string.h>
#include "config.h"
#include "utility.h"
#include "ea.h"
#include "solver.h"
#include "meta.h"

#pragma mark FUNCTION: main()

// core initialization, overarching EA logic:
//
// this implementation is a parallelized (island model) meta ea
// evolving solutions for a benchmark optimization problem (solver)
// and a population of communication topologies, measured in duration (meta)
//
// note: functional dependency and potential fitness influence of the meta ea as it is
// applied to the communication (migration) pattern employed by the parallel solver ea
//
// note: only the solver ea employs island model (parallel) evolution
//

int main(int argc, const char * argv[]) {
    
    solver solver;
    meta meta(solver);
    
    if(config::ea_mode > 0) {
        meta_begin(meta, solver);
    } else {
        benchmark_topology(meta);
        solver_begin(meta, solver, meta.topologies.population[0]);
    }
    
    solver.ea::end(solver.solutions);
    meta.ea::end(meta.topologies);
    
}

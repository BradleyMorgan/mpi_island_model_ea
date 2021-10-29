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
#include "ea_solver.h"
#include "ea_meta.h"

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
    
    ea_solver solver;
    ea_meta meta(solver);
    
    ea_begin(meta, solver, config::ea_mode);

    //    if(config::ea_mode > 0) {
//        meta_begin(meta, solver, 1);
//    } else {
//        meta_begin(meta, solver, 0)
//        solver_begin(meta, solver, meta.topologies.population[0], 30, solver.solutions.max_evo_cycles);
//    }
    
    solver.ea::end(solver.solutions);
    meta.ea::end(meta.topologies);
    
}

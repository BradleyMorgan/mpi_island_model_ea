//
//  ea.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 3/30/20.
//  Copyright © 2020 Bradley Morgan. All rights reserved.
//

#ifndef ea_h
#define ea_h

#include <sys/time.h>
#include "dtype_heirarchy.h"
#include "dtype_ea.h"
#include "dtype_solver.h"
#include "dtype_meta.h"
#include "island.h"
#include "topology.h"
#include "stats.h"

template<typename genome> genome parent(objective<genome> &o) {
    
    genome p;

    o.cpd();
    
    // implementation uses the single armed roulette wheel approach to select
    // an individual from the population
    
    int i = 1;
    
    // random double uniformly distributed between 0 and 1
    
    double r = ((double)rand()/(double)RAND_MAX);
    
    // spin the wheel
    
    while (o.aggregate.value.cpd[i] < r ) { i++; }
    
    p = o.population[i];
    
    return p;
        
}

//void solver_init(ea &solver) {
//
//    solver.solutions.mu = config::ea_1_mu;
//    solver.solutions.lambda = config::ea_1_lambda;
//    solver.solutions.max_runs = config::ea_1_runs;
//    solver.solutions.mutation_rate = config::ea_1_mutation_rate;
//    solver.solutions.max_evo_evals = config::ea_1_max_evo_evals;
//    solver.solutions.max_fit_evals = config::ea_1_max_fit_evals;
//    
//    if(solver.variant.isle.id == 0) {
//        solver.offsets = generate_offsets(-2.5, 2.5, .5);
//    }
//    
//    MPI_Bcast(&solver.offsets, DIM, MPI_DOUBLE, 0, solver.variant.tcomm);
//    
//    // collect the time consumed by all islands in this initialization ...
//
//    // TODO: the init time gather segfaults on higher core count runs, so may need to debug at some point, but currently the init duration is somewhat insigificant
//    // MPI_Gather(&local_init_duration, 1, MPI_DOUBLE, &multi.run.stats.init_duration, 1, MPI_DOUBLE, 0, multi.meta.isle.tcomm);
//    
//    LOG(2, solver.variant.isle.id, 0, "world size: %d\r\n", solver.variant.islands);
//    LOG(2, solver.variant.isle.id, 0, "mu mode: %d\r\n", config::mu_mode);
//    LOG(2, solver.variant.isle.id, 0, "global mu %s %d\r\n", config::mu_msg, config::ea_1_mu);
//    LOG(2, solver.variant.isle.id, 0, "island mu %s %d\r\n", config::subpop_msg, config::island_mu);
//    LOG(2, solver.variant.isle.id, 0, "island lambda %s %d\r\n\r\n\r\n", config::lambda_msg, stoi(config::items["island_lambda"]));
//    
//    double init_end = MPI_Wtime();
//    
//    solver.init_duration = (init_end - solver.init_start);
//    
//}
//
//void meta_init(ea &meta, ea &solver) {
//    
//    meta.variant = solver.variant;
//    
//    meta.start = MPI_Wtime();
//    meta.init_start = MPI_Wtime();
//    meta.variant.start = MPI_Wtime();
//    meta.run.id = 1;
//    meta.run.eval.id = 1;
//    meta.run.stats.init();
//    meta.run.eval.stats.init();
//    
//    meta.topologies.mu = config::ea_2_mu;
//    meta.topologies.lambda = config::ea_2_lambda;
//    meta.topologies.max_runs = config::ea_2_runs;
//    meta.topologies.mutation_rate = config::ea_2_mutation_rate;
//    meta.topologies.max_evo_evals = config::ea_2_max_evo_evals;
//    meta.topologies.max_fit_evals = config::ea_2_max_fit_evals;
//    meta.topologies.islands = meta.variant.islands;
//    
//    meta.init_duration += (MPI_Wtime() - meta.init_start);
//    
//}

void ea_end(ea &solver, objective<solution> &obj) {
    
    LOG(5, 0, 0, "--- END EA (island %d)\r\n", solver.variant.isle.id);
    
}

void ea_end(ea &meta, objective<topology> &obj) {

    LOG(5, 0, 0, "--- END EA (island %d)\r\n", meta.variant.isle.id);

    if(meta.variant.isle.id == 0) {

        char canary[128];

        sprintf(canary, "%s/end.txt", config::logs_subpath);
        FILE *eaend = fopen(canary, "w");

        fprintf(eaend, "ended at %lu", time(0));

        fclose(eaend);

    }

    MPI_Finalize();
   
}


#endif /* ea_h */

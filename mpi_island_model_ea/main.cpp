//
//  main.cpp
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 2/26/20.
//  Copyright © 2020 Bradley Morgan. All rights reserved.
//

#include <iostream>
#include <mpi.h>
#include <array>
#include <cmath>
#include <sys/time.h>
#include <vector>
#include <fstream>
#include <ctime>
#include <time.h>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include "config.h"
#include "utility.h"
#include "dtype_heirarchy.h"
#include "dtype_solution.h"
#include "dtype_island.h"
#include "dtype_topology.h"
#include "ea.h"

#pragma mark FUNCTION: main()

// core initialization, overarching EA logic:
// this implementation is a parallelized (island model) multi-objective ea
// evolving solutions for a minimization problem (offset rastrigin, objective 1)
// and a population of communication topologies, measured in duration (objective 2)
//
// note: functional dependency and potential fitness influence of objective 1 as it is
// applied to the communication (migration) pattern employed by the island model objective 2
//
// note: only objective (1) employs island model (parallel) evolution
//

int main(int argc, const char * argv[]) {
    
    // perform mpi and global variable initialization ...
    
    ea multi = ea_init();
    
    // begin experimental runs ...
    
    for(multi.run.id = 1; multi.run.id <= config::runs; multi.run.id++) {  // n (param) experimental runs

        // clear any current populations and reset counters ...
        
        multi.run_init();
        
        // generate random objective (solution) population with fitness calculated
        // and distribute (mu/n) subpopulations to each island
        
        multi.populate(multi.solutions, solution_populate);
        multi.distribute(multi.solutions, solution_scatter);

        if(config::ea_mode > 0) { // multi-objective, evolve topologies

            if(multi.run.id == 1) {
                multi.populate(multi.topologies, topologies_populate);
            }

            // evaluate initial topology population by applying each topology and using it for n solution evals

            for(int i=0; i<multi.topologies.mu; i++) {
                multi.evaluate(multi.topologies, multi.topologies.population[i], topology_evaluate);
            }

            // perform the remaining solution evolution cycles indirectly through topology evolution
            // topology evolution also triggers evolutionary cycles of the solution population in order
            // to determine topology fitness

            while(multi.solutions.eval <= config::evals) {
                multi.evolve(multi.solutions, topology_evolve);
            }

        } else {  // benchmark, use static ring topology

            // in the benchmark case, we only need to apply the topology once ...

            if(multi.run.id == 1) {
                benchmark_topology(multi);
                multi.topologies.population[0].apply(multi.meta.isle, multi.topologies.population[0]);
            }

            // perform solver evolution ...

            while(multi.solutions.eval <= config::evals) {
                multi.coevolve(multi.solutions, solutions_evolve, multi.topologies.population[0]);
            }

        }

        // perform any value resets and write the run statistics to file output ...
        
        multi.run_end();
        
    } // run end

    multi.ea_end();


//    logic test ...
    
//    int o1run = 1;
//    int o2run = 1;
//
//    while(o1run <= config::objective_1_runs && o2run <= config::objective_2_runs) {
//
//        // objective 1 (o1): topology
//        // generate initial topology population
//
//        // objective 2 (o2): solver
//        // generate initial solution population
//
//        // T1:run=1
//
//        int t = 0; // topology index
//
//        int o1_fit_eval = 1;
//        int o1eval = 1;
//        int o2eval = 1;
//
//        // (o1) apply topology t
//
//        while(t <= config::topo_mu) {
//
//            // (o1) calculate t.fitness as migration time from (o2)
//            // (o2) evolution cycle
//
//            // T2:t=7
//
//            for(int i=0; i<=config::lambda; i++) {
//                // (o2) parent selection
//                // (o2) crossover
//                // (o2) mutation
//            }
//
//            // (o2) survival selection
//            // (o2) migration
//            // (o1) fitness += migration time
//
//            o1_fit_eval++;
//            o2eval++;
//
//            if(o1_fit_eval%config::objective_1_max_fit_evals == 0) {
//                t++;
//                // (o1) apply topology t
//            }
//
//            if(o2eval >= config::objective_2_max_evo_evals) {
//                o2run++;
//                // (o2) init run
//            }
//
//            // T3:o1eval=100
//            // T4:o2eval=200
//
//        }
//
//        // (o1) initial population fitness calculated
//
//        while (o1eval <= config::objective_1_max_evo_evals) {
//
//            // (o1) evolution cycle
//
//            for(int i=0; i<=config::topo_lambda; i++) {
//                // (o1) parent selection
//                // (o1) crossover
//                // (o1) mutation
//            }
//
//            // (o1) calculate child fitness -> (o2)
//            // (o1) survival selection
//
//            o1eval++;
//
//        }
//
//        o2run++;
//
//    }
    
}

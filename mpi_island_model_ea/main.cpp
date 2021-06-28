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
#include <algorithm>
#include "config.h"
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
    
    for(multi.run.id = 1; multi.run.id <= config::runs; multi.run.id++) {

        // clear any current populations and reset counters ...
        
        multi.run_init();
        
        // generate random objective (solution) population with fitness calculated
        
        multi.populate(solution_populate);
        
        MPI_Barrier(multi.meta.tcomm);
        
        // assign (mu/n) subpopulations to each island
        
        multi.distribute(solution_scatter, multi.solutions);
        
        MPI_Barrier(multi.meta.tcomm);
        
        if(config::ea_mode > 0) { // multi-objective, evolve topologies
    
            multi.populate(topologies_populate);
        
            // evaluate initial topology population by applying each topology and using it for n solution evals
            
            for(int i=0; i<multi.topologies.mu; i++) {
                multi.evaluate(topology_evaluate, multi.topologies.population[i]);
            }
            
            // perform the remaining solution evolution cycles indirectly through topology evolution
            // topology evolution also triggers evolutionary cycles of the solution population in order
            // to determine topology fitness
            
            while(multi.solutions.eval_id <= config::evals) {
                multi.evolve(topology_evolve);
            }
                  
        } else {  // benchmark, use static ring topology
            
            // in the benchmark case, we only need to apply the topology once ...
            
            if(multi.run.id == 1) {
                benchmark_topology(multi);
                multi.topologies.population[0].apply(multi.meta.isle, multi.topologies.population[0]);
            }
            
            // perform solver evolution ...
            
            while(multi.solutions.eval_id <= config::evals) {
                multi.evolve(solutions_evolve, multi.topologies.population[0]);
                //multi.evaluate(solution_evaluate, multi.topologies.population[0]);
            }
            
        }
        
        // perform any value resets and write the run statistics to file output ...
        
        multi.run_end();
        
    }

    MPI_Finalize();
         
 }

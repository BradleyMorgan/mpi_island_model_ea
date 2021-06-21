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
        
    ea multi = ea_init();
    
    LOG(6, 0, 0, "beginning experimental runs ...\r\n");
    
    for(multi.run.id = 1; multi.run.id <= config::runs; multi.run.id++) {

        multi.run_init();
        
        // generate random objective (solution) population with fitness calculated
        
        multi.populate(solution_populate, multi.solutions);
        
        // assign (mu/n) subpopulations to each island
        
        multi.distribute(solution_scatter, multi.solutions);
        
        if(config::ea_mode > 0) { // multi-objective, evolve topologies
    
            multi.populate(topologies_populate, multi.topologies);
        
            // evaluate initial topology population by applying each topology and using it for n solution evals
            
            for(int i=0; i<multi.topologies.mu; i++) {
                multi.evaluate(topology_evaluate, multi.topologies, multi.topologies.population[i]);
            }
            
            // perform the remaining solution evolution cycles indirectly through topology evolution
            
            while(multi.solutions.eval_id <= config::evals) {
                topology_evolve(multi);
            }
                  
        } else {  // benchmark, use static ring topology
            
            benchmark_topology(multi);
            
            multi.topologies.population[0].apply(multi.meta.isle, multi.topologies.population[0]);
            
            for(int i=1; i <= config::evals; i++) {
                multi.evaluate(solution_evaluate, multi.solutions, multi.topologies.population[0]);
            }
            
        }
        
        multi.run_end();
        
    }

    MPI_Finalize();
         
 }

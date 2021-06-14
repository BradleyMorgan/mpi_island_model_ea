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

    multi.offsets = generate_offsets(-2.5, 2.5, .5);
    
    for(multi.run.id = 1; multi.run.id <= config::runs; multi.run.id++) {

        LOG(4, multi.meta.isle.id, 0, "initializing objective (solution) population ...\r\n");

        multi.solutions.population.clear();
        multi.meta.isle.population.clear();
        
        multi.eval.stats.init();
        multi.solutions.eval_id = 0;
        multi.topologies.eval_id = 0;
        
        multi.populate(solution_populate, multi.solutions);
        multi.distribute(solution_scatter, multi.solutions);
        
        LOG(4, multi.meta.isle.id, 0, "initialized objective (solution) population: total fitness = %f\r\n", multi.solutions.total_fitness);
        
        LOG(4, multi.meta.isle.id, 0, "initializing objective (topology) population ...\r\n");

        multi.populate(topologies_populate, multi.topologies);
        
        LOG(4, multi.meta.isle.id, 0, "initialized objective (topology) population: %lu\r\n", multi.topologies.population.size());
        LOG(4, multi.meta.isle.id, 0, "total channels = %d\r\n",  multi.topologies.population[0].channel_count);

        // evaluate initial topology population by applying each topology and using it for n solution evals ...
        
        for(int i=0; i<multi.topologies.mu; i++) {

            multi.evaluate(topology_evaluate, multi.topologies, i);

        }
        
        // perform the remaining solution evolution cycles indirectly through topology evolution ...
        
        while(multi.solutions.eval_id <= config::evals) {
        
            topology_evolve(multi);
            
        }
        
    }

    MPI_Finalize();
         
 }

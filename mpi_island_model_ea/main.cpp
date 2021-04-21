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

    //for(multi.run.id = 1; multi.run.id <= config::runs; multi.run.id++) {

    LOG(4, multi.meta.isle.id, 0, "initializing objective (solution) population ...\r\n");

    multi.populate(solution_populate, multi.solutions);
    multi.distribute(solution_scatter, multi.solutions);
    multi.offsets = generate_offsets(-2.5, 2.5, .5);
    
    LOG(4, multi.meta.isle.id, 0, "initialized objective (solution) population: total fitness = %f\r\n", multi.solutions.total_fitness);
    
    LOG(4, multi.meta.isle.id, 0, "initializing objective (topology) population ...\r\n");

    multi.populate(topologies_populate, multi.topologies);

    LOG(4, multi.meta.isle.id, 0, "initialized objective (topology) population:");
    LOG(4, multi.meta.isle.id, 0, "total channels = %d\r\n",  multi.topologies.population[0].channel_count);

    for(int i=0; i<multi.topologies.mu; i++) {
        
        multi.evaluate(topology_evaluate, multi.topologies, i);
        
    }
    
        //topology_evaluate(multi, multi.topologies.population[0]);
        

    
//        for(int i = 1; i <= config::evals; i++) {
//
//            //std::vector<topology> children = topo_gen(multi);
//
//            multi.topologies.children = topo_gen(multi);
//
//            for( int k=0; k<config::topo_lambda; k++ ) {
//
//                topology t = multi.topologies.children[k];
//
//                LOG(3, multi.meta.isle.id, 0, "evaluating child topology %d comm size %lu\r\n", t.id, t.island_comms.size());
//
//                // rank authoritative source parses and queues the distribution of
//                // island-specific migration channels to be assigned
//
//                LOG(3, multi.meta.isle.id, 0, "distributing topology %d\r\n", t.id);
//
//                topology_distribute(t, multi);
//
//                // commit the queued distribution of island senders and receivers
//
//                LOG(3, 0, 0, "rank %d applying topology %d\r\n", multi.meta.isle.id, t.id);
//
//                topology_apply(t, multi);
//
//                topology_eval(t, multi);
//
//                //multi.topologies.population.push_back(t);
//
//                if(multi.meta.isle.id == 0) {
//                    if(i == 0) {
//                        multi.eval.stats.best_topology = t;
//                    } else if(t.fitness > multi.eval.stats.best_topology.fitness) {
//                        multi.eval.stats.best_topology = t;
//                    }
//
//                }
//
//            }
//
//            if(multi.meta.isle.id == 0) {
//
//                multi.topologies.population.insert(multi.topologies.population.end(), multi.topologies.children.begin(), multi.topologies.children.end());
//
//                std::sort(multi.topologies.population.begin(), multi.topologies.population.end(), compare_topo_fitness);
//                std::reverse(multi.topologies.population.begin(), multi.topologies.population.end());
//
//                multi.topologies.population.erase(multi.topologies.population.begin()+config::topo_mu, multi.topologies.population.end());
//
//                LOG(4, multi.meta.isle.id, 0, "truncating topologies size %lu at %d...\r\n", multi.topologies.population.size(), multi.eval.id);
//
//                std::vector<topology>::iterator min = std::min_element(multi.topologies.population.begin(), multi.topologies.population.end(), compare_topo_fitness);
//                std::replace_if(multi.topologies.population.begin(), multi.topologies.population.end(), is_zero, *min);
//
//                LOG(6, multi.meta.isle.id, 0, "replacing 0 fit with %2.10f\r\n", min->fitness);
//
//                std::sort(multi.topologies.population.begin(), multi.topologies.population.end(), compare_topo_fitness);
//                std::reverse(multi.topologies.population.begin(), multi.topologies.population.end());
//
//                multi.topologies.population.erase(multi.topologies.population.begin()+config::topo_mu, multi.topologies.population.end());
//
//                LOG(3, multi.meta.isle.id, 0, "removed %d topologies, new size = %lu\r\n", config::topo_lambda, multi.topologies.population.size());
//                LOG(10, multi.meta.isle.id, 0, "topology 0 fitness = %2.10f\r\n", multi.topologies.population[0].fitness);
//
//            }
//
    //        for(int i=0; i<config::topo_mu; i++) {
    //            topology_eval(multi.topologies.population[i], multi);
    //        }
            
            
            // for(multi.eval.id = 1; multi.eval.id<=config::evals; multi.eval.id++) {
            // topologies must be evaluated for fitness, apply to (1) island migration pattern
            // and measure fitness in communication time
            
            //for(multi.eval.id = 1; multi.eval.id <= config::runs; multi.eval.id++) {
                
    //            for(int i=0; i<config::topo_mu; i++) {
    //
    //                // rank authoritative source parses and queues the distribution of
    //                // island-specific migration channels to be assigned
    //
    //                topology_distribute(i, multi);
    //
    //                // commit the queued distribution of island senders and receivers
    //
    //                topology_apply(i, multi);
    //
    //                topology_eval(i, multi);
    //
    //             }
    
            // }
        

            
    MPI_Finalize();
         
 }

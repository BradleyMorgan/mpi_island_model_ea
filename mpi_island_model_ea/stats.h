//
//  stats.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 3/30/20.
//  Copyright © 2020 Bradley Morgan. All rights reserved.
//

#ifndef stats_h
#define stats_h

#include <mpi.h>

struct rstats {
    
    individual solution;
    
    double init_duration;
    
};

struct estats {
    
    double init_duration = 0.0;
    double eval_start = 0.0;
    double eval_duration = 0.0;
    double total_scatter_time = 0.0;
    double total_gather_time = 0.0;
    double total_migrate_time = 0.0;
    double topo_migrate_time = 0.0;
    double local_best_fitness = 0.0;
    double global_best_fitness = 0.0;
    double average_local_best_fitness = 0.0;
    double average_global_best_fitness = 0.0;
    
    double local_best_topo_fitness = 0.0;
    double global_best_topo_fitness = 0.0;
    double average_local_best_topo_fitness = 0.0;
    double average_global_best_topo_fitness = 0.0;
    
    topology best_topology;
    
    std::vector<double> average_local_best_fitnesses;
    std::vector<double> average_global_best_fitnesses;
    std::vector<double> average_local_best_topo_fitnesses;
    std::vector<double> average_global_best_topo_fitnesses;
    
    individual solution;
    
};

void log_fn_eval_stats(std::vector<individual> &population, std::vector<topology> &topologies, int &run, int &eval, estats &eval_stats, rstats &run_stats) {
    
    eval_stats.local_best_fitness = population[0].fitness;
    eval_stats.average_local_best_fitnesses.push_back(population[0].fitness);
    
    eval_stats.local_best_topo_fitness = topologies[0].fitness;
    eval_stats.average_local_best_topo_fitnesses.push_back(topologies[0].fitness);
    
    if(population[0].fitness > eval_stats.global_best_fitness) {
        eval_stats.solution = population[0];
        eval_stats.average_global_best_fitnesses.push_back(population[0].fitness);
        eval_stats.global_best_fitness = population[0].fitness;
    }
    
    if(topologies[0].fitness > eval_stats.global_best_topo_fitness || eval_stats.global_best_topo_fitness == 0.0) {
        eval_stats.best_topology = topologies[0];
        eval_stats.average_global_best_topo_fitnesses.push_back(topologies[0].fitness);
        eval_stats.global_best_topo_fitness = topologies[0].fitness;
    }
    
    double total_fitness = 0.0;
    double total_topo_fitness = 0.0;
    
    for(int i=0; i<population.size(); i++) {
        total_fitness += population[i].fitness;
    }
    
    for(int i=0; i<topologies.size(); i++) {
        total_topo_fitness += topologies[i].fitness;
    }
    
    LOG(6, 0, 0, "total fit: %2.10f\r\n", total_fitness);
    
    eval_stats.average_local_best_fitness = std::accumulate(eval_stats.average_local_best_fitnesses.begin(), eval_stats.average_local_best_fitnesses.end(), 0.0) / eval_stats.average_local_best_fitnesses.size();
    eval_stats.average_global_best_fitness = std::accumulate(eval_stats.average_global_best_fitnesses.begin(), eval_stats.average_global_best_fitnesses.end(), 0.0) / eval_stats.average_global_best_fitnesses.size();
    eval_stats.average_local_best_topo_fitness = std::accumulate(eval_stats.average_local_best_topo_fitnesses.begin(), eval_stats.average_local_best_topo_fitnesses.end(), 0.0) / eval_stats.average_local_best_topo_fitnesses.size();
    eval_stats.average_global_best_topo_fitness = std::accumulate(eval_stats.average_global_best_topo_fitnesses.begin(), eval_stats.average_global_best_topo_fitnesses.end(), 0.0) / eval_stats.average_global_best_topo_fitnesses.size();
    
    double average_fitness = total_fitness / population.size();
    double average_gather_time = eval_stats.total_gather_time / eval;
    double average_scatter_time = eval_stats.total_scatter_time / eval;
    double average_migrate_time = eval_stats.total_migrate_time / eval;
    double average_topo_fitness = total_topo_fitness / topologies.size();
    
    double eval_duration = ( std::clock() - eval_stats.eval_start ) / (double) CLOCKS_PER_SEC;
    
    std::fprintf(config::stats_out, "%d,%d,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%d\r\n", run, eval, average_fitness, eval_stats.local_best_fitness, eval_stats.global_best_fitness, eval_stats.average_local_best_fitness, eval_stats.average_global_best_fitness, average_scatter_time, average_gather_time, average_migrate_time, eval_stats.init_duration, eval_duration, average_topo_fitness, eval_stats.local_best_topo_fitness, eval_stats.global_best_topo_fitness, eval_stats.average_local_best_topo_fitness, eval_stats.average_global_best_topo_fitness, topologies[0].round_fitness, topologies[0].rounds);
    
    LOG(2, 0, 0, "%d,%d,%.10f,%.10f,%.10f,%.10f,%.10f,%2.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%d\r\n", run, eval, average_fitness, eval_stats.local_best_fitness, eval_stats.global_best_fitness, eval_stats.average_local_best_fitness, eval_stats.average_global_best_fitness, average_scatter_time, average_gather_time, average_migrate_time, run_stats.init_duration, eval_duration, average_topo_fitness, eval_stats.local_best_topo_fitness, eval_stats.global_best_topo_fitness, eval_stats.average_local_best_topo_fitness, eval_stats.average_global_best_topo_fitness, topologies[0].round_fitness, topologies[0].rounds);
    
}

#endif /* stats_h */

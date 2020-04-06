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

double total_scatter_time = 0.0;
double total_gather_time = 0.0;
double total_migrate_time = 0.0;
double global_best_fitness = 0.0;
double average_local_best_fitness = 0.0;
double average_global_best_fitness = 0.0;

std::vector<double> average_local_best_fitnesses;
std::vector<double> average_global_best_fitnesses;

void init_stats() {
    
    total_scatter_time = 0.0;
    total_gather_time = 0.0;
    total_migrate_time = 0.0;
    global_best_fitness = 0.0;
    average_local_best_fitness = 0.0;
    average_global_best_fitness = 0.0;
    
}

void log_fn_eval_stats(std::vector<individual> &population, int &run, int &eval, individual &solution, MPI_Wtime &eval_start) {
    
    double local_best_fitness = population[0].fitness;
    average_local_best_fitnesses.push_back(population[0].fitness);
    
    if(population[0].fitness > global_best_fitness) {
        solution = population[0];
        average_global_best_fitnesses.push_back(population[0].fitness);
        global_best_fitness = population[0].fitness;
    }
    
    double total_fitness = 0.0;
    
    for(int i=0; i<population.size(); i++) {
        total_fitness += population[i].fitness;
    }
    
    average_local_best_fitness = std::accumulate(average_local_best_fitnesses.begin(), average_local_best_fitnesses.end(), 0.0) / average_local_best_fitnesses.size();
    average_global_best_fitness = std::accumulate(average_global_best_fitnesses.begin(), average_global_best_fitnesses.end(), 0.0) / average_global_best_fitnesses.size();
    
    double average_fitness = total_fitness / population.size();
    double average_gather_time = total_gather_time / eval;
    double average_scatter_time = total_scatter_time / eval;
    double average_migrate_time = total_migrate_time / eval;
    
    double eval_duration = ( std::clock() - eval_start ) / (double) CLOCKS_PER_SEC;
    
    std::fprintf(config::stats_out, "%d,%d,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f\r\n", run, eval, average_fitness, local_best_fitness, global_best_fitness, average_local_best_fitness, average_global_best_fitness, average_scatter_time, average_gather_time, average_migrate_time, init_duration, eval_duration);
    
    LOG(4, "%d,%d,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f\r\n", run, eval, average_fitness, local_best_fitness, global_best_fitness, average_local_best_fitness, average_global_best_fitness, average_scatter_time, average_gather_time, average_migrate_time, init_duration, eval_duration);
    
}

#endif /* stats_h */

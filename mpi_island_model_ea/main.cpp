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
#include "evolution.h"
#include "config.h"

// populate an initial population with random inputs ...

std::vector<individual> initial_population(std::array<double, DIM> &offsets) {
    
    std::vector<individual> population;
    
    for(int i=0; i<config::mu; i++) {
        
        individual p;
        
        for (int j = 0; j < DIM; j++) {
            p.input[j] = drand(-5.12, 5.12);
        }
        
        p.fitness = offset_rastrigin(p.input, offsets);
        
        population.push_back(p);
        
    }
    
    return population;
    
}

int main(int argc, const char * argv[]) {

    // initialize timer
    
    std::clock_t start;
    
    double init_duration;
    double eval_duration;

    start = std::clock();
    
    // initialize MPI environment ...
    
    MPI_Init(NULL, NULL);
    
    init_duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // load configuration items ...
    
    config::load("config.txt", world_size, world_rank);
    
    int subpopulation_size = config::mu / world_size;
    
    MPI_Datatype individual_type;
    
    // MPI derived datatype for population individual ...
    
    int lengths[4] = { DIM, 1, 1, 1};
    MPI_Aint displacements[4] = { 0, sizeof(double)*DIM, sizeof(double)*(DIM+1) };
    MPI_Datatype types[4] = { MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE };
    MPI_Type_create_struct(4, lengths, displacements, types, &individual_type);
    MPI_Type_commit(&individual_type);
    
    // evolve the populations ...
    
    for(int run=1; run<=config::runs; run++) {
        
        double run_start = MPI_Wtime();
        
        double total_scatter_time = 0.0;
        double total_gather_time = 0.0;
        double total_migrate_time = 0.0;
        double global_best_fitness = 0.0;
        double average_local_best_fitness = 0.0;
        double average_global_best_fitness = 0.0;
        
        std::vector<double> average_local_best_fitnesses;
        std::vector<double> average_global_best_fitnesses;
        
        // each MPI process will maintain its own population, so we define an island for each ...
        
        std::vector<individual> population;
        
        island isle;
        isle.id = world_rank;
        isle.population.resize(subpopulation_size);
        
        std::array<double, DIM> offsets = generate_offsets(-2.5, 2.5, .5);
        
        // only the root process will create the full initial population ...
        
        if(world_rank == 0) {
            
            population = initial_population(offsets);
            
            std::sort(population.begin(), population.end(), compare_fitness);
            std::reverse(population.begin(), population.end());
            
            global_best_fitness = population[0].fitness;
            
        }
        
        // separate the single full population from the root process to subpopulations across all processes ...
        
        double scatter_start = MPI_Wtime();
        MPI_Scatter(&population[0], subpopulation_size, individual_type, &isle.population[0], subpopulation_size, individual_type, 0, MPI_COMM_WORLD);
        double scatter_end = MPI_Wtime();
        double scatter_time = scatter_end - scatter_start;
        
        total_scatter_time += scatter_time;
        
        printf("rank %d starting population size is %lu\r\n", world_rank, isle.population.size());
        
        // build a topology by assigning send \ receive neighbors ...
        
        create_topology(isle, world_size);
    
        // begin evolution ...
        
        for(int eval=1; eval<=config::evals; eval++) {
            
            double eval_start = std::clock();
        
            isle.calc_cpd();
            
            std::vector<individual> children = crossover(isle, offsets);
            
            select_survivors(isle, children, config::mu/world_size);

            double migrate_start = MPI_Wtime();
            isle.send_migrant();
            isle.receive_migrant();
            double migrate_end = MPI_Wtime();
            double migrate_time = migrate_end - migrate_start;
            
            total_migrate_time += migrate_time;
            
            population.clear();
            population.resize(config::mu);
            
            double gather_start = MPI_Wtime();
            MPI_Gather(&isle.population[0], subpopulation_size, individual_type, &population[0], subpopulation_size, individual_type, 0, MPI_COMM_WORLD);
            double gather_end = MPI_Wtime();
            double gather_time = gather_start - gather_end;
            
            total_gather_time += gather_time;
            
            if(world_rank == 0 && eval % 100 == 0) {
                
                std::sort(population.begin(), population.end(), compare_fitness);
                std::reverse(population.begin(), population.end());
                
                double local_best_fitness = population[0].fitness;
                average_local_best_fitnesses.push_back(population[0].fitness);
                
                if(population[0].fitness > global_best_fitness) {
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
                
                eval_duration = ( std::clock() - eval_start ) / (double) CLOCKS_PER_SEC;
                
                std::fprintf(config::stats_out, "%d,%d,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f\r\n", run, eval, average_fitness, local_best_fitness, global_best_fitness, average_local_best_fitness, average_global_best_fitness, average_scatter_time, average_gather_time, average_migrate_time, init_duration, eval_duration);
                
            }
            
            if(eval % 100 == 0) {
                
                printf("rank %d eval %d average population is %2.8f\r\n", world_rank, eval, isle.average_fitness());
                
            }
            
        }
        
        if(world_rank == 0) {
            
            double run_end = MPI_Wtime();
            
            std::fprintf(config::run_stats_out, "%d,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%d,%d\r\n", run, global_best_fitness, average_local_best_fitness, average_global_best_fitness, scatter_time, run_end - run_start, init_duration, world_size, subpopulation_size);
            
        }
        
        fflush(config::run_stats_out);
        fflush(config::stats_out);
        
    }

    fclose(config::stats_out);
    
    MPI_Finalize();
    
}

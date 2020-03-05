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

    // initialize MPI environment ...
    
    MPI_Init(NULL, NULL);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // load configuration items ...
    
    config::load("config.txt", world_size, world_rank);
    
    int num_islands = config::mu / world_size;
    
    MPI_Datatype individual_type;
    
    // MPI derived datatype for population individual ...
    
    int lengths[4] = { DIM, 1, 1, 1};
    MPI_Aint displacements[4] = { 0, sizeof(double)*DIM, sizeof(double)*(DIM+1) };
    MPI_Datatype types[4] = { MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE };
    MPI_Type_create_struct(4, lengths, displacements, types, &individual_type);
    MPI_Type_commit(&individual_type);
    
    // evolve the populations ...
    
    for(int run=1; run<=config::runs; run++) {

        double total_scatter_time = 0.0;
        double total_gather_time = 0.0;
        double total_migrate_time = 0.0;
        
        // each MPI process will maintain its own population, so we define an island for each ...
        
        std::vector<individual> population;
        
        island isle;
        isle.id = world_rank;
        
        isle.population.resize(num_islands);
        
        std::array<double, DIM> offsets = generate_offsets(-2.5, 2.5, .5);
        
        // only the root process will create the full initial population ...
        
        if(world_rank == 0) { population = initial_population(offsets); }
        
        // separate the single full population from the root process to subpopulations across all processes ...
        
        double scatter_start = MPI_Wtime();
        MPI_Scatter(&population[0], num_islands, individual_type, &isle.population[0], num_islands, individual_type, 0, MPI_COMM_WORLD);
        double scatter_end = MPI_Wtime();
        
        double scatter_time = scatter_end - scatter_start;
        
        total_scatter_time += scatter_time;
        
        printf("rank %d starting population size is %lu\r\n", world_rank, isle.population.size());
        
        // build a topology by assigning send \ receive neighbors ...
        
        create_topology(isle, world_size);
    
        // begin evolution ...
        
        for(int eval=1; eval<=config::evals; eval++) {
        
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
            MPI_Gather(&isle.population[0], num_islands, individual_type, &population[0], num_islands, individual_type, 0, MPI_COMM_WORLD);
            double gather_end = MPI_Wtime();
            
            double gather_time = gather_start - gather_end;
            
            total_gather_time += gather_time;
            
            if(world_rank == 0 && eval % 100 == 0) {
                
                double total_fitness = 0.0;
                
                for(int i=0; i<population.size(); i++) {
                    total_fitness += population[i].fitness;
                }
                
                double average_fitness = total_fitness / population.size();
                double average_gather_time = total_gather_time / eval;
                double average_scatter_time = total_scatter_time / eval;
                double average_migrate_time = total_migrate_time / eval;
                
                std::fprintf(config::stats_out, "%d,%d,%2.10f,%2.10f,%2.10f,%2.10f\r\n", run, eval, average_fitness, average_scatter_time, average_gather_time, average_migrate_time);
                printf("%d global average fitness %2.10f\r\n", eval, average_fitness);
                printf("%d average gather time %2.10f\r\n", eval, average_gather_time);
                printf("%d average scatter time %2.10f\r\n", eval, average_scatter_time);
                printf("%d average migrate time %2.10f\r\n", eval, average_migrate_time);
                
            }
            
            if(eval % 100 == 0) {
                
                printf("rank %d eval %d average population is %2.8f\r\n", world_rank, eval, isle.average_fitness());
                
            }
            
        }
        
        fflush(config::stats_out);
        
    }

    fclose(config::stats_out);
    
    MPI_Finalize();
    
}

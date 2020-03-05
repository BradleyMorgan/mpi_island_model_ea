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
    
    config::load("config.txt", world_size);
    
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
    
        // each MPI process will maintain its own population, so we define an island for each ...
        
        std::vector<individual> population;
        
        island isle;
        isle.id = world_rank;
        
        isle.population.resize(num_islands);
        
        std::array<double, DIM> offsets = generate_offsets(-2.5, 2.5, .5);
        
        // only the root process will create the full initial population ...
        
        if(world_rank == 0) { population = initial_population(offsets); }
        
        // separate the single full population from the root process to subpopulations across all processes ...
        
        MPI_Scatter(&population[0], num_islands, individual_type, &isle.population[0], num_islands, individual_type, 0, MPI_COMM_WORLD);
        
        printf("rank %d starting population size is %lu\r\n", world_rank, isle.population.size());
        
        // build a topology by assigning send \ receive neighbors ...
        
        create_topology(isle, world_size);
    
        // begin evolution ...
        
        for(int eval=1; eval<=config::evals; eval++) {
        
            isle.calc_cpd();
            
            std::vector<individual> children = crossover(isle, offsets);
            
            select_survivors(isle, children, config::mu/world_size);

            isle.send_migrant();
            isle.receive_migrant();
            
            population.clear();
            population.resize(config::mu);
            
            MPI_Gather(&isle.population[0], num_islands, individual_type, &population[0], num_islands, individual_type, 0, MPI_COMM_WORLD);
            
            if(world_rank == 0 && eval % 100 == 0) {
                
                double total_fitness = 0.0;
                
                for(int i=0; i<population.size(); i++) {
                    total_fitness += population[i].fitness;
                }
                
                double average_fitness = total_fitness / population.size();
                
                std::fprintf(config::stats_out, "%d %d %2.10f\r\n", run, eval, average_fitness);
                printf("%d global average fitness %2.10f\r\n", eval, average_fitness);
                
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

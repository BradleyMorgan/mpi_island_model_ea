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
#include "utility.h"

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

    individual solution;
    
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
    
    // TODO: round for non-integer sizes ...
    
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
        
        MPI_Barrier(MPI_COMM_WORLD);
        
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
        population.clear();
        population.resize(config::mu);
        
        island isle;
        isle.id = world_rank;
        isle.population.resize(subpopulation_size);
        
        std::array<double, DIM> offsets = generate_offsets(-2.5, 2.5, .5);
        
        // only the root process will create the full initial population ...
        
        std::vector<int> island_ids;
        
        for(int i=0; i<world_size; i++) {
         
            island_ids.push_back(i);
            
        }
        
        std::vector<group> topology;
        topology.clear();
        topology.resize(world_size);
        
        int send_size = 0;
        int rec_size = 0;
        
        if(world_rank == 0) {
        
            population = initial_population(offsets);
            
            LOG(10, "received %lu isle ids\r\n", island_ids.size());
            
            std::sort(population.begin(), population.end(), compare_fitness);
            std::reverse(population.begin(), population.end());
            
            global_best_fitness = population[0].fitness;
            
            topology = create_dyn_topology(&island_ids);
            
            fprintf(config::topo_out, "---- RUN %d TOPOLOGY ----\r\n\r\n", run);
            
            for(int i=0; i<island_ids.size(); i++) {
                
                for(int k=0; k<topology[i].senders.size(); k++) {
                    
                    printf("%d -> ", topology[i].senders[k]);
                    fprintf(config::topo_out, "%d -> ", topology[i].senders[k]);
                    
                }
                
                printf("[%d] -> ", topology[i].node);
                fprintf(config::topo_out, "[%d] -> ", topology[i].node);
                
                for(int k=0; k<topology[i].receivers.size(); k++) {
                    
                    printf("%d -> ", topology[i].receivers[k]);
                    fprintf(config::topo_out, "%d -> ", topology[i].receivers[k]);
                    
                }
                
                
                printf("\r\n");
                fprintf(config::topo_out, "\r\n");
                
                fflush(config::topo_out);
                
            }
            
            fprintf(config::topo_out, "\r\n------------\r\n");
            
            for(int i=0; i<world_size; i++) {
                
                LOG(10, "sending topology to %d ...\r\n", i);
                
                send_size = (int)topology[i].senders.size();
                rec_size = (int)topology[i].receivers.size();
                
                LOG(10, "sending %d senders to rank %d ...\r\n", send_size, i);
                MPI_Send(&send_size, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
                MPI_Send(&topology[i].senders[0], send_size, MPI_INT, i, 2, MPI_COMM_WORLD);
                
                LOG(10, "sending %d receivers to rank %d ...\r\n", rec_size, i);
                MPI_Send(&rec_size, 1, MPI_INT, i, 3, MPI_COMM_WORLD);
                MPI_Send(&topology[i].receivers[0], rec_size, MPI_INT, i, 4, MPI_COMM_WORLD);
                
            }
            
            
        }
        
        MPI_Recv(&send_size, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        isle.senders.resize(send_size);
        
        MPI_Recv(&isle.senders[0], send_size, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        LOG(10, "rank %d received %lu senders: ", world_rank, isle.senders.size());
        for(int i=0; i<isle.senders.size(); i++) { LOG(10, "%d ", isle.senders[i]); }
        LOG(10, "\r\n");
        
        MPI_Recv(&rec_size, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        isle.receivers.resize(rec_size);
        
        MPI_Recv(&isle.receivers[0], rec_size, MPI_INT, 0, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        LOG(10, "rank %d got %lu receivers: ", world_rank, isle.receivers.size());
        for(int i=0; i<isle.receivers.size(); i++) { LOG(10, "%d ", isle.receivers[i]); }
        LOG(10, "\r\n");

        MPI_Comm tcomm = MPI_COMM_WORLD;

        int procnum, proclen;
        char procname[MPI_MAX_PROCESSOR_NAME];

        MPI_Get_processor_name(procname, &proclen);
        
        GETCPU(procnum);
        
        LOG(4, "** Process %d in world size %d on core %d processor %s\r\n", world_rank, world_size, procnum, procname);
        
        // separate the single full population from the root process to subpopulations across all processes ...
        
        LOG(10, "scattering population ...\r\n");
        double scatter_start = MPI_Wtime();
        MPI_Scatter(&population[0], subpopulation_size, individual_type, &isle.population[0], subpopulation_size, individual_type, 0, tcomm);
        double scatter_end = MPI_Wtime();
        double scatter_time = scatter_end - scatter_start;
        LOG(10, "population scattered...\r\n");
        
        total_scatter_time += scatter_time;
        
        // begin evolution ...
        
        for(int eval=1; eval<=config::evals; eval++) {
            
            MPI_Barrier(MPI_COMM_WORLD);
            
            double eval_start = std::clock();
        
            isle.calc_cpd();
            
            std::vector<individual> children = crossover(isle, offsets);
            
            select_survivors(isle, children, config::mu/world_size);

            double migrate_start = MPI_Wtime();
            isle.send_migrant(tcomm);
            isle.receive_migrant(tcomm);
            double migrate_end = MPI_Wtime();
            double migrate_time = migrate_end - migrate_start;
            
            total_migrate_time += migrate_time;
            
            LOG(10, "gathering population, subpopulation %d size %d ...\r\n", world_rank, subpopulation_size);
            double gather_start = MPI_Wtime();
            MPI_Gather(&isle.population[0], subpopulation_size, individual_type, &population[0], subpopulation_size, individual_type, 0, tcomm);
            double gather_end = MPI_Wtime();
            double gather_time = gather_end - gather_start;
            LOG(10, "population gathered ...\r\n");
            
            total_gather_time += gather_time;
            
            if(world_rank == 0 && eval % 100 == 0) {
            
                LOG(8, "population size %lu, member = %2.10f\r\n", population.size(), population[population.size()-1].fitness);
                
                std::sort(population.begin(), population.end(), compare_fitness);
                std::reverse(population.begin(), population.end());
                
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
                
                eval_duration = ( std::clock() - eval_start ) / (double) CLOCKS_PER_SEC;
                
                std::fprintf(config::stats_out, "%d,%d,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f\r\n", run, eval, average_fitness, local_best_fitness, global_best_fitness, average_local_best_fitness, average_global_best_fitness, average_scatter_time, average_gather_time, average_migrate_time, init_duration, eval_duration);
                
                LOG(4, "%d,%d,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f\r\n", run, eval, average_fitness, local_best_fitness, global_best_fitness, average_local_best_fitness, average_global_best_fitness, average_scatter_time, average_gather_time, average_migrate_time, init_duration, eval_duration);
                
            }
                
        }
        
        if(world_rank == 0) {
            
            double run_end = MPI_Wtime();
            
            std::fprintf(config::run_stats_out, "%d,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%d,%d\r\n", run, global_best_fitness, average_local_best_fitness, average_global_best_fitness, total_scatter_time, total_gather_time, total_migrate_time, run_end - run_start, init_duration, world_size, subpopulation_size);
            
            LOG(2, "%d,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%d,%d\r\n", run, global_best_fitness, average_local_best_fitness, average_global_best_fitness, total_scatter_time, total_migrate_time, total_gather_time, run_end - run_start, init_duration, world_size, subpopulation_size);
            
        }
        
        fflush(config::run_stats_out);
        fflush(config::stats_out);
        
    }

    if(world_rank == 0) {
        
        for(int i=0; i<DIM; i++) {
            std::fprintf(config::solution_out, "%2.10f,", solution.input[i]);
        }
        
        fclose(config::stats_out);
        fclose(config::run_stats_out);
        fclose(config::log_out);
        fclose(config::solution_out);
        fclose(config::topo_out);
        
    }
    
    MPI_Finalize();
        
}

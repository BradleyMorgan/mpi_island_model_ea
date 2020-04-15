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
#include "stats.h"

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

// populate an initial population with random inputs ...

std::vector<topology> initial_topo_population(std::vector<int> &ids) {
    
    LOG(10, 0, 0, "initializing topology population for %lu islands\r\n", ids.size());
    std::vector<topology> population;
    
    for(int i=0; i<config::topo_mu; i++) {
        
        topology t;
        
        LOG(10, 0, 0, "adding topology %d\r\n", i);
        t.comm = create_dyn_topology(ids);
        population.push_back(t);
        
    }
    
    LOG(10, 0, 0, "initialized population for %lu islands\r\n", ids.size());
    return population;
    
}

int main(int argc, const char * argv[]) {

    // initialize timer
    
    std::clock_t start = std::clock();
    
    rstats run_stats;
    
    // initialize MPI environment ...
    
    MPI_Init(NULL, NULL);
    
    run_stats.init_duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

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
    
        LOG(4, world_rank, 0, "**** RUN %d ****\r\n", run);
        MPI_Barrier(MPI_COMM_WORLD);
        
        double run_start = MPI_Wtime();

        estats eval_stats;
        
        // each MPI process will maintain its own population, so we define an island for each ...
        
        std::vector<individual> population;
        population.clear();
        population.resize(config::mu);
        
        std::vector<topology> topologies;
        topologies.clear();
        topologies.resize(config::topo_mu);
        
        island isle;
        isle.id = world_rank;
        isle.population.resize(subpopulation_size);
        
        std::array<double, DIM> offsets = generate_offsets(-2.5, 2.5, .5);
        
        // only the root process will create the full initial population ...
        
        std::vector<int> island_ids;
        
        for(int i=0; i<world_size; i++) {
         
            island_ids.push_back(i);
            
        }
        
        int send_size = 0;
        int rec_size = 0;
        
        if(world_rank == 0) {
        
            population = initial_population(offsets);
            topologies = initial_topo_population(island_ids);
            
            LOG(6, 0, 0, "received %lu isle ids for %lu topologies\r\n", island_ids.size(), topologies.size());
            
            std::sort(population.begin(), population.end(), compare_fitness);
            std::reverse(population.begin(), population.end());
            
            eval_stats.global_best_fitness = population[0].fitness;
            eval_stats.best_topology = topologies[0];
            
            
        }

        MPI_Comm tcomm = MPI_COMM_WORLD;

        int procnum, proclen;
        char procname[MPI_MAX_PROCESSOR_NAME];

        MPI_Get_processor_name(procname, &proclen);
        
        GETCPU(procnum);
        
        LOG(6, 0, 0, "** Process %d in world size %d on core %d processor %s\r\n", world_rank, world_size, procnum, procname);
        
        // separate the single full population from the root process to subpopulations across all processes ...
        
        LOG(10, world_rank, 0, "scattering population ...\r\n");
        double scatter_start = MPI_Wtime();
        MPI_Scatter(&population[0], subpopulation_size, individual_type, &isle.population[0], subpopulation_size, individual_type, 0, tcomm);
        double scatter_end = MPI_Wtime();
        double scatter_time = scatter_end - scatter_start;
        LOG(10, world_rank, 0, "population scattered...\r\n");
        
        eval_stats.total_scatter_time += scatter_time;
        
        // begin evolution ...
        
        for(int eval=1; eval<=config::evals; eval++) {
            
            int t = (eval-1)%(config::topo_mu);
        
            if(world_rank == 0) {
        
                LOG(10, world_rank, 0, "T = (%d)mod(%d) = %d\r\n", eval-1, config::topo_mu, t);
                
                for(int i=0; i<world_size; i++) {
                    
                    LOG(10, 0, 0, "sending topology %d to %d ...\r\n", t, i);
                    
                    send_size = (int)topologies[t].comm[i].senders.size();
                    rec_size = (int)topologies[t].comm[i].receivers.size();
                    
                    LOG(10, 0, 0, "sending %d senders to rank %d ...\r\n", send_size, i);
                    MPI_Send(&send_size, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
                    MPI_Send(&topologies[t].comm[i].senders[0], send_size, MPI_INT, i, 2, MPI_COMM_WORLD);
                    
                    LOG(10, 0, 0, "sending %d receivers to rank %d ...\r\n", rec_size, i);
                    MPI_Send(&rec_size, 1, MPI_INT, i, 3, MPI_COMM_WORLD);
                    MPI_Send(&topologies[t].comm[i].receivers[0], rec_size, MPI_INT, i, 4, MPI_COMM_WORLD);
                    
                }
                
                if(t == 0) {
                
                    LOG(10, world_rank, 0, "adding children...\r\n");
                    std::vector<topology> children = topo_gen(topologies, world_size);
                    topologies.insert(topologies.end(), children.begin(), children.end());
                    
                    LOG(8, world_rank, 0, "added %lu children to topologies, new topologies size %lu\r\n", children.size(), topologies.size());
                    
                    std::vector<topology>::iterator it;
                    
                    for(it = children.begin(); it != children.end(); ++it) {
                        
                        for(int i=0; i<world_size; i++) {
                            
                            LOG(10, 0, 0, "sending child topology %d to %d ...\r\n", t, i);
                            
                            send_size = (int)it->comm[i].senders.size();
                            rec_size = (int)it->comm[i].receivers.size();
                            
                            LOG(10, 0, 0, "sending %d child senders to rank %d ...\r\n", send_size, i);
                            MPI_Send(&send_size, 1, MPI_INT, i, 5, MPI_COMM_WORLD);
                            MPI_Send(&it->comm[i].senders[0], send_size, MPI_INT, i, 6, MPI_COMM_WORLD);
                            
                            LOG(10, 0, 0, "sending %d child receivers to rank %d ...\r\n", rec_size, i);
                            MPI_Send(&rec_size, 1, MPI_INT, i, 7, MPI_COMM_WORLD);
                            MPI_Send(&it->comm[i].receivers[0], rec_size, MPI_INT, i, 8, MPI_COMM_WORLD);
                            
                        }
                        
                    }
                    
                }
                
            }
            
            MPI_Recv(&send_size, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            isle.senders.resize(send_size);
            
            MPI_Recv(&isle.senders[0], send_size, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            LOG(10, 0, 0, "rank %d received %lu senders: ", world_rank, isle.senders.size());
            for(int i=0; i<isle.senders.size(); i++) { LOG(10, 0, 0, "%d ", isle.senders[i]); }
            LOG(10, 0, 0, "\r\n");
            
            MPI_Recv(&rec_size, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            isle.receivers.resize(rec_size);
            
            MPI_Recv(&isle.receivers[0], rec_size, MPI_INT, 0, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            LOG(10, 0, 0, "rank %d got %lu receivers: ", world_rank, isle.receivers.size());
            for(int i=0; i<isle.receivers.size(); i++) { LOG(10, 0, 0, "%d ", isle.receivers[i]); }
            LOG(10, 0, 0, "\r\n");
            
            if(t == 0) {
                
                MPI_Recv(&send_size, 1, MPI_INT, 0, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                isle.senders.resize(send_size);
                
                MPI_Recv(&isle.senders[0], send_size, MPI_INT, 0, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                LOG(10, 0, 0, "rank %d received %lu child senders: ", world_rank, isle.senders.size());
                for(int i=0; i<isle.senders.size(); i++) { LOG(10, 0, 0, "%d ", isle.senders[i]); }
                LOG(10, 0, 0, "\r\n");
                
                MPI_Recv(&rec_size, 1, MPI_INT, 0, 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                isle.receivers.resize(rec_size);
                
                MPI_Recv(&isle.receivers[0], rec_size, MPI_INT, 0, 8, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                LOG(10, 0, 0, "rank %d got %lu child receivers: ", world_rank, isle.receivers.size());
                for(int i=0; i<isle.receivers.size(); i++) { LOG(10, 0, 0, "%d ", isle.receivers[i]); }
                LOG(10, 0, 0, "\r\n");
                    
            }
            
            eval_stats.eval_start = std::clock();
        
            isle.calc_cpd();
            
            std::vector<individual> children = crossover(isle, offsets);
            
            select_survivors(isle, children, config::mu/world_size);

            double migrate_start = MPI_Wtime();
            isle.send_migrant(tcomm);
            isle.receive_migrant(tcomm);
            double migrate_end = MPI_Wtime();
            double migrate_time = migrate_end - migrate_start;
            
            eval_stats.total_migrate_time += migrate_time;
            
            MPI_Reduce(&migrate_time, &topologies[t].fitness, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

            topologies[t].rounds++;
            
            if(topologies[t].fitness >= 0.0) {
                topologies[t].fitness = topologies[t].fitness * -1;
                topologies[t].round_fitness += (topologies[t].fitness * -1);
            } else {
                LOG(1, 0, 0, "FOUND NEGATIVE FITNESS: %2.10f in topology %d\r\n", topologies[t].fitness, t);
                topologies[t].fitness = topologies[t].fitness;
                topologies[t].round_fitness += topologies[t].fitness;
            }
            
            LOG(10, world_rank, 0, "rank %d migrate time = %2.10f sending to topology %d of %lu size %lu fitness %2.10f rounds %d\r\n", world_rank, migrate_time, t, topologies.size(), topologies[t].comm.size(), topologies[t].fitness, topologies[t].rounds);
            
            LOG(10, 0, 0, "gathering population, subpopulation %d size %d ...\r\n", world_rank, subpopulation_size);
            double gather_start = MPI_Wtime();
            MPI_Gather(&isle.population[0], subpopulation_size, individual_type, &population[0], subpopulation_size, individual_type, 0, tcomm);
            double gather_end = MPI_Wtime();
            double gather_time = gather_end - gather_start;
            LOG(10, world_rank, 0, "population gathered ...\r\n");
            
            eval_stats.total_gather_time += gather_time;
            
            if(world_rank == 0 && eval % 100 == 0) {
            
                LOG(8, 0, 0, "population size %lu, member = %2.10f\r\n", population.size(), population[population.size()-1].fitness);
                
                std::sort(population.begin(), population.end(), compare_fitness);
                std::reverse(population.begin(), population.end());
                
                std::vector<topology>::iterator min = std::min_element(topologies.begin(), topologies.end(), compare_topo_fitness);
                std::replace_if(topologies.begin(), topologies.end(), is_zero, *min);
                
                std::sort(topologies.begin(),topologies.end(), compare_topo_fitness);
                std::reverse(topologies.begin(), topologies.end());
                
                log_fn_eval_stats(population, topologies, run, eval, eval_stats, run_stats);
                
            }
            
            if (world_rank == 0) {
                
                LOG(10, 0, 0, "run %d eval %d topology %d fitness = %2.10f (round %2.10f)\r\n", run, eval, t, topologies[t].fitness, topologies[t].round_fitness);

                eval_stats.topo_migrate_time = 0.0;
                
                if(topologies[t].fitness <= eval_stats.best_topology.fitness) {
                    
                    eval_stats.best_topology = topologies[t];
                    
                    fprintf(config::topo_out, "---- RUN %d EVAL %d TOPOLOGY BEST FITNESS %3.10f ----\r\n\r\n", run, eval, topologies[t].fitness);

                    for(int i=0; i<island_ids.size(); i++) {
                      
                      for(int j=0; j<island_ids.size(); j++) {
                          
                          if(std::find(topologies[t].comm[i].receivers.begin(), topologies[t].comm[i].receivers.end(), j) != topologies[t].comm[i].receivers.end()) {
                              LOG(8, 0, 0, "1, ");
                              fprintf(config::topo_out,"1, ");
                          } else {
                              LOG(8, 0, 0, "0, ");
                              fprintf(config::topo_out,"0, ");
                          }
                          
                      }
                      
                      LOG(8, 0, 0, "\r\n");
                      fprintf(config::topo_out, "\r\n");
                      fflush(config::topo_out);
                      
                    }
                    
                }
                
                if(eval%(config::topo_mu+config::topo_lambda) == 0) {

                    LOG(10, world_rank, 0, "truncating ...\r\n");

                    std::vector<topology>::iterator min = std::min_element(topologies.begin(), topologies.end(), compare_topo_fitness);
                    std::replace_if(topologies.begin(), topologies.end(), is_zero, *min);
                    
                    std::sort(topologies.begin(),topologies.end(), compare_topo_fitness);
                    std::reverse(topologies.begin(), topologies.end());
                    
                    topologies.erase(topologies.begin()+config::topo_mu, topologies.end());

                    LOG(6, world_rank, 0, "topology 0 fitness = %2.10f\r\n", topologies[0].fitness);

                }
                
            }
                
        }
        
        if(world_rank == 0) {
            
            double run_end = MPI_Wtime();
            
            std::fprintf(config::run_stats_out, "%d,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%d,%d,%2.10f,%2.10f,%2.10f\r\n", run, eval_stats.global_best_fitness, eval_stats.average_local_best_fitness, eval_stats.average_global_best_fitness, eval_stats.total_scatter_time, eval_stats.total_gather_time, eval_stats.total_migrate_time, run_end - run_start, run_stats.init_duration, world_size, subpopulation_size, eval_stats.global_best_topo_fitness, eval_stats.average_local_best_topo_fitness, eval_stats.average_global_best_topo_fitness);
            
            LOG(2, 0, 0, "%d,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%d,%d,%2.10f,%2.10f,%2.10f\r\n", run, eval_stats.global_best_fitness, eval_stats.average_local_best_fitness, eval_stats.average_global_best_fitness, eval_stats.total_scatter_time, eval_stats.total_migrate_time, eval_stats.total_gather_time, run_end - run_start, run_stats.init_duration, world_size, subpopulation_size, eval_stats.global_best_topo_fitness, eval_stats.average_local_best_topo_fitness, eval_stats.average_global_best_topo_fitness);
            
        }
        
        fflush(config::run_stats_out);
        fflush(config::stats_out);
        
    }

    if(world_rank == 0) {
        
        for(int i=0; i<DIM; i++) {
            std::fprintf(config::solution_out, "%2.10f,", run_stats.solution.input[i]);
        }
        
        fclose(config::stats_out);
        fclose(config::run_stats_out);
        fclose(config::log_out);
        fclose(config::solution_out);
        fclose(config::topo_out);
        
    }
    
    MPI_Finalize();
        
}

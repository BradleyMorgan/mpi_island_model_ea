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

#pragma mark FUNCTION: initial_population()

// returns a vector of a randomly generated offset rastrigin solution population                        |
// accepts an array of floating point numbers of n-dimensions @offsets[param:0]                         |
// used to calculate the solution fitness.                                                              |

std::vector<individual> initial_population(std::array<double, DIM> &offsets) {
    
    // the "individual" datatype represents a single offset rastrigin solution as
    // an array of size DIM = @dim[config.txt:12] holding the solution's randomly
    // generated gene values.
    
    std::vector<individual> population;
    
    for(int i=0; i<config::mu; i++) {
        
        individual p;
        
        for (int j = 0; j < DIM; j++) {
            p.input[j] = drand(-5.12, 5.12); // rastrigin says: x[i] ∈ [-5.12,5.12]
        }
        
        p.fitness = offset_rastrigin(p.input, offsets);
        
        population.push_back(p);
        
    }
    
    // we have our initial primary population with the calculated fitnesses
    
    return population;
    
}

#pragma mark FUNCTION: initial_topo_population()

// returns a collection of randomly generated adjaceny matrices, representing                           |
// an island (communication) topology.  accepts a reference to a list of island                         |
// (process) identifiers @ids[param:0] to use as indices.                                               |

std::vector<topology> initial_topo_population(std::vector<int> &ids) {
    
    LOG(8, 0, 0, "initializing topology population for %lu islands\r\n", ids.size());
    
    std::vector<topology> population;
    
    for(int i=0; i<config::topo_mu+config::topo_lambda; i++) {
        
        topology t;
        
        LOG(10, 0, 0, "adding topology %d\r\n", i);

        t.comm = create_dyn_topology(ids);
        
        population.push_back(t);
        
    }
    
    LOG(8, 0, 0, "initialized topology population for %lu islands\r\n", ids.size());
    
    return population;
    
}

#pragma mark FUNCTION: main()

// core initialization, overarching EA logic.

// initialize globals -> experimental run [n] -> initialize populations ->                              |
// iterate topology[t] evaluation * @config.txt[min_evals] -> apply topology[t] to rastrigin[r] ->      |
// iterate rastrigin[r] evaluation * @config.txt[min_evals] -> next r                                   |

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
    
        // start of experimental run
        
        LOG(2, world_rank, 0, "**** RUN %d ****\r\n", run);
        
        double run_start = MPI_Wtime();

        estats eval_stats;
        
        LOG(4, 0, 0, "initializing run variables...\r\n");
        
        // each MPI process will maintain its own population, so we define an island for each ...
        
        LOG(4, 0, 0, "initializing rastrigin population...\r\n");
        
        std::vector<individual> population;
        
        LOG(4, 0, 0, "clearing rastrigin population...\r\n");
        
        population.clear();
        
        LOG(4, 0, 0, "resizing rastrigin population...\r\n");
        
        population.resize(config::mu);
        
        LOG(4, 0, 0, "initializing topology population...\r\n");
        
        std::vector<topology> topologies;
        topologies.clear();
        topologies.resize(config::topo_mu+config::topo_lambda);
        
        LOG(4, 0, 0, "initializing run variables...\r\n");
        
        island isle;
        isle.id = world_rank;
        isle.senders.clear();
        isle.receivers.clear();
        isle.population.clear();
        isle.population.resize(subpopulation_size);
        isle.cpd.clear();
        
        std::array<double, DIM> offsets = generate_offsets(-2.5, 2.5, .5);
        
        LOG(4, 0, 0, "creating initial populations...\r\n");
        
        // only the root process will create the full initial population ...
        
        std::vector<int> island_ids;
        
        for(int i=0; i<world_size; i++) {
         
            island_ids.push_back(i);
            
        }
        
        int send_size = 0;
        int rec_size = 0;
        
        if(world_rank == 0) {
            
            // create two populations, one for the objective function and one for the island topologies
            
            population = initial_population(offsets);
            topologies = initial_topo_population(island_ids);
            
            LOG(6, 0, 0, "received %lu isle ids for %lu topologies\r\n", island_ids.size(), topologies.size());
            
            std::sort(population.begin(), population.end(), compare_fitness);
            std::reverse(population.begin(), population.end());
            
            eval_stats.global_best_fitness = population[0].fitness;
            eval_stats.best_topology = topologies[0];
            
        }

        LOG(4, 0, 0, "duplication communicator...\r\n");
        
        MPI_Comm tcomm;
        MPI_Comm_dup(MPI_COMM_WORLD, &tcomm);
        
        // separate the single full population from the root process to subpopulations across all processes ...
        
        LOG(4, 0, 0, "scattering population ...\r\n");
        double scatter_start = MPI_Wtime();
        MPI_Scatter(&population[0], subpopulation_size, individual_type, &isle.population[0], subpopulation_size, individual_type, 0, tcomm);
        double scatter_end = MPI_Wtime();
        double scatter_time = scatter_end - scatter_start;
        LOG(4, 0, 0, "population scattered...\r\n");
        
        eval_stats.total_scatter_time += scatter_time;
            
        // begin evolution ...
        
        int rindex = 0;
        
        LOG(4, 0, 0, "starting evals...\r\n");
        
        // TODO: finish new loop logic stub
        
        for(int topo_idx = 0; topo_idx <= config::topo_mu; topo_idx++) {
            
            // distribute topologies[topo_idx]
            
            for(int rast_idx = 0; rast_idx <= config::mu; rast_idx++) {
                
                //
                
            }
            
        }
        
        // TODO: Remove the OLD multi-objective logic
        // needs significant changes, borrow from this for the new logic stub (above)
        
        for(int eval=1; eval<=config::evals; eval++) {  // eval start
                    
            // all topologies must have been evaluated and assigned a fitness before we can perform topology evolution cycle operations,
            // so track that interval with var t, and also use it as an index increment for the topology population so that topology[t]
            // will be employed by the islands for the migration pattern within the corresponding evaluation ...
            
            int t = (eval-1)%(config::topo_mu);
            int c = (eval-1)%(config::topo_mu+config::topo_lambda);
            int m = (eval-1)%((config::topo_mu+config::topo_lambda)*config::topo_evals);
            int tindex = (eval-1)%(config::topo_mu+config::topo_lambda);
            
            // determine the topology individual to evaluate from the current eval number and distribute it if needed ...
            
            #pragma mark LOGIC START: topology distribution
            
            // rank 0 will determine the topology to be distributed, if needed and initiate send operations  ...
            
            if(world_rank == 0) {  // start topology send logic
                
                // rotate back to the first individual if we have reached the min required evals ...
                
                if(rindex > 0 && eval%(rindex*config::topo_evals) == 0) {
                    rindex = 0;
                }
                
                // find a topology that has not been evaluated min evals ...
                
                if(rindex > (config::topo_evals - 1)) {
                    rindex = 0;
                    while(topologies[rindex].rounds >= config::topo_evals) {
                        rindex++;
                    }
                }
                
                // if for some reason we didn't find any topologies needed to evaluate, start over ...
                
                if(rindex > (config::topo_evals - 1)) {
                    rindex = 0;
                }
                
                LOG(10, world_rank, 0, "C = (%d)mod(%d) = %d | T = (%d)mod(%d) = %d\r\n", eval, (config::topo_mu+config::topo_lambda), c, eval-1, config::topo_mu, t);
                
                // create a new topology generation ...
                
                if(eval != 1 && m == 0) {
                    
                    LOG(4, world_rank, 0, "adding children...\r\n");
                    std::vector<topology> children = topo_gen(topologies, world_size);
                    topologies.insert(topologies.end(), children.begin(), children.end());
                    
                    LOG(4, world_rank, 0, "added %lu children to topologies, new topologies size %lu\r\n", children.size(), topologies.size());
                    
                }
                
                // if we have evaluated this topology min times, increment the topology index and distribute it ...
                
                if(topologies[rindex].rounds >= config::topo_evals || eval == 1) {
    
                    rindex++;
                    
                    for(int i=0; i<world_size; i++) {
                        
                        // set the senders and receivers for each island, so that we use topology[rindex] for the next rastrigin evaluations ...
                        
                        LOG(10, world_rank, 0, "TOPOLOGY INDEX %d | C = %d | T = (%d)mod(%d) = %d\r\n", rindex, c, eval-1, config::topo_mu, t);
                        
                        send_size = (int)topologies[rindex].comm[i].senders.size();
                        rec_size = (int)topologies[rindex].comm[i].receivers.size();
                        
                        LOG(4, 0, 0, "sending topology %d (senders = %d, receivers = %d) to %d ...\r\n", rindex, send_size, rec_size, i);
                        
                        LOG(4, 0, 0, "sending topology %d, %d senders to rank %d ...\r\n", rindex, send_size, i);
                        MPI_Send(&send_size, 1, MPI_INT, i, 1, tcomm);
                        MPI_Send(&topologies[rindex].comm[i].senders[0], send_size, MPI_INT, i, 2, tcomm);
                        
                        LOG(4, 0, 0, "sending topology %d, %d receivers to rank %d ...\r\n", rindex, rec_size, i);
                        MPI_Send(&rec_size, 1, MPI_INT, i, 3, tcomm);
                        MPI_Send(&topologies[rindex].comm[i].receivers[0], rec_size, MPI_INT, i, 4, tcomm);
                        
                    }
                    
                }
                
            }  // end topology send logic
            
            // rank 0, distribute the topology metadata ...
            
            MPI_Bcast(&rindex, 1, MPI_INT, 0, tcomm);
            MPI_Bcast(&topologies[rindex].rounds, 1, MPI_INT, 0, tcomm);
            
            if(topologies[rindex].rounds >= config::topo_evals || eval == 1) {  // begin topology receive logic
            
                // previous topology was evaluated min times, so the islands will receive a new one ...
                
                LOG(4, 0, 0, "rank %d receiving topology %d\r\n", world_rank, rindex);
                
                MPI_Recv(&send_size, 1, MPI_INT, 0, 1, tcomm, MPI_STATUS_IGNORE);
                isle.senders.clear();
                isle.senders.resize(send_size);
                
                MPI_Recv(&isle.senders[0], send_size, MPI_INT, 0, 2, tcomm, MPI_STATUS_IGNORE);
                LOG(4, 0, 0, "rank %d received %lu senders: ", world_rank, isle.senders.size());
                for(int i=0; i<isle.senders.size(); i++) { LOG(10, 0, 0, "%d ", isle.senders[i]); }
                LOG(4, 0, 0, "\r\n");
                
                MPI_Recv(&rec_size, 1, MPI_INT, 0, 3, tcomm, MPI_STATUS_IGNORE);
                isle.receivers.clear();
                isle.receivers.resize(rec_size);
                
                MPI_Recv(&isle.receivers[0], rec_size, MPI_INT, 0, 4, tcomm, MPI_STATUS_IGNORE);
                LOG(10, 0, 0, "rank %d got %lu receivers: ", world_rank, isle.receivers.size());
                for(int i=0; i<isle.receivers.size(); i++) { LOG(10, 0, 0, "%d ", isle.receivers[i]); }
                LOG(10, 0, 0, "\r\n");
            
            }  // end topology receive logic
            
            
            #pragma mark LOGIC START: primary (rastrigin) population evolution
            
            
            eval_stats.eval_start = std::clock();
        
            isle.calc_cpd();
            
            // parent selection and recombination ...
            
            std::vector<individual> children = crossover(isle, offsets);
            
            // selection ...
            
            select_survivors(isle, children, config::mu/world_size);

            
            #pragma mark CORE START: island migrations
            

            double migrate_start = MPI_Wtime();
            
            isle.send_migrant(tcomm);
            isle.receive_migrant(tcomm);
            
            double migrate_end = MPI_Wtime();
            double migrate_time = migrate_end - migrate_start;
            
            LOG(8, 0, 0, "migrate start = %3.10f, migrate end = %3.10f, migrate time = %3.10f\r\n", migrate_start, migrate_end, migrate_time);
            
            eval_stats.total_migrate_time += migrate_time;
            
            // aggregate the migration time from each island to get the total topology fitness ...
            
            MPI_Reduce(&migrate_time, &topologies[rindex].fitness, 1, MPI_DOUBLE, MPI_SUM, 0, tcomm);

            LOG(5, world_rank, 0, "total migration time -> %2.10f\r\n", topologies[rindex].fitness);
            
            if(world_rank == 0) {
                topologies[rindex].rounds++;
            }
                
            LOG(4, 0, 0, "rank %d rindex = %d rounds = %d eval = %d\r\n", world_rank, rindex, topologies[rindex].rounds, eval);
            
            // workaround if for some reason we received an anomalous topology fitness value
            // this shouldn't happen so log it for review ...
            
            if(topologies[rindex].fitness >= 0.0) {
                topologies[rindex].fitness = topologies[rindex].fitness * -1;
                topologies[rindex].round_fitness += topologies[rindex].fitness;
            } else {
                LOG(6, 0, 0, "FOUND NEGATIVE FITNESS: %2.10f in topology %d\r\n", topologies[rindex].fitness, tindex);
            }
            
            LOG(10, world_rank, 0, "rank %d migrate time = %2.10f sending to topology %d of %lu size %lu fitness %2.10f rounds %d\r\n", world_rank, migrate_time, tindex, topologies.size(), topologies[rindex].comm.size(), topologies[rindex].fitness, topologies[rindex].rounds);
            
            
            #pragma mark CORE END: island migrations
            
            
            // aggregate primary (rastrigin) population results from islands ...
            
            LOG(10, 0, 0, "gathering population, subpopulation %d size %d ...\r\n", world_rank, subpopulation_size);
            
            double gather_start = MPI_Wtime();
            
            MPI_Gather(&isle.population[0], subpopulation_size, individual_type, &population[0], subpopulation_size, individual_type, 0, tcomm);
            
            double gather_end = MPI_Wtime();
            double gather_time = gather_end - gather_start;
            
            LOG(10, world_rank, 0, "population gathered ...\r\n");
            
            eval_stats.total_gather_time += gather_time;
            
            if(eval%100 == 0) {
            
                LOG(8, 0, 0, "population size %lu, member = %2.10f\r\n", population.size(), population[population.size()-1].fitness);
                
                std::sort(population.begin(), population.end(), compare_fitness);
                std::reverse(population.begin(), population.end());
                
                LOG(10, 0, 0, "topologies size %lu, member = %2.10f\r\n", topologies.size(), topologies[topologies.size()-1].fitness);
                
                std::vector<topology>::iterator min = std::min_element(topologies.begin(), topologies.end(), compare_topo_fitness);
                std::replace_if(topologies.begin(), topologies.end(), is_zero, *min);
                
                std::sort(topologies.begin(),topologies.end(), compare_topo_fitness);
                std::reverse(topologies.begin(), topologies.end());
            
            }
            
            // calculate and log the stats from the current population states ...
            
            if(world_rank == 0 && eval % 100 == 0) {
                
                log_fn_eval_stats(population, topologies, run, eval, eval_stats, run_stats);
                
            }
            
            if (world_rank == 0) {
                
                // log the best topology's adjacency matrix in csv format ...
                
                LOG(10, 0, 0, "run %d eval %d topology %d fitness = %2.10f (round %2.10f)\r\n", run, eval, tindex, topologies[tindex].fitness, topologies[tindex].round_fitness);

                eval_stats.topo_migrate_time = 0.0;  // ?? what
                
                if(topologies[tindex].fitness >= eval_stats.best_topology.fitness) {
                    
                    eval_stats.best_topology = topologies[tindex];
                    
                    fprintf(config::topo_out, "---- RUN %d EVAL %d TOPOLOGY BEST FITNESS %3.10f ----\r\n\r\n", run, eval, topologies[tindex].fitness);

                    for(int i=0; i<island_ids.size(); i++) {
                      
                      for(int j=0; j<island_ids.size(); j++) {
                          
                          if(std::find(topologies[tindex].comm[i].receivers.begin(), topologies[tindex].comm[i].receivers.end(), j) != topologies[tindex].comm[i].receivers.end()) {
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
                
                // perform survival selection for topology population ...
                
                if(eval != 1 && c == 0) {

                    LOG(4, world_rank, 0, "truncating topologies size %lu at %d...\r\n", topologies.size(), eval);

                    std::vector<topology>::iterator min = std::min_element(topologies.begin(), topologies.end(), compare_topo_fitness);
                    std::replace_if(topologies.begin(), topologies.end(), is_zero, *min);

                    LOG(6, world_rank, 0, "replacing 0 fit with %2.10f\r\n", min->fitness);

                    std::sort(topologies.begin(),topologies.end(), compare_topo_fitness);
                    std::reverse(topologies.begin(), topologies.end());

                    topologies.erase(topologies.begin()+config::topo_mu, topologies.end());

                    LOG(4, world_rank, 0, "removed %d topologies, new size = %lu\r\n", config::topo_lambda, topologies.size());
                    LOG(10, world_rank, 0, "topology 0 fitness = %2.10f\r\n", topologies[0].fitness);
                    
                    if(t == 0) {
                        tindex = 0;
                    } else {
                        tindex++;
                    }

                }
                
            }
                
        } // eval end
        
        // log aggregated run stats ...
        
        if(world_rank == 0) {
            
            double run_end = MPI_Wtime();
            
            std::fprintf(config::run_stats_out, "%d,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%d,%d,%2.10f,%2.10f,%2.10f\r\n", run, eval_stats.global_best_fitness, eval_stats.average_local_best_fitness, eval_stats.average_global_best_fitness, eval_stats.total_scatter_time, eval_stats.total_gather_time, eval_stats.total_migrate_time, run_end - run_start, run_stats.init_duration, world_size, subpopulation_size, eval_stats.global_best_topo_fitness, eval_stats.average_local_best_topo_fitness, eval_stats.average_global_best_topo_fitness);
            
            LOG(2, 0, 0, "%d,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%d,%d,%2.10f,%2.10f,%2.10f\r\n", run, eval_stats.global_best_fitness, eval_stats.average_local_best_fitness, eval_stats.average_global_best_fitness, eval_stats.total_scatter_time, eval_stats.total_migrate_time, eval_stats.total_gather_time, run_end - run_start, run_stats.init_duration, world_size, subpopulation_size, eval_stats.global_best_topo_fitness, eval_stats.average_local_best_topo_fitness, eval_stats.average_global_best_topo_fitness);
            
            rindex = 0;
            
        }
        
        fflush(config::run_stats_out);
        fflush(config::stats_out);
        
    }  // run end
  
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

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

#define RED   "\x1B[31m"
#define GRN   "\x1B[32m"
#define YEL   "\x1B[33m"
#define BLU   "\x1B[34m"
#define MAG   "\x1B[35m"
#define CYN   "\x1B[36m"
#define WHT   "\x1B[37m"
#define RESET "\x1B[0m"

struct rstats {
    
    solution sol;
    
    double init_duration;
    
};

struct estats {
    
    int topo_best_count = 0;
    
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
    
    solution sol;
    
    void init() {
        
        this->init_duration = 0.0;
        this->eval_start = 0.0;
        this->eval_duration = 0.0;
        this->total_scatter_time = 0.0;
        this->total_gather_time = 0.0;
        this->total_migrate_time = 0.0;
        this->topo_migrate_time = 0.0;
        this->local_best_fitness = 0.0;
        this->average_local_best_fitness = 0.0;
        this->local_best_topo_fitness = 0.0;
        this->average_local_best_topo_fitness = 0.0;
        this->average_local_best_fitnesses.clear();
        this->average_local_best_fitnesses.resize(0);
        this->average_local_best_topo_fitnesses.clear();
        this->average_local_best_topo_fitnesses.resize(0);
        
    }
    
};

void log_topology_matrix(int world_size, topology &t, int count) {
    
    sprintf(config::topo_fname, "%s/topo_%d_%d_%ld.py", config::topos_subpath, count, world_size, time(0));
    config::topo_out = fopen(config::topo_fname, "w");
    
    fprintf(config::topo_out,"matrix = [");
    
    for(int i=0; i<world_size; i++) {
    
        fprintf(config::topo_out,"[");
      
        for(int j=0; j<world_size; j++) {
            
            if(std::find(t.channels[i].receivers.begin(), t.channels[i].receivers.end(), j) != t.channels[i].receivers.end()) {
                LOG(8, 0, 0, "1,");
                fprintf(config::topo_out,"1");
                
            } else {
                LOG(8, 0, 0, "0, ");
                fprintf(config::topo_out,"0");
                
            }
            
            if(j <= world_size-2) {
                fprintf(config::topo_out,",");
            } else {
                fprintf(config::topo_out,"]");
            }
            
        }
      
        if(i <= world_size-2) {
            fprintf(config::topo_out,",\r\n          ");
        } else {
            fprintf(config::topo_out,"]\r\n");
        }
        
    }
    
    fflush(config::topo_out);
    
}

void log_fn_eval_stats(std::vector<solution> &population, std::vector<topology> &topologies, int &run, int &eval, estats &eval_stats, rstats &run_stats, topology &t) {
    
    // only consider evaluated topologies for stats
    
    std::vector<topology> filtered;
    
    std::copy_if( topologies.begin(), topologies.end(), std::back_inserter(filtered), [](const topology &item) { return item.fitness != 0.0 && item.rounds >= config::topo_evals; });
    std::sort(filtered.begin(), filtered.end(), compare_topo_fitness);
    std::reverse(filtered.begin(), filtered.end());
    
    LOG(3, 0, 0, "STATS (run %d, eval %d): solutions size %lu, mem[0] fit = %f \r\n", run, eval, population.size(), population[0].fitness);
    
    eval_stats.local_best_fitness = population[0].fitness;
    eval_stats.average_local_best_fitnesses.push_back(population[0].fitness);
    
    eval_stats.local_best_topo_fitness = filtered[0].fitness;
    eval_stats.average_local_best_topo_fitnesses.push_back(filtered[0].fitness);
    
    LOG(3, 0, 0, "STATS (run %d, eval %d): local best = %f, locat best fitnesses tracked = %lu\r\n", run, eval, eval_stats.local_best_fitness, eval_stats.average_local_best_fitnesses.size());
    
    if(population[0].fitness > eval_stats.global_best_fitness) {
        LOG(3, 0, 0, "STATS (run %d, eval %d): found new global best solution %f > %f\r\n", run, eval, population[0].fitness, eval_stats.global_best_fitness);
        eval_stats.sol = population[0];
        eval_stats.average_global_best_fitnesses.push_back(population[0].fitness);
        eval_stats.global_best_fitness = population[0].fitness;
    }
    
    if(filtered[0].fitness > 0) {
        LOG(1, 0, 0, "FOUND NEGATIVE FITNESS: %2.10f in topology\r\n", topologies[0].fitness);
    }
    
    if(filtered[0].fitness > eval_stats.global_best_topo_fitness || eval_stats.global_best_topo_fitness == 0.0) {
        LOG(3, 0, 0, "STATS (run %d, eval %d): found new global best topology %f > %f\r\n", run, eval, topologies[0].fitness, eval_stats.global_best_topo_fitness);
        eval_stats.topo_best_count++;
        eval_stats.best_topology = filtered[0];
        eval_stats.average_global_best_topo_fitnesses.push_back(filtered[0].fitness);
        eval_stats.global_best_topo_fitness = filtered[0].fitness;
        log_topology_matrix(6, filtered[0], eval_stats.topo_best_count);
    }
    
    double total_fitness = 0.0;
    double total_topo_fitness = 0.0;
    
    for(int i=0; i<population.size(); i++) {
        total_fitness += population[i].fitness;
    }
    
    for(int i=0; i<filtered.size(); i++) {
        total_topo_fitness += filtered[i].fitness;
    }
    
    LOG(3, 0, 0, "STATS (run %d, eval %d): solution population size %lu total fitness = %2.10f\r\n", run, eval, population.size(), total_fitness);
    
    if(eval_stats.average_local_best_fitnesses.size() > 1) {
        eval_stats.average_local_best_fitness = std::accumulate(eval_stats.average_local_best_fitnesses.begin(), eval_stats.average_local_best_fitnesses.end(), 0.0) / eval_stats.average_local_best_fitnesses.size();
    } else {
        eval_stats.average_local_best_fitness = population[0].fitness;
    }
    
    eval_stats.average_global_best_fitness = std::accumulate(eval_stats.average_global_best_fitnesses.begin(), eval_stats.average_global_best_fitnesses.end(), 0.0) / eval_stats.average_global_best_fitnesses.size();
   
    eval_stats.average_local_best_topo_fitness = std::accumulate(eval_stats.average_local_best_topo_fitnesses.begin(), eval_stats.average_local_best_topo_fitnesses.end(), 0.0) /
        eval_stats.average_local_best_topo_fitnesses.size();
    eval_stats.average_global_best_topo_fitness = std::accumulate(eval_stats.average_global_best_topo_fitnesses.begin(), eval_stats.average_global_best_topo_fitnesses.end(), 0.0) / eval_stats.average_global_best_topo_fitnesses.size();
    
    double average_fitness = total_fitness / population.size();
    double average_gather_time = eval_stats.total_gather_time / eval;
    double average_scatter_time = eval_stats.total_scatter_time / eval;
    double average_migrate_time = eval_stats.total_migrate_time / eval;
    double average_topo_fitness = filtered.size() > 0 ? (total_topo_fitness / filtered.size()) : 0.0;
    
    double eval_duration = ( std::clock() - eval_stats.eval_start ) / (double) CLOCKS_PER_SEC;
    
    //std::fprintf(config::stats_out, "%d,%d,%4.10f,%4.10f,%4.10f,%4.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%d,%2.10f,%d,%d,%d,%d,%d,%2.10f\r\n", run, eval, average_fitness, eval_stats.local_best_fitness, eval_stats.global_best_fitness, eval_stats.average_local_best_fitness, eval_stats.average_global_best_fitness, average_scatter_time, average_gather_time, average_migrate_time, eval_stats.init_duration, eval_duration, average_topo_fitness, eval_stats.local_best_topo_fitness, eval_stats.global_best_topo_fitness, eval_stats.average_local_best_topo_fitness, eval_stats.average_global_best_topo_fitness, filtered[0].id, filtered[0].round_fitness, filtered[0].rounds, filtered[0].channel_count, t.id, t.rounds, t.channel_count, t.fitness);
    
    
    std::fprintf(config::stats_out, "%d," "%d,"  "%3.10f,"         "%3.10f,"                      "%3.10f,"                       "%13.10f,"                             "%3.10f,",
                                    run,  eval,  average_fitness,  eval_stats.local_best_fitness, eval_stats.global_best_fitness, eval_stats.average_local_best_fitness, eval_stats.average_global_best_fitness);

    
    std::fprintf(config::stats_out, "%3.10f,"             "%3.10f,"            "%3.10f,"             "%3.10f,"               "%3.10f,",
                                    average_scatter_time, average_gather_time, average_migrate_time, run_stats.init_duration, eval_duration);
    
    std::fprintf(config::stats_out, "%3.10f,"             "%d,"            "%d,"               "%d,"                       "%3.10f,"                   "%3.10f,"            "%3.10f,"                            "%3.10f,"                           "%3.10f," ,
                                    average_topo_fitness, filtered[0].id,  filtered[0].rounds, filtered[0].channel_count,  filtered[0].round_fitness,  filtered[0].fitness, eval_stats.local_best_topo_fitness,  eval_stats.global_best_topo_fitness, eval_stats.average_local_best_topo_fitness);
    
    std::fprintf(config::stats_out, "%3.10f,"                                    "%d,"       "%d,"       "%d,"              "%3.10f\r\n",
                                    eval_stats.average_global_best_topo_fitness, t.id,       t.rounds,   t.channel_count,    t.fitness);
    
    
    LOG(3, 0, 0, "%-2s,%-4s,%-13s,%-13s,%-13s,%-13s,%-13s,%-12s,%-12s,%-12s,%-12s,%-12s,%-13s,%-13s,%-13s,%-13s,%-13s,%-13s,%-13s\r\n", "r", "e", "avg_fitness", "lb_fitness", "gb_fitness", "alb_fitness", "agb_fitness", "avg_scatter", "avg_gather", "avg_migrate", "init_time", "eval_time", "avgt_fit", "lbt_fit", "gbt_fit", "albt_fit", "agbt_fit", "round_fit", "rounds");
    
    
    
    LOG(2, 0, 0, "%3d," "%5d,"  "%013.10f,"         "%3.10f,"                      "%3.10f,"                       "%13.10f,"                              "%3.10f,"                               "%3.10f,"             "%3.10f,"            "%3.10f,"             "%3.10f,"               "%013.10f"     "%s",
                 run,   eval,   average_fitness,    eval_stats.local_best_fitness, eval_stats.global_best_fitness, eval_stats.average_local_best_fitness, eval_stats.average_global_best_fitness, average_scatter_time, average_gather_time, average_migrate_time, run_stats.init_duration, eval_duration, "|");
        
    LOG(2, 0, 0, "%3.10f,"             "%03d,"          "%05d,"             "%03d,"                     "%3.10f,"                   "%3.10f,"            "%3.10f,"                            "%s"   "%3.10f,"                           "%s"  "%3.10f"                                    "%s,",
                 average_topo_fitness, filtered[0].id,  filtered[0].rounds, filtered[0].channel_count,  filtered[0].round_fitness,  filtered[0].fitness, eval_stats.local_best_topo_fitness,  BLU,   eval_stats.global_best_topo_fitness, YEL, eval_stats.average_local_best_topo_fitness, RESET);
    
    LOG(2, 0, 0, "%3.10f,"                                    "%03d,"     "%05d,"     "%03d,"             "%013.10f,"    "%013.10f,"          "%013.10f\r\n",
                 eval_stats.average_global_best_topo_fitness, t.id,       t.rounds,   t.channel_count,    t.fitness,        t.total_migration_time, eval_stats.total_migrate_time);
    
    
    
}

#endif /* stats_h */

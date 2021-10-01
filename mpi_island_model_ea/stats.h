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
    
    topology *best_topology;
    
    std::vector<double> average_local_best_fitnesses;
    std::vector<double> average_global_best_fitnesses;
    std::vector<double> average_local_best_topo_fitnesses;
    std::vector<double> average_global_best_topo_fitnesses;
    
    solution *sol;
    
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
    
    sprintf(config::topo_fname, "%s/topo_%03d_%d_%ld.py", config::topos_subpath, count, world_size, time(0));
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
            fprintf(config::topo_out,",\r\n");
        } else {
            fprintf(config::topo_out,"]\r\n");
        }
        
    }
    
    //fflush(config::topo_out);
    fclose(config::topo_out);
    
}

void log_pop_stats(int &run, int &eval, std::vector<solution> &solutions, island &isle, MPI_Datatype &visa_type) {
    
    std::vector<visa> visas;
    unsigned long int vrecsize = 0;
    unsigned long int vsndsize = isle.visas.size();

    if(isle.id==0) {
        for(int i=1; i<config::world_size; i++) {
            MPI_Recv(&vrecsize, 1, MPI_INT, i, i, isle.tcomm, MPI_STATUS_IGNORE);
            std::vector<visa> v;
            v.resize(vrecsize);
            MPI_Recv(&v[0], int(vrecsize), visa_type, i, i*10, isle.tcomm, MPI_STATUS_IGNORE);
            visas.insert(visas.end(), v.begin(), v.end());
        }
        visas.insert(visas.end(), isle.visas.begin(), isle.visas.end());
    } else {
        MPI_Send(&vsndsize, 1, MPI_INT, 0, isle.id, isle.tcomm);
        MPI_Send(&isle.visas[0], int(vsndsize), visa_type, 0, isle.id*10, isle.tcomm);        
    }
    
    if(isle.id==0) {
            
        std::map<long long int, int> groups_e5;
        std::map<long long int, int> groups_e4;
        std::map<long int, int> groups_e3;
        
        for (auto it = solutions.begin(); it != solutions.end(); ++it) {
            
            long long int e5 = (it->group * 1000000000);
            ++groups_e5[e5];
            
            long long int e4 = (it->group * 1000000);
            ++groups_e4[e4];
            
            long int e3 = (it->group * 1000);
            ++groups_e3[e3];
            
            //std::sort(sol_stats.begin(), sol_stats.end());
            
            char sol[DIM*sizeof(double)];
            
            int offs = 0;
            
            for(int i=0; i<DIM; i++) {
                char delim = (i == DIM-1) ? '\0' : ';';
                offs += snprintf(sol+offs, sizeof(sol)>offs?sizeof(sol)-offs:0, "%f%c", it->input[i], delim);
            }
    
    //            char vis[it->data.visas.size()];
    //
    //            int offv = 0;
    //
    //            for(int i=0; i<it->data.visas.size(); i++) {
    //            //for(int island : it->visas) {
    //                char delim = (i == it->data.visas.size()-1) ? '\0' : ';';
    //                //char prefix = (i==0) ? '"' : ' ';
    //                //offv += snprintf(vis+offv, sizeof(vis)-offv, "%c%d%c", prefix, it->visas[i], delim);
    //                //offv += snprintf(vis+offv, sizeof(sol)>offv?sizeof(sol)-offv:0, "%c%d%c", prefix, it->visas[i], delim);
    //                offv += snprintf(vis+offv, sizeof(vis)>offv?sizeof(vis)-offv:0, "%d%c", it->data.visas[i], delim);
    //            }
                
                std::fprintf(config::solpop_out, "%d,"  "%d,"    "%s,"      "%d,"       "%d,"       "%s,"            "%s,"             "%d,"          "%d,",
                                                  run,   eval,    it->id,   it->source, it->locale, it->parents[0],   it->parents[1],  it->selected,  it->survival);
                             
                             
                             
                std::fprintf(config::solpop_out, "%lld,"  "%d,"          "%lld,"   "%d,"           "%ld,"   "%d,",
                                                 e5,      groups_e5[e5],  e4,      groups_e4[e4],   e3,     groups_e3[e3]);
                             
                             
                std::fprintf(config::solpop_out, "%f,"        "%f,"                      "%d,"             "%s,"            "%s\r\n",
                                                 it->fitness, it->selection_distribution, it->migrations, "placeholder",    sol);
            
        }
        
    }
    
}

void log_fn_topology_stats(std::vector<topology> &topologies, int &run, int &eval, estats &eval_stats, rstats &run_stats, topology &t) {
    
    LOG(3, 0, 0, "STATS (run %d, eval %d): topologies size %lu, mem[0] fit = %f \r\n", run, eval, topologies.size(), topologies[0].fitness);
    
    std::vector<topology> filtered_top;
    
    std::copy_if( topologies.begin(), topologies.end(), std::back_inserter(filtered_top), [](const topology &item) { return item.fitness != 0.00000000000 && item.rounds >= config::objective_1_max_fit_evals; });
    std::sort(filtered_top.begin(), filtered_top.end(), compare_topo_fitness);
    std::reverse(filtered_top.begin(), filtered_top.end());
    
    if(filtered_top.size() == 0) {
        filtered_top.push_back(topologies[0]);
    }
    
    LOG(3, 0, 0, "STATS (run %d, eval %d): filtered topologies size %lu, mem[0] fit = %f \r\n", run, eval, filtered_top.size(), filtered_top[0].fitness);
    
    eval_stats.local_best_topo_fitness = filtered_top[0].fitness;
    eval_stats.average_local_best_topo_fitnesses.push_back(filtered_top[0].fitness);

    LOG(3, 0, 0, "STATS (run %d, eval %d): best topology fit = %f \r\n", run, eval, eval_stats.local_best_topo_fitness);
    
    if(filtered_top[0].fitness > eval_stats.global_best_topo_fitness || eval_stats.global_best_topo_fitness == 0.0) {
        LOG(3, 0, 0, "STATS (run %d, eval %d): found new global best topology %f > %f\r\n", run, eval, topologies[0].fitness, eval_stats.global_best_topo_fitness);
        eval_stats.topo_best_count++;
        eval_stats.average_global_best_topo_fitnesses.push_back(filtered_top[0].fitness);
        eval_stats.global_best_topo_fitness = filtered_top[0].fitness;
        log_topology_matrix(t.world_size, filtered_top[0], eval_stats.topo_best_count);
    }
    
    double total_topo_fitness = 0.0;
    
    for(int i=0; i<filtered_top.size(); i++) {
        total_topo_fitness += filtered_top[i].fitness;
    }
    
    LOG(3, 0, 0, "STATS (run %d, eval %d): total topology fit = %f \r\n", run, eval, total_topo_fitness);
    
    eval_stats.average_local_best_topo_fitness = std::accumulate(eval_stats.average_local_best_topo_fitnesses.begin(), eval_stats.average_local_best_topo_fitnesses.end(), 0.0) /
        eval_stats.average_local_best_topo_fitnesses.size();
    eval_stats.average_global_best_topo_fitness = std::accumulate(eval_stats.average_global_best_topo_fitnesses.begin(), eval_stats.average_global_best_topo_fitnesses.end(), 0.0) / eval_stats.average_global_best_topo_fitnesses.size();
    
    double average_topo_fitness = filtered_top.size() > 0 ? (total_topo_fitness / filtered_top.size()) : 0.0;
    
    std::fprintf(config::topo_stats_out, "%3.10f,"             "%d,"            "%d,"               "%d,"                       "%3.10f,"                   "%3.10f,"            "%3.10f,"                            "%3.10f,"                           "%3.10f," ,
                 average_topo_fitness, filtered_top[0].id,  filtered_top[0].rounds, filtered_top[0].channel_count,  filtered_top[0].round_fitness,  filtered_top[0].fitness, eval_stats.local_best_topo_fitness,  eval_stats.global_best_topo_fitness, eval_stats.average_local_best_topo_fitness);

    std::fprintf(config::topo_stats_out, "%3.10f,"                                    "%d,"       "%d,"       "%d,"              "%3.10f\r\n",
                                    eval_stats.average_global_best_topo_fitness, t.id,       t.rounds,   t.channel_count,    t.fitness);
    

    LOG(2, 0, 0, "\r\n\r\n------ TOPOLOGY %d FITNESS=%013.10f EVALUATED %d TIMES AT %d ------", t.id, t.fitness, t.rounds, eval);
    LOG(2, 0, 0, "\r\n\r\n%5s %9s %10s %13s %13s %13s", "t_id", "t_evals", "t_channels", "t_fitness", "t_mig_time", "sum_mig_time");
    
    LOG(2, 0, 0, "\r\n%05d,"     "%09d,"     "%010d,"             "%013.10f,"    "%013.10f,"            "%013.10f,",
                 t.id,       t.rounds,   t.channel_count,    t.fitness,     t.total_migration_time, eval_stats.total_migrate_time);
    
    LOG(2, 0, 0, "\r\n\r\n%13s %12s %13s %13s %13s %13s %13s %13s %13s %13s", "avg_t_fit", "gbest_t_id", "gbest_t_rounds", "gbest_t_chan", "gbest_t_efit", "gbest_tfit1", "lbest_tfit", "gbest_tfit2", "avg_lbest_tfit", "avg_gbest_tfit");
    
    LOG(2, 0, 0, "\r\n%3.10f,"          "%012d,"               "%014d,"                 "%013d,"                         "%3.10f,"                       "%3.10f,"               "%3.10f,",
                 average_topo_fitness,   filtered_top[0].id,  filtered_top[0].rounds, filtered_top[0].channel_count,  filtered_top[0].round_fitness,  filtered_top[0].fitness, eval_stats.local_best_topo_fitness);
    
    LOG(2, 0, 0, "%s"   "%3.10f,"                           "%s,"  "%3.10f,"                                    "%3.10f,"                                   "%s",
                 BLU,   eval_stats.global_best_topo_fitness, YEL, eval_stats.average_local_best_topo_fitness, eval_stats.average_global_best_topo_fitness, RESET);
    
    LOG(2, 0, 0, "\r\n--------------------------------------------------------------------------");
    
    LOG(2, 0, 0, "\r\n\r\n%3s %5s %14s %14s %14s %13s %13s %12s %12s %12s %12s %12s", "run,","eval,","avg_fit,","lbest_fit,","gbest_fit,","avg_lbest_fit,","avg_gbest_fit,","avg_scatt_t,","avg_gath_t,","avg_migrt,","init_t,","eval_t");
    
    //LOG(2, 0, 0, "%3.10f,"                                    "%03d,"     "%05d,"     "%03d,"             "%013.10f,"    "%013.10f,"            "%013.10f\r\n",
    //             eval_stats.average_global_best_topo_fitness, t.id,       t.rounds,   t.channel_count,    t.fitness,     t.total_migration_time, eval_stats.total_migrate_time);

}

void log_fn_eval_stats(std::vector<solution> &solutions, std::vector<topology> &topologies, int &run, int &eval, estats &eval_stats, rstats &run_stats, topology &t) {
    
    // only consider evaluated topologies for stats
    
    std::vector<solution> filtered_sol;
    
    std::copy_if( solutions.begin(), solutions.end(), std::back_inserter(filtered_sol), [](const solution &item) { return item.fitness != 0.0; });
    std::sort(filtered_sol.begin(), filtered_sol.begin(), compare_fitness<solution>);
    std::reverse(filtered_sol.begin(), filtered_sol.end());
    
    LOG(3, 0, 0, "STATS (run %d, eval %d): solutions size %lu, mem[0] fit = %f \r\n", run, eval, solutions.size(), solutions[0].fitness);
    
    eval_stats.local_best_fitness = filtered_sol[0].fitness;
    eval_stats.average_local_best_fitnesses.push_back(filtered_sol[0].fitness);
        
    LOG(3, 0, 0, "STATS (run %d, eval %d): local best = %f, local best fitnesses tracked = %lu\r\n", run, eval, eval_stats.local_best_fitness, eval_stats.average_local_best_fitnesses.size());
    
    if(filtered_sol[0].fitness > eval_stats.global_best_fitness) {
        LOG(3, 0, 0, "STATS (run %d, eval %d): found new global best solution %f > %f\r\n", run, eval, filtered_sol[0].fitness, eval_stats.global_best_fitness);
        eval_stats.sol = &filtered_sol[0];
        eval_stats.average_global_best_fitnesses.push_back(filtered_sol[0].fitness);
        eval_stats.global_best_fitness = filtered_sol[0].fitness;
    }
    
    double total_fitness = 0.0;
    
    for(int i=0; i<filtered_sol.size(); i++) {
        total_fitness += filtered_sol[i].fitness;
    }
    
    LOG(3, 0, 0, "STATS (run %d, eval %d): solution population size %lu total fitness = %2.10f\r\n", run, eval, filtered_sol.size(), total_fitness);
    
    if(eval_stats.average_local_best_fitnesses.size() > 1) {
        eval_stats.average_local_best_fitness = std::accumulate(eval_stats.average_local_best_fitnesses.begin(), eval_stats.average_local_best_fitnesses.end(), 0.0) / eval_stats.average_local_best_fitnesses.size();
    } else {
        eval_stats.average_local_best_fitness = filtered_sol[0].fitness;
    }
    
    eval_stats.average_global_best_fitness = std::accumulate(eval_stats.average_global_best_fitnesses.begin(), eval_stats.average_global_best_fitnesses.end(), 0.0) / eval_stats.average_global_best_fitnesses.size();
    
    double average_fitness = total_fitness / filtered_sol.size();
    double average_gather_time = eval_stats.total_gather_time / eval;
    double average_scatter_time = eval_stats.total_scatter_time / eval;
    double average_migrate_time = eval_stats.total_migrate_time / eval;
    
    double eval_duration = ( std::clock() - eval_stats.eval_start ) / (double) CLOCKS_PER_SEC;
    
    //std::fprintf(config::stats_out, "%d,%d,%4.10f,%4.10f,%4.10f,%4.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%d,%2.10f,%d,%d,%d,%d,%d,%2.10f\r\n", run, eval, average_fitness, eval_stats.local_best_fitness, eval_stats.global_best_fitness, eval_stats.average_local_best_fitness, eval_stats.average_global_best_fitness, average_scatter_time, average_gather_time, average_migrate_time, eval_stats.init_duration, eval_duration, average_topo_fitness, eval_stats.local_best_topo_fitness, eval_stats.global_best_topo_fitness, eval_stats.average_local_best_topo_fitness, eval_stats.average_global_best_topo_fitness, filtered[0].id, filtered[0].round_fitness, filtered[0].rounds, filtered[0].channel_count, t.id, t.rounds, t.channel_count, t.fitness);
    
    
    std::fprintf(config::sol_stats_out, "%d," "%d,"  "%3.10f,"         "%3.10f,"                      "%3.10f,"                       "%13.10f,"                             "%3.10f,",
                                    run,  eval,  average_fitness,  eval_stats.local_best_fitness, eval_stats.global_best_fitness, eval_stats.average_local_best_fitness, eval_stats.average_global_best_fitness);

    
    std::fprintf(config::sol_stats_out, "%3.10f,"             "%3.10f,"            "%3.10f,"             "%3.10f,"               "%3.10f\r\n",
                                    average_scatter_time, average_gather_time, average_migrate_time, run_stats.init_duration, eval_duration);
    
    
    LOG(3, 0, 0, "%-2s,%-4s,%-13s,%-13s,%-13s,%-13s,%-13s,%-12s,%-12s,%-12s,%-12s,%-12s,%-13s,%-13s,%-13s,%-13s,%-13s,%-13s,%-13s\r\n", "r", "e", "avg_fitness", "lb_fitness", "gb_fitness", "alb_fitness", "agb_fitness", "avg_scatter", "avg_gather", "avg_migrate", "init_time", "eval_time", "avgt_fit", "lbt_fit", "gbt_fit", "albt_fit", "agbt_fit", "round_fit", "rounds");
    
    
    
    LOG(2, 0, 0, "\r\n%03d," "%05d,"  "%013.10f,"         "%3.10f,"                      "%3.10f,"                       "%13.10f,"                              "%3.10f,"                               "%3.10f,"             "%3.10f,"            "%3.10f,"             "%3.10f,"               "%013.10f"     "%s",
                 run,   eval,   average_fitness,    eval_stats.local_best_fitness, eval_stats.global_best_fitness, eval_stats.average_local_best_fitness, eval_stats.average_global_best_fitness, average_scatter_time, average_gather_time, average_migrate_time, run_stats.init_duration, eval_duration, "|");
    
    
}

#endif /* stats_h */

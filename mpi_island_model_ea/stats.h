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

void log_pop_stats(solver &solver, island &isle, MPI_Datatype &visa_type) {
    
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
        
        // top 20 individuals to cut log size
        
        for (auto it = solver.solutions.population.begin(); it !=solver.solutions.population.begin() + 20; ++it) {
            
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
    
            char sol_id[sizeof(it->id)*sizeof(char)];
            strcpy(sol_id, it->id);
            
            std::vector<visa> sol_visas;
            
            std::copy_if(visas.begin(), visas.end(), std::back_inserter(sol_visas), [&](visa &v) { return strcmp(sol_id, v.genome_id) == 0; });
            
            char vis[sizeof(int)*sizeof(sol_visas)+sizeof(char)*sizeof(sol_visas)];
            int offv = 0;
            
            for(int v=0; v<sol_visas.size(); v++) {
                char delim = (v == sol_visas.size()-1) ? '\0' : ';';
                offv += snprintf(vis+offv, sizeof(vis)>offv?sizeof(vis)-offv:0, "%d%c", sol_visas[v].destination, delim);
            }
            
            std::fprintf(config::solpop_out, "%d," "%d," "%s," "%d," "%d," "%s," "%s," "%d," "%d,", solver.solutions.run.id, solver.solutions.run.eval.id, it->id, it->source, it->locale, it->parents[0], it->parents[1], it->selected, it->survival);
                         
            std::fprintf(config::solpop_out, "%lld," "%d,"  "%lld," "%d," "%ld," "%d,", e5, groups_e5[e5], e4, groups_e4[e4], e3, groups_e3[e3]);
                         
            std::fprintf(config::solpop_out, "%f," "%f," "%d," "%s," "%s\r\n", it->fitness, it->selection_distribution, it->migrations, vis, sol);
            
        }
        
    }
    
}

void log_fn_topology_stats(solver &solver, meta &meta, topology &t) {
    
    char msg[24];
    
    sprintf(msg, config::ea_mode > 0 ? "EXPERIMENTAL" : "BENCHMARK");
    
    LOG(2, solver.variant.isle.id, 0, "\r\n\r\n    --- META GENOME %d (%s) EVALUATION %d of %d AT SOLVER[%d,%d] migrations = %d | duration = %f | ", t.id, msg, t.rounds, meta.topologies.max_fit_evals, solver.solutions.run.id, solver.solutions.run.eval.id, t.rounds, meta.topologies.run.eval.stats.eval_duration);
    
    LOG(3, 0, 0, "STATS (run %d, eval %d): topologies size %lu, mem[0] fit = %f \r\n", meta.topologies.run.id, meta.topologies.run.eval.id, meta.topologies.population.size(), meta.topologies.population[0].fitness);
    
    std::vector<topology> filtered_top;
    
    std::copy_if( meta.topologies.population.begin(), meta.topologies.population.end(), std::back_inserter(filtered_top), [](const topology &item) { return item.fitness != 0.00000000000 && item.rounds >= config::ea_2_max_fit_evals; });
    std::sort(filtered_top.begin(), filtered_top.end(), compare_fitness<topology>);
    std::reverse(filtered_top.begin(), filtered_top.end());
    
    if(filtered_top.size() == 0) { filtered_top.push_back(meta.topologies.population[0]); }
    
    LOG(3, 0, 0, " *** STATS (run %d, eval %d): filtered topologies size %lu, mem[0] fit = %f \r\n", meta.topologies.run.id, meta.topologies.run.eval.id, filtered_top.size(), filtered_top[0].fitness);
    
    meta.run.eval.stats.local_best_topo_fitness = filtered_top[0].fitness;
    meta.run.eval.stats.average_local_best_topo_fitnesses.push_back(filtered_top[0].fitness);

    LOG(3, 0, 0, " *** STATS (run %d, eval %d): best topology fit = %f", meta.topologies.run.id, meta.topologies.run.eval.id, meta.run.eval.stats.local_best_topo_fitness);
    
    if(filtered_top[0].fitness > meta.run.eval.stats.global_best_topo_fitness || meta.run.eval.stats.global_best_topo_fitness == 0.0) {
        LOG(3, 0, 0, " *** STATS (run %d, eval %d): found new global best topology %f > %f", meta.topologies.run.id, meta.topologies.run.eval.id, meta.topologies.population[0].fitness, meta.run.eval.stats.global_best_topo_fitness);
        meta.run.eval.stats.best_topology = filtered_top[0];
        meta.run.eval.stats.topo_best_count++;
        meta.run.eval.stats.average_global_best_topo_fitnesses.push_back(filtered_top[0].fitness);
        meta.run.eval.stats.global_best_topo_fitness = filtered_top[0].fitness;
        log_topology_matrix(solver.variant.islands, filtered_top[0], meta.run.eval.stats.topo_best_count);
    }
    
    double total_topo_fitness = 0.0;
    
    for(int i=0; i<filtered_top.size(); i++) {
        total_topo_fitness += filtered_top[i].fitness;
    }
    
    LOG(3, 0, 0, " *** STATS (run %d, eval %d): total topology fit = %f", meta.topologies.run.id, meta.topologies.run.eval.id, total_topo_fitness);
    
    meta.run.eval.stats.average_local_best_topo_fitness = std::accumulate(meta.run.eval.stats.average_local_best_topo_fitnesses.begin(), meta.run.eval.stats.average_local_best_topo_fitnesses.end(), 0.0) / meta.run.eval.stats.average_local_best_topo_fitnesses.size();
    
    meta.run.eval.stats.average_global_best_topo_fitness = std::accumulate(meta.run.eval.stats.average_global_best_topo_fitnesses.begin(), meta.run.eval.stats.average_global_best_topo_fitnesses.end(), 0.0) / meta.run.eval.stats.average_global_best_topo_fitnesses.size();
    
    double average_topo_fitness = filtered_top.size() > 0 ? (total_topo_fitness / filtered_top.size()) : 0.0;
    
    std::fprintf(config::topo_stats_out, "%3.10f,"             "%d,"            "%d,"               "%d,"                       "%3.10f,"                   "%3.10f,"            "%3.10f,"                            "%3.10f,"                           "%3.10f," ,
                 average_topo_fitness, filtered_top[0].id,  filtered_top[0].rounds, filtered_top[0].channel_count,  filtered_top[0].round_fitness,  filtered_top[0].fitness, meta.run.eval.stats.local_best_topo_fitness,  meta.run.eval.stats.global_best_topo_fitness, meta.run.eval.stats.average_local_best_topo_fitness);

    std::fprintf(config::topo_stats_out, "%3.10f,"                                    "%d,"       "%d,"       "%d,"              "%3.10f\r\n",
                 meta.run.eval.stats.average_global_best_topo_fitness, t.id,       t.rounds,   t.channel_count,    t.fitness);
        
    LOG(2, 0, 0, "fitness=%f | channels=%d | time=%f || global best=%d | channels=%d | fitness=%f\r\n\r\n", t.fitness, t.channel_count, t.total_migration_time, filtered_top[0].id, filtered_top[0].channel_count, filtered_top[0].fitness);
    
    LOG(3, 0, 0, "\r\n\r\n%5s %9s %10s %13s %13s %13s", "t_id", "t_evals", "t_channels", "t_fitness", "t_mig_time", "sum_mig_time");
    
    LOG(3, 0, 0, "\r\n%0d,"     "%9d,"     "%10d,"             "%16.11f,"    "%16.11f,"            "%16.11f,",
                 t.id,       t.rounds,   t.channel_count,    t.fitness,     t.total_migration_time, meta.topologies.run.eval.stats.total_migrate_time);
    
    LOG(4, 0, 0, "\r\n\r\n%13s %12s %13s %13s %13s %13s %13s %13s %13s %13s", "avg_t_fit", "gbest_t_id", "gbest_t_rounds", "gbest_t_chan", "gbest_t_efit", "gbest_tfit1", "lbest_tfit", "gbest_tfit2", "avg_lbest_tfit", "avg_gbest_tfit");
    
    LOG(3, 0, 0, "\r\n%3.10f,"          "%012d,"               "%014d,"                 "%013d,"                         "%3.10f,"                       "%3.10f,"               "%3.10f,",
                 average_topo_fitness,   filtered_top[0].id,  filtered_top[0].rounds, filtered_top[0].channel_count,  filtered_top[0].round_fitness,  filtered_top[0].fitness, meta.topologies.run.eval.stats.local_best_topo_fitness);
    
    LOG(3, 0, 0, "%s"   "%3.10f,"                           "%s,"  "%3.10f,"                                    "%3.10f,"                                   "%s",
                 BLU,   meta.run.eval.stats.global_best_topo_fitness, YEL, meta.run.eval.stats.average_local_best_topo_fitness, meta.run.eval.stats.average_global_best_topo_fitness, RESET);
    
}


void log_fn_eval_stats(solver &solver, meta &meta, topology &t) {
    
    // only consider evaluated topologies for stats
    
    std::vector<solution> filtered_sol;
    
    std::copy_if(solver.solutions.population.begin(), solver.solutions.population.end(), std::back_inserter(filtered_sol), [](const solution &item) { return item.fitness != 0.0; });
    std::sort(filtered_sol.begin(), filtered_sol.begin(), compare_fitness<solution>);
    std::reverse(filtered_sol.begin(), filtered_sol.end());
    
    LOG(3, 0, 0, "STATS (run %d, eval %d): solutions size %lu, mem[0] fit = %f \r\n", solver.solutions.run.id, solver.solutions.run.eval.id, solver.solutions.population.size(), solver.solutions.population[0].fitness);
    
    solver.run.eval.stats.local_best_fitness = filtered_sol[0].fitness;
    solver.run.eval.stats.average_local_best_fitnesses.push_back(filtered_sol[0].fitness);
        
    LOG(3, 0, 0, "STATS (run %d, eval %d): local best = %f, local best fitnesses tracked = %lu\r\n", solver.solutions.run.id, solver.solutions.run.eval.id, solver.run.eval.stats.local_best_fitness, solver.run.eval.stats.average_local_best_fitnesses.size());
    
    if(filtered_sol[0].fitness >solver.run.eval.stats.global_best_fitness) {
        LOG(3, 0, 0, "STATS (run %d, eval %d): found new global best solution %f > %f\r\n", solver.solutions.run.id, solver.solutions.run.eval.id, filtered_sol[0].fitness,solver.run.eval.stats.global_best_fitness);
       solver.run.eval.stats.sol = &filtered_sol[0];
       solver.run.eval.stats.average_global_best_fitnesses.push_back(filtered_sol[0].fitness);
       solver.run.eval.stats.global_best_fitness = filtered_sol[0].fitness;
    }
    
    double total_fitness = 0.0;
    
    for(int i=0; i<filtered_sol.size(); i++) {
        total_fitness += filtered_sol[i].fitness;
    }
    
    LOG(3, 0, 0, "STATS (run %d, eval %d): solution population size %lu total fitness = %2.10f\r\n", solver.solutions.run.id, solver.solutions.run.eval.id, filtered_sol.size(), total_fitness);
    
    if(solver.run.eval.stats.average_local_best_fitnesses.size() > 1) {
       solver.run.eval.stats.average_local_best_fitness = std::accumulate(solver.run.eval.stats.average_local_best_fitnesses.begin(),solver.run.eval.stats.average_local_best_fitnesses.end(), 0.0) /solver.run.eval.stats.average_local_best_fitnesses.size();
    } else {
       solver.run.eval.stats.average_local_best_fitness = filtered_sol[0].fitness;
    }
    
   solver.run.eval.stats.average_global_best_fitness = std::accumulate(solver.run.eval.stats.average_global_best_fitnesses.begin(),solver.run.eval.stats.average_global_best_fitnesses.end(), 0.0) /solver.run.eval.stats.average_global_best_fitnesses.size();
    
    double average_fitness = total_fitness / filtered_sol.size();
    double average_gather_time = solver.run.eval.stats.total_gather_time / solver.solutions.run.eval.id;
    double average_scatter_time = solver.run.eval.stats.total_scatter_time / solver.solutions.run.eval.id;
    double average_migrate_time = solver.run.eval.stats.total_migrate_time / solver.solutions.run.eval.id;
    
    std::fprintf(config::sol_stats_out, "%d," "%d," "%3.10f," "%3.10f," "%3.10f," "%13.10f," "%3.10f,", solver.solutions.run.id, solver.solutions.run.eval.id,  average_fitness, solver.run.eval.stats.local_best_fitness, solver.run.eval.stats.global_best_fitness, solver.run.eval.stats.average_local_best_fitness, solver.run.eval.stats.average_global_best_fitness);

    
    std::fprintf(config::sol_stats_out, "%f," "%f," "%f," "%f," "%f\r\n", average_scatter_time, average_gather_time, average_migrate_time, solver.init_duration, solver.run.eval.stats.eval_duration);
        
    
//    LOG(2, 0, 0, "\r\n%5d," "%6d,"  "%16.11f," "%16.11f," "%16.11f," "%16.11f," "%16.11f," "%16.11f," "%16.11f,""%16.11f," "%16.11f", solver.solutions.run.id, solver.solutions.run.eval.id, average_fitness, solver.run.eval.stats.local_best_fitness, solver.run.eval.stats.global_best_fitness, solver.run.eval.stats.average_local_best_fitness, solver.run.eval.stats.average_global_best_fitness, average_scatter_time, average_gather_time, average_migrate_time, solver.run.eval.stats.eval_duration);
    
    LOG(2, 0, 0, "\r\n%5d," "%6d,"  "%16.11f," "%16.11f," "%16.11f," "%16.11f," "%16.11f," "%16.11f," "%16.11f,""%16.11f," "%16.11f", solver.solutions.run.id, solver.solutions.run.eval.id, average_fitness, solver.run.eval.stats.local_best_fitness, solver.run.eval.stats.global_best_fitness, solver.run.eval.stats.average_local_best_fitness, solver.run.eval.stats.average_global_best_fitness, average_scatter_time, average_gather_time, average_migrate_time, solver.solutions.run.eval.stats.eval_duration);
    
    
}

void log_run_stats(meta &meta, objective<topology> &o) {

    if(meta.variant.isle.id != 0) { return; }
    
    LOG(2, meta.variant.isle.id, 0, "\r\n\r\n--- END META RUN %d ---\r\n\r\n", meta.topologies.run.id);

    o.run.end();
    
    fprintf(config::topo_run_stats_out, "average_topo_fitness, global_best_topo_id, global_best_topo_rounds, global_best_topo_channels, global_best_topo_round_fitness, global_best_topo_fitness1, local_best_topo_fitness, global_best_topo_fitness2, average_local_best_topo_fitness, average_global_best_topo_fitness, t_id, t_rounds, t_channels, t_fitness\r\n");

    std::fprintf(config::topo_run_stats_out, "%d,%f,%f,%f,%f,%f,%d", meta.topologies.run.id, meta.topologies.run.stats.run_duration, meta.run.eval.stats.average_local_best_topo_fitness, meta.run.eval.stats.average_global_best_topo_fitness, meta.run.eval.stats.global_best_topo_fitness, meta.run.eval.stats.total_migrate_time, meta.run.stats.total_channels);

    fflush(config::topo_run_stats_out);

}

#endif /* stats_h */

//
//  dtype_stats.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 10/13/21.
//  Copyright © 2021 Bradley Morgan. All rights reserved.
//

#ifndef dtype_stats_h
#define dtype_stats_h

#include <mpi.h>

#define RED   "\x1B[31m"
#define GRN   "\x1B[32m"
#define YEL   "\x1B[33m"
#define BLU   "\x1B[34m"
#define MAG   "\x1B[35m"
#define CYN   "\x1B[36m"
#define WHT   "\x1B[37m"
#define RESET "\x1B[0m"

struct estats {
    
    int topo_best_count = 0;
    
    double init_duration = 0.0;
    double eval_start = 0.0;
    double eval_duration = 0.0;
    double cycle_duration = 0.0;
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
    
    solution *sol;
    
    void init() {
        
        LOG(3, 0, 0, "INIT EA STATS\r\n");
        
        this->init_duration = 0.0;
        this->eval_start = 0.0;
        this->eval_duration = 0.0;
        this->total_scatter_time = 0.0;
        this->total_gather_time = 0.0;
        this->total_migrate_time = 0.0;
        this->topo_migrate_time = 0.0;
        this->local_best_fitness = 0.0;
        this->global_best_fitness = 0.0;
        this->average_local_best_fitness = 0.0;
        this->average_global_best_fitness = 0.0;
        this->local_best_topo_fitness = 0.0;
        this->average_local_best_topo_fitness = 0.0;
        this->average_local_best_fitnesses.clear();
        this->average_local_best_fitnesses.resize(0);
        this->average_local_best_topo_fitnesses.clear();
        this->average_local_best_topo_fitnesses.resize(0);
        
    }
    
//    template<typename f, typename e, typename o, typename genome> void log(e &ea, e &meta, o &obj, genome &g, f function) { function(ea, meta, g); }
//
//    template<typename f, typename e, typename o> void log(e &ea, o &obj, f function) { function(ea, obj, *this); }
    
};

struct rstats {

    double run_start;
    double run_duration;
    double init_duration;
    
    int total_channels;
    
    void init() {
        
        LOG(3, 0, 0, "INIT RUN STATS\r\n");
        
        this->run_start = 0.0;
        this->run_duration = 0.0;
        this->init_duration = 0.0;
        this->total_channels = 0;
        
    }
    
    //template<typename f, typename o, typename e> void log(e &ea, o &obj, f function) { function(ea, obj); }
    
};

#pragma mark DATATYPE: @ea_eval{}

// track eval stats, etc. per island

struct ea_eval {
  
    int id;
    estats stats;
  
    ea_eval(): id(0) {};
    
    void begin() {
        LOG(3, 0, 0, "BEGIN EA EVAL %d -> ", this->id+1);
        this->stats.eval_start = MPI_Wtime();
    }
    
    void end() {
        double eval_end = MPI_Wtime();
        this->stats.eval_duration = eval_end - this->stats.eval_start;
        LOG(3, 0, 0, "END EA EVAL %d\r\n", this->id);
    }
    
    //template<typename e> void log(e &ea) { }
    
};

#pragma mark EA::DATATYPE: @ea_run{}

// track run stats, etc. per island

struct ea_run {
    
    int id;
    rstats stats;
    ea_eval eval;
  
    ea_run(): id(0) {};
    
    void begin() {
        LOG(3, 0, 0, "BEGIN EA RUN %d -> ", this->id);
        this->stats.run_duration = 0.0;
        this->stats.run_start = MPI_Wtime();
    }
    
    void end() {
        LOG(3, 0, 0, "END EA RUN %d\r\n", this->id);
        double run_end = MPI_Wtime();
        this->stats.run_duration = run_end - this->stats.run_start;
    }
    
};

#pragma mark EA::OBJECTIVE::DATATYPE: @objective_eval{}

// track eval stats, etc. per island

struct objective_eval {
  
    int id = 0;
    estats stats;
    
    objective_eval(): id(0) {};
    
    void begin() {
        LOG(3, 0, 0, "BEGIN OBJECTIVE EVAL %d\r\n", this->id);
        //this->stats.init();
        this->stats.eval_duration = 0.0;
        this->stats.eval_start = MPI_Wtime();
    }
    
    void end() {
        LOG(3, 0, 0, "END OBJECTIVE EVAL %d\r\n", this->id);
        double eval_end = MPI_Wtime();
        this->stats.eval_duration = eval_end - this->stats.eval_start;
    }
    
};

#pragma mark EA::OBJECTIVE::DATATYPE: @objective_run{}

// track run stats, etc. per island

struct objective_run {
    
    int id = 0;
    
    rstats stats;
    objective_eval eval;

    objective_run(): id(0) {};
    
    void begin() {
        LOG(3, 0, 0, "BEGIN OBJECTIVE RUN %d -> ", this->id);
        this->eval.stats.init();
        this->stats.run_duration = 0.0;
        this->stats.run_start = MPI_Wtime();
    }
    
    void end() {
        LOG(3, 0, 0, "END OBJECTIVE RUN %d\r\n", this->id);
        double run_end = MPI_Wtime();
        this->stats.run_duration = run_end - this->stats.run_start;
    }
    
};

struct evolution_cycle {
    
    int id = 0;
    double start = 0.0;
    double duration = 0.0;
    
    void begin() {
        LOG(3, 0, 0, "BEGIN OBJECTIVE CYCLE %d -> ", this->id);
        this->duration = 0.0;
        this->start = MPI_Wtime();
    }
    
    void end() {
        LOG(3, 0, 0, "END OBJECTIVE CYCLE %d\r\n", this->id);
        double cycle_end = MPI_Wtime();
        this->duration = cycle_end - this->start;
    }
    
    
};

#endif /* dtype_stats_h */


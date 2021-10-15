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
    
    template<typename f, typename e, typename o, typename genome> void log(e &ea, e &meta, o &obj, genome &g, f function) { function(ea, meta, g); }
    template<typename f, typename e, typename o> void log(e &ea, o &obj, f function) { function(ea, obj, *this); }
    
};

struct rstats {

    double run_start;
    double run_duration;
    double init_duration;
    
    int total_channels;
    
    void init() {
        
        this->run_start = 0.0;
        this->run_duration = 0.0;
        this->init_duration = 0.0;
        this->total_channels = 0;
        
    }
    
    template<typename f, typename o, typename e> void log(e &ea, o &obj, f function) { function(ea, obj); }
    
};

#pragma mark DATATYPE: @ea_eval{}

// track eval stats, etc. per island

struct ea_eval {
  
    int id;
    double start;
    estats stats;
  
    ea_eval(): id(0) {};
    
    void begin() {
        this->start = MPI_Wtime();
    }
    
    void end() {
        double eval_end = MPI_Wtime();
        this->stats.eval_duration = eval_end - this->stats.eval_start;
    }
    
//    template<typename o, typename e> void begin(e &ea, o &obj) { this->begin_eval(); }
//    template<typename o, typename e> void end(e &ea, o &obj) { this->end_eval(); }
    
//    template<typename f, typename o, typename e> void begin(e &ea, o &obj, f function) { this->begin_eval(); this->stats.log(ea, obj, function); }
//    template<typename f, typename o, typename e> void end(e &ea, o &obj, f function) { this->end_eval(); this->stats.log(ea, obj, function); }
    

    
    template<typename e> void log(e &ea) { }
    
};

#pragma mark DATATYPE: @ea_run{}

// track run stats, etc. per island

struct ea_run {
    
    int id;
    double start;
    rstats stats;
    ea_eval eval;
  
    ea_run(): id(0) {};
    
    void begin_t() {
        this->start = MPI_Wtime();
        this->stats.init();
        this->eval.stats.init();
    }
    
    void end_t() {
        double run_end = MPI_Wtime();
        this->stats.run_duration = run_end - this->stats.run_start;
    }
    
    template<typename o, typename e> void begin(e &ea, o &obj) { this->begin_t(); }
    template<typename o, typename e> void end(e &ea, o &obj) { this->end_t(); }
    
};

#pragma mark DATATYPE: @objective_eval{}

// track eval stats, etc. per island

struct objective_eval {
  
    int id = 0;
    double start;
    
    objective_eval(): id(0) {};
    
    template<typename f, typename e, typename g> void log(e &ea, e &meta, g &genome, f function) { function(ea, meta, genome); }
    
};

#pragma mark DATATYPE: @objective_run{}

// track run stats, etc. per island

struct objective_run {
    
    int id = 0;
    double start;
    objective_eval eval;
    
    objective_run(): id(0) {};
        
    template<typename f, typename o, typename e, typename genome> void log(e &ea, e &meta, o &obj, genome &g, f function) { function(obj, ea, meta, g); }
//    template<typename f, typename e, typename o> void log(e &ea, e &meta, o &obj, f function) { function(ea, meta, obj); }
//    template<typename f, typename o, typename e> void log(e &ea, o &obj, f function) { function(ea, obj); }
    
};


#endif /* dtype_stats_h */


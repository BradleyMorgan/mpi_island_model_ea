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

struct mpi_local {
    double value;
    int island = mpi.id;
    
    void init() { mpi_local(); }
    mpi_local() { value=0.0; island=mpi.id; }
};

struct time_stats {
    
    double start = 0.0;
    double duration = 0.0;
    double init_duration = 0.0;
    
    mpi_local local_t;
    mpi_local min_t;
    mpi_local max_t;
    mpi_local sum_t;
    double avg_t = 0.0;
    
};

struct parallel_stats {
    
    mpi_local local_migration_t;
    mpi_local min_migration_t;
    mpi_local max_migration_t;
    mpi_local sum_migration_t;
    
    double avg_migration_t = 0.0;
    
    mpi_local local_gather_t;
    mpi_local min_gather_t;
    mpi_local max_gather_t;
    mpi_local sum_gather_t;
   
    double avg_gather_t = 0.0;
    
    mpi_local local_scatter_t;
    mpi_local min_scatter_t;
    mpi_local max_scatter_t;
    mpi_local sum_scatter_t;
    
    double avg_scatter_t = 0.0;
    
};

struct fitness_stats {
    
    std::vector<double> best;
    
    double total_fitness = 0.0;
    double avg_fitness = 0.0;
    double min_fitness = 0.0;
    double max_fitness = 0.0;
    double sum_fitness = 0.0;
    
    double best_fitness = 0.0;
    double avg_best_fitness = 0.0;
    
};

struct eval_stats: time_stats, fitness_stats {
    
    double init_duration = 0.0;
    double start = 0.0;
    double duration = 0.0;
    
    void init() {
        
        LOG(5, 0, 0, "INIT EA STATS ");
        
        this->min_t.value = 0.0;
        this->max_t.value = 0.0;
        this->sum_t.value = 0.0;
        this->avg_t = 0.0;
        
        this->init_duration = 0.0;
        this->start = 0.0;
        this->duration = 0.0;
        
    }
    
};

struct cycle_stats : time_stats, parallel_stats, fitness_stats {

    std::vector<double> best;
    
    double best_fitness = 0.0;
    double avg_best_fitness = 0.0;
    
    void init() {
        
        LOG(5, 0, 0, "INIT CYCLE STATS ");
        
        this->start = 0.0;
        this->duration = 0.0;
        
    }
    
};

struct run_stats : time_stats, parallel_stats, fitness_stats {
    
    std::vector<double> best;
    
    double best_fitness = 0.0;
    double avg_best_fitness = 0.0;
    
    void init() {
        
        LOG(5, 0, 0, "INIT RUN STATS ");
        
        this->start = 0.0;
        this->duration = 0.0;
        this->init_duration = 0.0;
        
    }
    
};


#endif /* dtype_stats_h */


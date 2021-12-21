//
//  dtype_meta.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 10/17/21.
//  Copyright © 2021 Bradley Morgan. All rights reserved.
//

#ifndef dtype_meta_h
#define dtype_meta_h


template<> struct obj_eval_stats<topology> : time_stats, parallel_stats, fitness_stats {
    
    //obj_cycle_stats<topology> *cycle_stats;
    
    double min = 0.0;
    double max = 0.0;
    double sum = 0.0;
    
    void begin_header() {
        LOG(2, mpi.id, 0, "\n%s %s %s\n\n", std::string(6,'~').c_str(), "BEGIN META EA OBJECTIVE EVAL", std::string(98,'~').c_str());
    }
    
    void end_header() {
        LOG(2, mpi.id, 0, "\n\n%s %s %s\n", std::string(6,'~').c_str(), "END META EA OBJECTIVE EVAL", std::string(100,'~').c_str());
        LOG(2, mpi.id, 0, "\n%4s %6s %8s %11s %11s %11s %11s %11s %11s %14s %14s %11s\n", "r", "c", "E", "avg_fit", "gb_fit", "sm_scat", "sm_gat", "sm_mig", "sum_t", "min_t", "max_t", "avg_t");
    }
    
    obj_eval_stats<topology>() : time_stats(), parallel_stats(), fitness_stats() {}
    
};

template<> struct obj_cycle_stats<topology> : time_stats, parallel_stats, fitness_stats {
    
    //obj_run_stats<topology> *run_stats;
    
    topology *global_best;
    
    int total_channels = 0;
    
    void init() {
        
        LOG(5, 0, 0, "INIT RUN STATS\r\n");
        
        this->total_channels = 0;
        this->start = 0.0;
        this->duration = 0.0;
        
    }
    
    void begin_header() {
        LOG(2, mpi.id, 0, "\n%s %s %s\n", std::string(6,'~').c_str(), "BEGIN META EA OBJECTIVE CYCLE", std::string(104,'~').c_str());
        //LOG(2, mpi.id, 0, "\n%4s %6s %8s %11s %11s %11s %11s %11s %11s %14s %14s %11s\n", "R", "c", "e", "avg_fit", "gb_fit", "sm_scat", "sm_gat", "sm_mig", "sum_t", "min_t", "max_t", "avg_t");
    }
    
    void end_header() {
        LOG(2, mpi.id, 0, "\n%s %s %s\n", std::string(6,'~').c_str(), "END META EA OBJECTIVE CYCLE", std::string(104,'~').c_str());
    }
    
    obj_cycle_stats<topology>() : time_stats(), parallel_stats(), fitness_stats() {}
    
};

template<> struct obj_run_stats<topology> : time_stats, parallel_stats, fitness_stats {
    
    int total_channels = 0;
    
    topology *global_best;
    
    std::vector<double> best;
    
    void begin_header() {
        LOG(2, mpi.id, 0, "\n%s %s %s\n", std::string(6,'~').c_str(), "BEGIN META EA OBJECTIVE RUN", std::string(104,'~').c_str());
        //LOG(2, mpi.id, 0, "\n%4s %6s %8s %11s %11s %11s %11s %11s %11s %14s %14s %11s\n", "R", "c", "e", "avg_fit", "gb_fit", "sm_scat", "sm_gat", "sm_mig", "sum_t", "min_t", "max_t", "avg_t");
    }
    
    void end_header() {
        LOG(2, mpi.id, 0, "\n%s %s %s\n", std::string(6,'~').c_str(), "END META EA OBJECTIVE RUN", std::string(104,'~').c_str());
    }
    
    obj_run_stats<topology>() : time_stats(), parallel_stats(), fitness_stats() {}
    
};

struct ea_meta : ea<ea_meta> {
    
    objective<topology> topologies;
    
    template<typename variant> ea_meta(variant &target) {
        
        this->model = target.model;
        
        sprintf(this->name, "%s", config::ea_2_name);
        sprintf(this->topologies.name, "%s", config::ea_2_o1_name);
        
        LOG(2, mpi.id, 0, "%s EA INIT\r\n", this->name);
        
        // termination
        
        this->run.max = config::ea_2_max_runs;
        this->run.cycle.max = config::ea_2_max_cycles;
        this->run.cycle.eval.max = config::ea_1_max_evals;
        
        // logging
        
        this->run.log_interval = config::ea_2_log_run_interval;
        this->run.cycle.log_interval = config::ea_2_log_cycle_interval;
        this->run.cycle.eval.log_interval = config::ea_2_log_eval_interval;
        
        this->topologies.run.max = config::ea_2_o1_max_runs;
        this->topologies.run.cycle.max = config::ea_2_o1_max_cycles;
        this->topologies.run.cycle.eval.max = config::ea_2_o1_max_evals;
        this->topologies.log_interval = config::ea_2_log_population_interval;
        
        this->topologies.run.stats_out = config::ea_2_run_out;
        this->topologies.run.cycle.stats_out = config::ea_2_cycle_out;
        this->topologies.run.cycle.eval.stats_out = config::ea_2_eval_out;
        
        this->topologies.run.cycle.eval.log_stdout = true;
        this->topologies.run.log_interval = config::ea_2_o1_log_run_interval;
        this->topologies.run.cycle.log_interval = config::ea_2_o1_log_cycle_interval;
        this->topologies.run.cycle.eval.log_interval = config::ea_2_o1_log_eval_interval;
        
        this->topologies.mu = config::ea_2_mu;
        this->topologies.lambda = config::ea_2_lambda;
        this->topologies.mutation_rate = config::ea_2_mutation_rate;
        
        this->stats.init_duration += (MPI_Wtime() - this->stats.start);
        
    }
    
    template<typename e> void begin(e &target);
    
};


#endif /* dtype_meta_h */

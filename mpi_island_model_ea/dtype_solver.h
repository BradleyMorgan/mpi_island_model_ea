//
//  dtype_solver.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 10/17/21.
//  Copyright © 2021 Bradley Morgan. All rights reserved.
//

#ifndef dtype_solver_h
#define dtype_solver_h

template<> struct obj_run_stats<solution> : time_stats, parallel_stats, fitness_stats {
    
    std::vector<double> best;
    
    template<typename stats> void cascade(double stats::*field, double val, char op) {
        
        switch (op) {
            case '<': this->*field = this->*field < val ? this->*field : val;
            case '>': this->*field = this->*field > val ? this->*field : val;
            case '+': this->*field += val;
            case '=': this->*field = val;
        }

    }
    
    void begin_header() {
        LOG(2, mpi.id, 0, "\n%s %s %s\n", std::string(6,'-').c_str(), "SOLVER EA OBJECTIVE RUN", std::string(103,'-').c_str());
        LOG(2, mpi.id, 0, "\n%4s %6s %8s %11s %11s %11s %11s %11s %11s %14s %14s %11s\n", "R", "c", "e", "avg_fit", "gb_fit", "sm_scat", "sm_gat", "sm_mig", "sum_t", "min_t", "max_t", "avg_t");
    }
    
    void end_header() {
        LOG(2, mpi.id, 0, "%s\n", std::string(134,'-').c_str());
    }
    
    obj_run_stats<solution>() : time_stats(), parallel_stats(), fitness_stats() {}
    
};

template<> struct obj_cycle_stats<solution> : time_stats, parallel_stats, fitness_stats {
    
    //obj_run_stats<solution> *run_stats;
    
    std::vector<double> best;
    
//    template<typename stats> void cascade(double stats::*field, double val, char op) {
//
//        switch (op) {
//            case '+': {
//                this->*field += val;
//                this->run_stats->*field += val;
//            }
//            case '=': {
//                this->*field = val;
//                this->run_stats->*field += val;
//            }
//        }
//
//    }
//
//    template<typename stats> void cascade(mpi_local stats::*field, mpi_local source, char op) {
//
//        switch (op) {
//                switch (op) {
//                    case '<':
//                        this->*field = (this->*field).value > source.value ? source : this->*field;
//                        this->run_stats->*field = (this->run_stats->*field).value > source.value ? source : this->run_stats->*field;
//                    case '>':
//                        this->*field = (this->*field).value < source.value ? source : this->*field;
//                        this->run_stats->*field = (this->run_stats->*field).value < source.value ? source : this->run_stats->*field;
//                    case '+':
//                        (this->*field).value += source.value;
//                        (this->run_stats->*field).value += source.value;
//                    case '=':
//                        (this->*field).value = source.value;
//                        (this->run_stats->*field).value = source.value;
//                }
//
//        }
//
//    }
    
    //obj_cycle_stats<solution>(obj_run_stats<solution> &r) : time_stats(), parallel_stats(), fitness_stats(), run_stats(&r) {}
    
    void begin_header() {
        //LOG(2, mpi.id, 0, "\n%s %s %s\n", std::string(6,'-').c_str(), "SOLVER EA OBJECTIVE CYCLE", std::string(103,'-').c_str());
        //LOG(2, mpi.id, 0, "\n%4s %6s %8s %11s %11s %11s %11s %11s %11s %14s %14s %11s\n", "R", "c", "e", "avg_fit", "gb_fit", "sm_scat", "sm_gat", "sm_mig", "sum_t", "min_t", "max_t", "avg_t");
    }
    
    void end_header() {
        //LOG(2, mpi.id, 0, "%s\n", std::string(134,'-').c_str());
    }
    
    obj_cycle_stats<solution>() : time_stats(), parallel_stats(), fitness_stats() {}
    
};

template<> struct obj_eval_stats<solution> : time_stats, parallel_stats, fitness_stats {
    
    //obj_cycle_stats<solution> *cycle_stats;
    
//    template<typename stats> void cascade(double stats::*field, double val, char op) {
//                    
//            switch (op) {
//                case '<': {
//                    this->*field = this->*field < val ? this->*field : val;
//                    double cval = this->cycle_stats->*field;
//                    cval = cval == 0.0 ? val : cval < val ? cval : val;
//                }
//                case '>': {
//                    this->*field = this->*field > val ? this->*field : val;
//                    double cval = this->cycle_stats->*field;
//                    cval = cval == 0.0 ? val : cval > val ? cval : val;
//                }
//                case '+': {
//                    this->*field += val;
//                    this->cycle_stats->*field += val;
//                }
//                case '=': {
//                    this->*field = val;
//                    this->cycle_stats->*field += val;
//                }
//            }
//        
//        this->cycle_stats->cascade(field, val, op);
//        
//    }
//
    
    void begin_header() {
        //LOG(2, mpi.id, 0, "\n%s %s %s\n", std::string(6,'-').c_str(), "SOLVER EA OBJECTIVE EVAL", std::string(103,'-').c_str());
        //LOG(2, mpi.id, 0, "\n%4s %6s %8s %11s %11s %11s %11s %11s %11s %14s %14s %11s\n", "R", "c", "e", "avg_fit", "gb_fit", "sm_scat", "sm_gat", "sm_mig", "sum_t", "min_t", "max_t", "avg_t");
    }
    
    void end_header() {
        //LOG(2, mpi.id, 0, "%s\n", std::string(134,'-').c_str());
    }
    
    obj_eval_stats<solution>() : time_stats(), parallel_stats(), fitness_stats() {}
    
};

struct ea_solver: ea<ea_solver> {
    
    std::array<double, DIM> offsets;
    
    solution loc;
    
    objective<solution> solutions = objective<solution>();
    
    ea_solver() {
        
        this->model.init();
    
        sprintf(this->name, "%s", config::ea_1_name);
        sprintf(this->solutions.name, "%s", config::ea_1_o1_name);
        this->solutions.aggregate = {};
        
        LOG(2, mpi.id, 0, "%s EA INIT\r\n", this->name);
        
        // ea global termination
        
        this->run.max = config::ea_1_max_runs;
        this->run.cycle.max = config::ea_1_max_cycles;
        this->run.cycle.eval.max = config::ea_1_max_evals;
        
        // ea global logging
        
        this->run.log_interval = config::ea_1_log_run_interval;
        this->run.cycle.log_interval = config::ea_1_log_cycle_interval;
        this->run.cycle.eval.log_interval = config::ea_1_log_eval_interval;
        
        // ea objective termination

        this->solutions.run.max = config::ea_1_o1_max_runs;
        this->solutions.run.cycle.max = config::ea_1_o1_max_cycles;
        this->solutions.run.cycle.eval.max = config::ea_1_o1_max_evals;
        
        // ea objective logging
        
        this->solutions.run.stats_out = config::ea_1_run_out;
        this->solutions.run.cycle.stats_out = config::ea_1_cycle_out;
        this->solutions.run.cycle.eval.stats_out = config::ea_1_eval_out;
        
        this->solutions.run.log_interval = config::ea_1_o1_log_run_interval;
        this->solutions.run.cycle.log_interval = config::ea_1_o1_log_cycle_interval;
        this->solutions.run.cycle.eval.log_interval = config::ea_1_o1_log_eval_interval;
        
        // ea objective evolution
        
        this->solutions.mu = config::ea_1_mu;
        this->solutions.lambda = config::ea_1_lambda;
        this->solutions.mutation_rate = config::ea_1_mutation_rate;
        this->solutions.log_interval = config::ea_1_log_population_interval;
        
        if(mpi.id == 0) {
            this->offsets = generate_offsets(-2.5, 2.5, .5);
        }
        
        MPI_Bcast(&this->offsets, DIM, MPI_DOUBLE, 0, this->model.tcomm);
        
        double init_end = MPI_Wtime();
        
        this->stats.init_duration = (init_end - this->stats.start);
        
    }
    
};

#endif /* dtype_solver_h */

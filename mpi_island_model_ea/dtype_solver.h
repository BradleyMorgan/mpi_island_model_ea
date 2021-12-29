//
//  dtype_solver.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 10/17/21.
//  Copyright © 2021 Bradley Morgan. All rights reserved.
//

#ifndef dtype_solver_h
#define dtype_solver_h

struct solution_stats : time_stats, parallel_stats, fitness_stats {
    
    bool log_head = false;
    bool log_tail = false;
    
    int head_interval = 0;
    int tail_interval = 0;
    
    solution_stats(bool head, bool tail) : log_head(head), log_tail(tail) {};
    
};

template<> struct obj_run_stats<solution> : solution_stats {
    
    void begin_header(char *head) {
        if (this->log_head == true) { log_separator(head, '='); }
        LOG(2, mpi.id, 0, "\n%4s %6s %8s %11s %11s %11s %11s %11s %11s %14s %14s %11s\n", "R", "c", "e", "avg_fit", "gb_fit", "sm_scat", "sm_gat", "sm_mig", "sum_t", "min_t", "max_t", "avg_t");
    }
    
    void end_header(char *tail) {
        if (this->log_tail == true) { log_separator(tail, '='); }
        LOG(2, mpi.id, 0, "\n%4s %6s %8s %11s %11s %11s %11s %11s %11s %14s %14s %11s\n", "R", "c", "e", "avg_fit", "gb_fit", "sm_scat", "sm_gat", "sm_mig", "sum_t", "min_t", "max_t", "avg_t");
    }
    
    obj_run_stats<solution>() : solution_stats(true, true) {}
    
};

template<> struct obj_cycle_stats<solution> : solution_stats {
    
    void begin_header(char *head) {
        if (this->log_head == true) { log_separator(head, '-'); }
    }
    
    void end_header(char *tail) {
        if (this->log_tail == true) { log_separator(tail, '-'); }
    }
    
    obj_cycle_stats<solution>() : solution_stats(false, false) {}
    
};

template<> struct obj_eval_stats<solution> : solution_stats {
    
    void begin_header(char *head) {
        if (this->log_head == true) { log_separator(head, '.'); }
    }
    
    void end_header(char *tail) {
        if (this->log_tail == true) { log_separator(tail, '.'); }
    }
    
    obj_eval_stats<solution>() : solution_stats(false, false) {}
    
};

struct ea_solver: ea<ea_solver> {
    
    std::array<double, DIM> offsets;

    objective<solution> solutions = objective<solution>(1);
    
    ea_solver() {
        
        this->model.init();
        this->model.isle.stats.log_interval = config::ea_1_o1_log_island_interval;
        
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
        this->solutions.run.cycle.eval.max_local = config::ea_1_o1_max_fitness_evals;
        
        // ea objective logging
        
        this->solutions.run.stats_out = config::ea_1_run_out;
        this->solutions.run.cycle.stats_out = config::ea_1_cycle_out;
        this->solutions.run.cycle.eval.stats_out = config::ea_1_eval_out;
        
        this->solutions.run.log_interval = config::ea_1_o1_log_run_interval;
        this->solutions.run.cycle.log_interval = config::ea_1_o1_log_cycle_interval;
        this->solutions.run.cycle.eval.log_interval = config::ea_1_o1_log_eval_interval;
        
        // log population every nth cycle
        this->solutions.run.cycle.log_population_interval = config::ea_1_o1_log_population_interval;
        
        // log current genome every nth run
        this->solutions.run.log_genome_interval = config::ea_1_o1_log_genome_interval;
        
        // ea objective evolution
        
        this->solutions.mu = config::ea_1_mu;
        this->solutions.lambda = config::ea_1_lambda;
        this->solutions.mutation_rate = config::ea_1_mutation_rate;
        
        if(mpi.id == 0) {
            this->offsets = generate_offsets(-2.5, 2.5, .5);
        }
        
        MPI_Bcast(&this->offsets, DIM, MPI_DOUBLE, 0, this->model.tcomm);
        
        double init_end = MPI_Wtime();
        
        this->stats.init_duration = (init_end - this->stats.start);
        
    }
    
    template<typename i> void log_population(i &interval);
    
};

template<> template<typename i> void objective<solution>::log_population(i &interval) {

    if(mpi.id==0) {
            
        // top 20 individuals to cut log size
        
        for (auto it = this->population.begin(); it != this->population.begin() + 20; ++it) {
            
            long rank = std::distance(this->population.begin(), it);
            
            std::fprintf(config::ea_1_population_out, "%d," "%d," "%d," "%s," "%ld," "%f," "%d," "%d," "%s," "%s," "%d," "%d," "%d," "%f\r\n", this->run.id, this->run.cycle.id, this->run.cycle.eval.id, it->id, rank, it->fitness, it->source, it->locale, it->parents[0], it->parents[1], it->selected, it->survival, it->migrations, it->selection_distribution);
            
        }
        
    }

}

#endif /* dtype_solver_h */

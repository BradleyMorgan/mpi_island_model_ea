//
//  dtype_meta.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 10/17/21.
//  Copyright © 2021 Bradley Morgan. All rights reserved.
//

#ifndef dtype_meta_h
#define dtype_meta_h

template<> struct obj_eval_stats<topology> : topology_stats {
    
    double min = 0.0;
    double max = 0.0;
    double sum = 0.0;
    
    int log_interval = 1;
    
    void begin_header(char *head) {
        if (this->log_head == true) {  log_separator(head, '_'); }
    }
    
    void end_header(char *tail) {
        if (this->log_tail == true) { log_separator(tail, '_'); }
        LOG(2, mpi.id, 0, "\n%4s %6s %8s %11s %11s %11s %11s %11s %11s %14s %14s %11s\n", "r", "c", "E", "avg_fit", "gb_fit", "sm_scat", "sm_gat", "sm_mig", "sum_t", "min_t", "max_t", "avg_t");
    }
    
    obj_eval_stats<topology>() : topology_stats(true, true) {}
    
};

template<> struct obj_cycle_stats<topology> : topology_stats {
    
    void init() {
        
        LOG(5, 0, 0, "INIT RUN STATS\r\n");
        
        this->total_channels = 0;
        this->start = 0.0;
        this->duration = 0.0;
        
    }
    
    void begin_header(char *head) {
        if (this->log_head == true) { log_separator(head, '*'); }
    }
    
    void end_header(char *tail) {
        if (this->log_tail == true) { log_separator(tail, '*'); }
    }
    
    obj_cycle_stats<topology>() : topology_stats(true, true) {}
    
};

template<> struct obj_run_stats<topology> : topology_stats {
    
    void begin_header(char *head) {
        if (this->log_head == true) { log_separator(head, '#'); }
        LOG(2, mpi.id, 0, "\n%4s %6s %8s %11s %11s %11s %11s %11s %11s %14s %14s %11s\n", "r", "c", "E", "avg_fit", "gb_fit", "sm_scat", "sm_gat", "sm_mig", "sum_t", "min_t", "max_t", "avg_t");
    }
    
    void end_header(char *tail) {
        if (this->log_tail == true) { log_separator(tail, '#'); }
    }
    
    obj_run_stats<topology>() : topology_stats(true, true) {}
    
};

struct ea_meta : ea<ea_meta> {
    
    objective<topology> topologies = objective<topology>(1);
    
    template<typename variant> ea_meta(variant &target) {
        
        this->model = target.model;
        this->model.isle.stats.log_interval = config::ea_2_o1_log_island_interval;
        
        sprintf(this->name, "%s", config::ea_2_name);
        sprintf(this->topologies.name, "%s", config::ea_2_o1_name);
        
        LOG(2, mpi.id, 0, "%s EA INIT\r\n", this->name);
        
        // ea global termination
        
        this->run.max = config::ea_2_max_runs;
        this->run.cycle.max = config::ea_2_max_cycles;
        this->run.cycle.eval.max = config::ea_1_max_evals;
        
        // ea objective termination
        
        this->topologies.run.max = config::ea_2_o1_max_runs;
        this->topologies.run.cycle.max = config::ea_2_o1_max_cycles;
        this->topologies.run.cycle.eval.max = config::ea_2_o1_max_evals;
        this->topologies.run.cycle.eval.max_local = config::ea_2_o1_max_fitness_evals;
        
        // ea global logging
        
        this->run.log_interval = config::ea_2_log_run_interval;
        this->run.cycle.log_interval = config::ea_2_log_cycle_interval;
        this->run.cycle.eval.log_interval = config::ea_2_log_eval_interval;
        
        // ea objective logging
        
        this->topologies.run.stats_out = config::ea_2_run_out;
        this->topologies.run.cycle.stats_out = config::ea_2_cycle_out;
        this->topologies.run.cycle.eval.stats_out = config::ea_2_eval_out;
        
        this->topologies.run.cycle.eval.log_stdout = true;
        this->topologies.run.log_interval = config::ea_2_o1_log_run_interval;
        this->topologies.run.cycle.log_interval = config::ea_2_o1_log_cycle_interval;
        this->topologies.run.cycle.eval.log_interval = config::ea_2_o1_log_eval_interval;
        
        // log population every nth cycle
        
        this->topologies.run.cycle.log_population_interval = config::ea_2_o1_log_population_interval;
        
        // log current genome every nth eval
        
        this->topologies.run.cycle.eval.log_genome_interval = config::ea_2_o1_log_genome_interval;
        
        // evolution properties
        
        this->topologies.mu = config::ea_2_mu;
        this->topologies.lambda = config::ea_2_lambda;
        this->topologies.mutation_rate = config::ea_2_mutation_rate;
        
        // calculate objective initialization time
        
        this->stats.init_duration += (MPI_Wtime() - this->stats.start);
        
    }
    
    template<typename e> void begin(e &target);
    template<typename i> void log_population(i &interval);
    
};

template<> std::vector<std::vector<topology>> objective<topology>::dominated_sort() {
    
    std::vector<std::vector<topology>> fronts(1);
    
    if(mpi.id != 0) { return fronts; }
    
    for(std::vector<topology>::iterator t1 = this->population.begin(); t1 != this->population.end(); ++t1) {
    
        long offset = std::distance(this->population.begin(), t1) + 1;
        
        for(std::vector<topology>::iterator t2 = this->population.begin() + offset; t2 != this->population.end(); ++t2) {
        
            if (t1 == t2) { continue; }
            
            if (t1->dominates(*t2)) {
                
                t1->dom_genomes.push_back(*t2);
                
                LOG(2, mpi.id, 0, "%d (%f,%f) dominates %d (%f,%f) and has %lu dom genomes\r\n", t1->id, t1->fitness_multi.first, t1->fitness_multi.second, t2->id, t2->fitness_multi.first, t2->fitness_multi.second, t1->dom_genomes.size());
                
            } else if(t2->dominates(*t1)) {
                
                t1->dom_count++;
                
                LOG(2, mpi.id, 0, "%d (%f,%f) dominates %d (%f,%f), dominated %d times\r\n", t2->id, t2->fitness_multi.first, t2->fitness_multi.second, t1->id, t1->fitness_multi.first, t1->fitness_multi.second, t1->dom_count);
                
            }
            
        }
        
        if(t1->dom_count == 0) {
            fronts[0].push_back(*t1);
            LOG(2, mpi.id, 0, "added %d (dominated by %d, dominates %lu) to fronts (size %lu)\r\n ", t1->id, t1->dom_count, t1->dom_genomes.size(), fronts[0].size());
        }

    }
    
    int idx = 0;
    
    for(std::vector<topology>::iterator t1 = fronts[idx].begin(); t1 != fronts[idx].end(); ++t1) {
        
        std::vector<topology> front;
        
        LOG(2, mpi.id, 0, "creating front %d from %lu genomes dominated by genome %d\r\n", idx, t1->dom_genomes.size(), t1->id);
        
        std::copy_if(t1->dom_genomes.begin(), t1->dom_genomes.end(), std::back_inserter(front), [&](topology &t) {
            t.dom_count -= 1;
            LOG(2, mpi.id, 0, "front %d rank %d genome %d dom_count = %d\r\n", idx, t1->id, t.id, t.dom_count);
            return t.dom_count == 0;
        });
        
       // for(std::vector<topology>::iterator t2 = front.dominated.begin(); t2 != front.dominated.end; ++t2)
        
    }

    return fronts;
    
}

template<> template<typename v> void objective<topology>::gather(v &variant) {
    
}

template<> template<typename i> void objective<topology>::log_population(i &interval) {

    if(mpi.id==0) {
        
        for (auto it = this->population.begin(); it != this->population.end(); ++it) {
            
            std::fprintf(config::ea_2_population_out, "%d," "%d," "%d," "%d," "%f," "%d," "%d," "%d," "%d," "%d," "%d," "%d," "%f\r\n", this->run.id, this->run.cycle.id, this->run.cycle.eval.id, it->id, it->fitness, it->stats.send_channels, it->stats.recv_channels, it->stats.total_channels, it->stats.arrivals, it->stats.departures, it->stats.migrations, it->stats.target_runs, it->selection_distribution);
            
        }
        
        fflush(config::ea_2_population_out);
        
    }

}

#endif /* dtype_meta_h */

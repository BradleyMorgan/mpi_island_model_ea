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

bool sort_o1(const std::pair<int,int> &a, const std::pair<int,int> &b) {
       return a.first > b.first;
}

bool sort_o2(const std::pair<int,int> &a, const std::pair<int,int> &b) {
       return a.second > b.second;
}

bool cmp_o1(const topology *lt, const topology *rt) {
  return lt->fitness_multi.first < rt->fitness_multi.first;
}

bool cmp_o2(const topology *lt, const topology *rt) {
    return lt->fitness_multi.second < rt->fitness_multi.second;
}

bool compare_multi(const topology *lt, const topology *rt) {
    return lt->dom_rank <= rt->dom_rank && (lt->distance > rt->distance);
}

struct PlotRange {
    double begin, end, count;
    double get_step() const { return (end - begin) / count; }
};

template <typename F>
void TextPlot2d(F func, PlotRange range_x, PlotRange range_y) {
    const auto step_x = range_x.get_step();
    const auto step_y = range_y.get_step();
    // multiply steps by iterated integer
    // to avoid accumulation of error in x and y
    for(int j = 0;; ++j) {
        auto y = range_y.begin + step_y * j;
        if(y >= range_y.end) break;
        for(int i = 0;; ++i) {
            auto x = range_x.begin + step_x * i;
            if(x >= range_x.end) break;
            auto z = func(x, y);
            if(z != z) { printf("?"); } // NaN outputs a '?'
            else       { printf("%c", z < 0 ? '#':'o'); }
        }
        printf("\n");
    }
}

void plot_fronts(std::vector<std::vector<topology*>> &fronts) {

    for(std::vector<std::vector<topology*>>::iterator it = fronts.begin(); it != fronts.end(); ++it) {
        
        std::pair<std::vector<topology*>::iterator, std::vector<topology*>::iterator> o1_minmax = std::minmax_element(it->begin(), it->end(), cmp_o1);
        
        std::pair<std::vector<topology*>::iterator, std::vector<topology*>::iterator> o2_minmax = std::minmax_element(it->begin(), it->end(), cmp_o2);
        
        //double o1_min = std::min_element(it->begin(), it->end(), cmp_o1);
       // double o1_max = std::max_element(it->begin(), it->end(), cmp_o1);
        
        for(std::vector<topology*>::iterator t = it->begin(); t != it->end(); ++t) {
        
            TextPlot2d(
                [](double x, double y){ return std::sin(x) + std::cos(y/2)*std::cos(y/2) - x/y; },
                {0.0, (*t)->fitness_multi.first, o1_minmax.second[0]->fitness_multi.first},
                {(*t)->fitness_multi.second, (*t)->fitness_multi.second, o1_minmax.second[0]->fitness_multi.second}
            );
            
        }
        
    }

}

double calculate_distance_o1(std::vector<topology*> &subjects, unsigned idx) {
    
    return subjects[idx]->distance + (subjects[idx+1]->fitness_multi.first - subjects[idx-1]->fitness_multi.first);
    
}
        
double calculate_distance_o2(std::vector<topology*> &subjects, unsigned idx) {
    
    return subjects[idx]->distance + (subjects[idx+1]->fitness_multi.second - subjects[idx-1]->fitness_multi.second);
    
}

// take the average distance of the two points on either side of this point along each of the objectives

void crowding_distance(std::vector<topology*> &front) {
    
    if(front.size() == 0) { return; }

    const double infinity = std::numeric_limits<double>::infinity();

    //for(unsigned m = 0; m < front.size(); ++m) {
        
        std::sort(front.begin(), front.end(), [](const topology *a, const topology *b) {
            return (*a).fitness_multi.first < (*b).fitness_multi.first;
        });
        
        front[0]->distance = infinity;
        front[front.size()-1]->distance = infinity;

        for(unsigned k = 1; k < front.size()-1; ++k) {
            
            front[k]->distance = calculate_distance_o1(front, k);
        
        }
    
    std::sort(front.begin(), front.end(), [](const topology *a, const topology *b) {
        return (*a).fitness_multi.second < (*b).fitness_multi.second;
    });
    
    front[0]->distance = infinity;
    front[front.size()-1]->distance = infinity;

    for(unsigned k = 1; k < front.size()-1; ++k) {
        
        front[k]->distance = calculate_distance_o2(front, k);
    
    }

    printf("\r\n");

}

template<> std::vector<std::vector<topology*>> objective<topology>::define_fronts() {
    
    // create an array instance to store pareto fronts
    // return an empty front if the parallel ea process is not root
    
    std::vector<std::vector<topology*>> fronts;
    std::vector<topology*> front_c;
    
    if(mpi.id != 0) { return fronts; }
    
    // for each solution calculate:
    // the number of times the solution @var{population[i]} is dominated
    // set of individuals @var{population[j..n]} dominated by the solution at @var{population[i]}
        
    for(std::vector<topology>::iterator t1 = this->population.begin(); t1 != this->population.end(); ++t1) {
        
        // @var{population[i]}
        
        for(std::vector<topology>::iterator t2 = this->population.begin(); t2 != this->population.end(); ++t2) {
        
            // @var{population[j]}
            
            if (t1 == t2) { continue; }
            
            if (t1->dominates(*t2)) {
                
                // add population[j] to set of @var{population[i]} dominated indivduals
                
                t1->dom_genomes.push_back(&*t2);
                
                LOG(2, mpi.id, 0, "%d (%f,%f) dominates %d (%f,%f) and has %lu dom genomes\r\n", t1->id, t1->fitness_multi.first, t1->fitness_multi.second, t2->id, t2->fitness_multi.first, t2->fitness_multi.second, t1->dom_genomes.size());
                
            } else if(t2->dominates(*t1)) {
                
                // @var{population[j]} is dominant, so increment @var{population[i]} domination count
                
                t1->dom_count++;
                
                LOG(2, mpi.id, 0, "%d (%f,%f) dominates %d (%f,%f), dominated %d times\r\n", t2->id, t2->fitness_multi.first, t2->fitness_multi.second, t1->id, t1->fitness_multi.first, t1->fitness_multi.second, t1->dom_count);
                
            }
            
        }
        
        // @var{population[i]} domination comparison complete
        // if it was not dominated, add it to @var{front_c}
        
        if(t1->dom_count == 0) {
            
            t1->dom_rank = 1;
            front_c.push_back(&*t1);
            
        }

    }
    
    crowding_distance(front_c);
    fronts.push_back(front_c);
    
    // for each solution in @var{front_c} we iterate over each solution @var{i}
    
    int idx = 0;
    
    while(fronts[idx].size() != 0) {
    
        std::vector<topology*> front_n;
    
        for(std::vector<topology*>::iterator i = fronts[idx].begin(); i != fronts[idx].end(); ++i) {
        
            // decrement the associated @var{dom_count} by one
            // if @var{dom_count} == 0 for any solution in the current front @var{front_c}
            // add it to set @var{front_n}
            
            std::copy_if((*i)->dom_genomes.begin(), (*i)->dom_genomes.end(), std::back_inserter(front_n), [&](topology *t) {
                (*t).dom_count -= 1;
                if((*t).dom_count == 0) {
                    (*t).dom_rank = idx + 1;
                }
                return(*t).dom_count == 0;
            });
        
        }
        
        idx++;

        crowding_distance(front_n);
       
        fronts.push_back(front_n);
        
    }
    
    return fronts;
    
}

template<> template<typename v> void objective<topology>::gather(v &variant) {
    
    // placeholder
    
}

template<> template<typename i> void objective<topology>::log_population(i &interval) {

    if(mpi.id==0) {
        
        for (auto it = this->population.begin(); it != this->population.end(); ++it) {
            
            std::fprintf(config::ea_2_population_out, "%d," "%d," "%d," "%d," "%f," "%f," "%d," "%f," "%d," "%lu," "%d," "%d," "%d," "%d," "%d," "%d," "%d," "%f\r\n", this->run.id, this->run.cycle.id, this->run.cycle.eval.id, it->id, it->fitness_multi.first, it->fitness_multi.second, it->dom_rank, it->distance, it->dom_count, it->dom_genomes.size(), it->stats.send_channels, it->stats.recv_channels, it->stats.total_channels, it->stats.arrivals, it->stats.departures, it->stats.migrations, it->stats.target_runs, it->selection_distribution);
            
        }
        
        fflush(config::ea_2_population_out);
        
    }

}

#endif /* dtype_meta_h */

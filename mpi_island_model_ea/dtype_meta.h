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
    void log_fronts(std::vector<std::vector<topology*>> &fronts);
    
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
    if (lt->dom_rank < rt->dom_rank) { return true; }
    return lt->distance > rt->distance;
}

void ea_meta::log_fronts(std::vector<std::vector<topology*>> &fronts) {
    
    if(mpi.id==0) {
        
        for (std::vector<std::vector<topology*>>::iterator front = fronts.begin(); front != fronts.end(); ++front) {
            
            for (std::vector<topology*>::iterator t = front->begin(); t != front->end(); ++t) {
                
                double d = (*t)->distance == std::numeric_limits<double>::infinity() ? INFINITY : (*t)->distance;
                
                std::fprintf(config::ea_2_multi_out, "%d," "%d," "%d," "%d," "%f," "%d," "%lu," "%f," "%f\r\n", this->topologies.run.id, this->topologies.run.cycle.id, this->topologies.run.cycle.eval.id, (*t)->dom_rank, d, (*t)->dom_count, (*t)->dom_genomes.size(), (*t)->fitness_multi.first, (*t)->fitness_multi.second);
                
            }
            
        }
        
        fflush(config::ea_2_multi_out);
        
    }
    
    
}

void plot_fronts(std::vector<std::vector<topology*>> &fronts) {

    std::vector<std::vector<int>> plane;
    
    plane.resize(fronts.size());
    
    for(std::vector<std::vector<topology*>>::iterator front = fronts.begin(); front != fronts.end(); ++front) {
        
 //       std::vector<int> y;
        
      // unsigned long y_len = front->size() - 1;
//
//        y.resize(y_len);
        
        std::pair<std::vector<topology*>::iterator, std::vector<topology*>::iterator> o1_minmax = std::minmax_element(front->begin(), front->end(), cmp_o1);
        
        double range_min = (*o1_minmax.first)->fitness_multi.first;
        double range_max = (*o1_minmax.second)->fitness_multi.first;
        double range_step = range_max - range_min / front->size();
        //double range_offset = range_step * 1.5;
        
        std::vector<int> y;
        
        for(std::vector<topology*>::iterator t = front->begin(); t != front->end(); ++t) {
  
        }
        

        // std::pair<std::vector<topology*>::iterator, std::vector<topology*>::iterator> o2_minmax =
        // std::minmax_element(front->begin(), front->end(), cmp_o2);
        // double o1_min = std::min_element(it->begin(), it->end(), cmp_o1);
        // double o1_max = std::max_element(it->begin(), it->end(), cmp_o1);
        
    }

}

double calculate_distance_o1(std::vector<topology*> &subjects, unsigned idx) {
    
    double result = subjects[idx]->distance + abs(subjects[idx+1]->fitness_multi.first - subjects[idx-1]->fitness_multi.first) / abs(subjects[subjects.size()-1]->fitness_multi.first - subjects[0]->fitness_multi.first);
    
    LOG(2, 0, 0, "front index %d o1 distance = %f\r\n", idx, result);
    
    return result;
    
}
        
double calculate_distance_o2(std::vector<topology*> &subjects, unsigned idx) {
    
    double result = subjects[idx]->distance + abs(subjects[idx+1]->fitness_multi.second - subjects[idx-1]->fitness_multi.second) / abs(subjects[subjects.size()-1]->fitness_multi.second - subjects[0]->fitness_multi.second);;
    
    LOG(2, 0, 0, "front index %d o2 distance = %f\r\n", idx, result);
    
    return result;
    
}

// take the average distance of the two points on either side of this point along each of the objectives

template<> void objective<topology>::crowding_distance(std::vector<std::vector<topology*>> &fronts) {
    
    if(mpi.id != 0 || fronts[0].size() == 0) { return; }

    for(std::vector<std::vector<topology*>>::iterator front = fronts.begin(); front != fronts.end(); ++front) {
        
        if(front->size() == 0) { return; }
        
        const double infinity = std::numeric_limits<double>::infinity();

        // sort front by objective 1 fitness
        // currently maximization of solution quality
            
        std::sort(front->begin(), front->end(), [](const topology *a, const topology *b) {
            return (*a).fitness_multi.first > (*b).fitness_multi.first;
        });
        
        // set o1 relative best and worst boundaries
        
        (*front)[0]->distance = infinity;
        (*front)[front->size()-1]->distance = infinity;

        // for each individual in the current front, measure the distance of o1
        // from its neighboring solutions
        // for maximization, sorted front will contain values in ascending order
        // so idx-1 will be the next worst and idx+1 will be the next best
        // presumably idx+1 > idx-1
        
        // D_i = F_best[+1] (o1 fitness of next worst) - F_best[-1] (o1 fitness of next best)
        
        for(unsigned k = 1; k < front->size()-1; ++k) {
            
            (*front)[k]->distance = calculate_distance_o1(*front, k);
        
        }

        std::sort(front->begin(), front->end(), [](const topology *a, const topology *b) {
            return (*a).fitness_multi.second > (*b).fitness_multi.second;
        });

        for(unsigned k = 1; k < front->size()-1; ++k) {
            
            (*front)[k]->distance = calculate_distance_o2(*front, k);
        
        }
        
    }

    printf("\r\n");

}

template<> std::vector<std::vector<topology*>> objective<topology>::define_fronts() {
    
    // create an array instance to store pareto fronts
    // return an empty front if the parallel ea process is not root
    
    std::vector<std::vector<topology*>> fronts;
    std::vector<topology*> front_c;
    
    unsigned long count = 0;
    
    if(mpi.id != 0) { return fronts; }
    
    // for each solution calculate:
    // the number of times the solution @var{population[i]} is dominated
    // set of individuals @var{population[j..n]} dominated by the solution at @var{population[i]}
    
    for(std::vector<topology>::iterator t1 = this->population.begin(); t1 != this->population.end(); ++t1) {
        t1->init_multi();
    }
    
    for(std::vector<topology>::iterator t1 = this->population.begin(); t1 != this->population.end(); ++t1) {
        
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
            
            LOG(2, 0, 0, "adding topology %d to initial front (size %lu) ...\r\n", t1->id, front_c.size());
            
            t1->dom_rank = 1;
            front_c.push_back(&*t1);
            
            LOG(2, 0, 0, "added topology %d, to first front (size %lu) ...\r\n", t1->id, front_c.size());
            
        }

    }
        
    LOG(2, 0, 0, "adding first front (size %lu) to fronts (size %lu) ...\r\n", front_c.size(), fronts.size());
    
    fronts.push_back(front_c);
    count += front_c.size();
    
    LOG(2, 0, 0, "added front %d (size %lu) to fronts (size %lu) ...\r\n", 0, front_c.size(), fronts.size());
    
    // for each solution in @var{front_c} we iterate over each solution @var{i}
    
    int idx = 0;
    
    LOG(2, 0, 0, "iterating front genomes ...\r\n");
    
    while(fronts[idx].size() != 0) {
    
        std::vector<topology*> front_n;
    
        LOG(2, 0, 0, "creating front %d  ...\r\n", idx);
        
        for(std::vector<topology*>::iterator i = fronts[idx].begin(); i != fronts[idx].end(); ++i) {
        
            // decrement the associated @var{dom_count} by one
            // if @var{dom_count} == 0 for any solution in the current front @var{front_c}
            // add it to set @var{front_n}
            
            std::copy_if((*i)->dom_genomes.begin(), (*i)->dom_genomes.end(), std::back_inserter(front_n), [&](topology *t) {
                LOG(2, 0, 0, "checking topology %d dom_count %d in front %d (size %lu)\r\n", (*t).id, (*t).dom_count, idx, front_n.size());
                (*t).dom_count -= 1;
                if((*t).dom_count == 0) {
                    (*t).dom_rank = idx + 1;
                    LOG(2, 0, 0, "assigned rank %d to topology %d dom_count %d in front %d (size %lu)\r\n", (*t).dom_rank, (*t).id, (*t).dom_count, idx, front_n.size());
                }
                return(*t).dom_count == 0;
            });
            
        }
        
        LOG(2, 0, 0, "created front %d (size %lu)\r\n", idx, front_n.size());
        
        idx++;
       
        LOG(2, 0, 0, "adding front %d (size %lu) to fronts (size %lu) ...\r\n", idx, front_n.size(), fronts.size());
        
        fronts.push_back(front_n);
        count += front_n.size();
        
        LOG(2, 0, 0, "added front %d (size %lu) to fronts (size %lu) ...\r\n", idx, front_n.size(), fronts.size());
        
    }
    
    
    LOG(2, 0, 0, "returning fronts (size %lu) with %lu individuals ...\r\n", fronts.size(), count);
    
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

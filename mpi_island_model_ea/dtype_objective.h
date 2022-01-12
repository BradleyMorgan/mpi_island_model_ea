//
//  objective.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 8/18/21.
//  Copyright © 2021 Bradley Morgan. All rights reserved.
//

#ifndef objective_h
#define objective_h

# pragma mark EA::OBJECTIVE::DATATYPE: @objective{}

template<typename genome> struct objective {
    
    char name[128];
    
    int id = 0;
    int mu = 0;
    int lambda = 0;
    
    FILE *population_out;
    FILE *genome_out;
    
    double mutation_rate = 0.0;

    std::vector<genome> population = {};
    
    #pragma mark DATATYPE: @metrics{}
    
    struct evaluation_interval;
    struct cycle_interval;
    struct run_interval;
    
    run_interval run;
    
    genome *local = run.cycle.eval.local;
    
    int genome_max_evals = 0;
    
    struct metrics {
        
        struct values {
            
            std::vector<double> cpd;
            double fitness = 0.0;
            
            values() : fitness(0.0) {}
            
        };
        
        metrics(void) {}
        
        values value;
        
    };
    
    #pragma mark EA::FUNCTION::TEMPLATES: metrics
    
    metrics aggregate;
    
    void cpd();
    void fitness();
    
    void init() {
        this->population.clear();
        this->local = {};
        this->aggregate = {};
        this->run.init();
    }
    
    std::vector<genome> filter(bool check_evals);
    
    #pragma mark EA::FUNCTION::TEMPLATES: initialization
    
    objective<genome>(int n) : id(n), local(run.cycle.eval.local), run(&population[0]) {};
    
    template<typename i> void end(i &interval, genome *current = NULL);
    template<typename i> void begin(i &interval, genome *current = NULL);
    
    #pragma mark EA::FUNCTION::TEMPLATES: initialization, logging

    // interval stats, mostly for timing metrics
    
    template<typename i> void log_begin(i &interval);
    template<typename i> void log_end(i &interval);
    template<typename i> void log_stats(i &interval);
    
    // the following functions handle inter-ea dependencies, e.g. meta ea
    // TODO: look for a cleaner way to do this

    template<typename source, typename target> void begin(target &interval, source &data_interval, genome *target_genome);
    template<typename source, typename target> void end(target &interval, source &data_interval, genome *target_genome);
    template<typename source, typename target> void log_begin(target &interval, source &data_interval, genome *target_genome);
    template<typename source, typename target> void log_end(target &interval, source &data_interval, genome *target_genome);
    template<typename sint, typename tint> void measure(tint &target, sint &source, genome *current);
    
    // various functions for calculating objective property values
    
    template<typename i> void measure(i &interval, genome *g);
    
    template<typename i> void log_population(i &interval);
    template<typename i> void log_genome(i &interval);
    
    template<typename i> double reduce(i &interval, mpi_local *locf, mpi_local *minf, mpi_local *maxf, mpi_local *sumf);
    
    void crowding_distance(std::vector<std::vector<topology*>> &fronts);
    std::vector<std::vector<topology*>> define_fronts();
    
    void minmax(mpi_local *field, mpi_local result, double const &(*func)(double const&, double const&));
    
    #pragma mark EA::FUNCTION::TEMPLATES: evolution cycle

    template<typename f> genome select(f function) { return function(*this); }
    template<typename f> std::pair<genome,genome> tournament(f function) { return function(*this); }
    
    template<typename f, typename v> genome crossover(v &variant, f function) { return function(*this, variant); }
    template<typename f, typename v> void populate(v &variant, f function) { function(variant); }
    template<typename f, typename v> void distribute(v &variant, f function) { function(variant); }
    template<typename v> void gather(v &variant);
    template<typename f, typename v, typename m> void evolve(v &variant, m &meta, f function) { function(variant, meta); }
    template<typename f, typename v, typename m, typename g> void evolve(v &variant, m &meta, g &individual, f function) { function(variant, meta, individual); }
    
};

template<typename genome> struct obj_eval_stats;
template<typename genome> struct obj_cycle_stats;
template<typename genome> struct obj_run_stats;

template<typename genome> struct objective<genome>::evaluation_interval {
    
    char name[5];
    
    int id = 0;
    int max = 0;
    int max_local = 0;
    
    int log_interval = 0;
    int log_population_interval = 0;
    int log_genome_interval = 0;
    
    bool log_stdout = false;
    bool log_fout = true;
    
    char log_head[160];
    char log_tail[160];
    
    FILE *stats_out;
    
    obj_eval_stats<genome> stats;
    
    cycle_interval *parent;
    
    genome *local;
    
    mpi_local stat_send;
    mpi_local stat_recv;
    
    void init() {
        this->local = {};
        this->stat_send = {};
        this->stat_recv = {};
    }
    
    void begin(genome *current) {
        this->id++;
        this->local = current;
        this->stats.start = MPI_Wtime();
        this->stats.begin_header(this->log_head);
    }

    void end(genome *current = NULL) {
        this->stats.local_t.value = MPI_Wtime() - this->stats.start;
        this->stats.end_header(this->log_tail);
    }
    
    evaluation_interval(cycle_interval *cycle): parent(cycle), stats() { strcpy(name, "EVAL"); };

};

template<typename genome> struct objective<genome>::cycle_interval {

    char name[6];
    
    int id = 1;
    int max = 0;
    
    const int max_local = 0;
    
    int log_interval = 0;
    int log_population_interval = 0;
    int log_genome_interval = 0;
    
    bool log_stdout = true;
    bool log_fout = true;
    
    char log_head[160];
    char log_tail[160];
    
    FILE *stats_out;
    
    mpi_local stat_send;
    mpi_local stat_recv;
    
    obj_cycle_stats<genome> stats;
    
    evaluation_interval eval;
    run_interval *parent;
    
    genome *local = eval.local;
    
    void init() {
        this->local = {};
        this->stats = {};
        this->eval.stats = {};
        this->stat_send = {};
        this->stat_recv = {};
        this->eval.init();
    }
    
    void begin(genome *current) {
        this->local = current;
        this->eval.stats = {};
        this->stats.start = MPI_Wtime();
        this->stats.begin_header(this->log_head);
    }

    void end(genome *current = NULL) {
        this->stats.local_t.value = MPI_Wtime() - this->stats.start;
        this->stats.end_header(this->log_tail);
    }
    
    cycle_interval(run_interval *run, genome *current): local(current), stats(), eval(this) { strcpy(name, "CYCLE"); };
    
};

template<typename genome> struct objective<genome>::run_interval {
    
    char name[4];
    
    int id = 1;
    int max = 0;
    
    const int max_local = 0;
    
    int log_interval = 0;
    int log_population_interval = 0;
    int log_genome_interval = 0;
    
    bool log_stdout = true;
    bool log_fout = true;

    char log_head[160];
    char log_tail[160];
    
    FILE *stats_out;
    
    obj_run_stats<genome> stats;
    
    cycle_interval cycle;
    run_interval *parent;
    
    genome *local;
    
    mpi_local stat_send;
    mpi_local stat_recv;
    
    void init() {
        this->local = {};
        this->stats = {};
        this->cycle.stats = {};
        this->stat_send = {};
        this->stat_recv = {};
        this->cycle.init();
    }
    
    void begin(genome *current) {
        this->local = current;
        this->stats = {};
        this->cycle.stats = {};
        this->stats.start = MPI_Wtime();
        this->stats.begin_header(this->log_head);
    }

    void end(genome *current = NULL) {
        this->stats.local_t.value = MPI_Wtime() - this->stats.start;
        this->stats.end_header(this->log_tail);
    }
    
    run_interval(genome *current): parent(this), local(current), stats(), cycle(this, current) { strcpy(name, "RUN"); };
    
};


#pragma mark EA::FUNCTION::TEMPLATE: <*>fitness{}

// calculates an arbitrary population's aggregate fitness

template<typename genome> void objective<genome>::fitness() {
    
     LOG(8, 0, 0, "calculating total fitness ...\r\n");
    
    this->aggregate.value.fitness = 0.0;
    
    for(typename std::vector<genome>::iterator it = this->population.begin(); it != this->population.end(); ++it) {
        LOG(10, 0, 0, "genome %d fitness = %f\r\n", this->id, it->fitness);
        this->aggregate.value.fitness += it->fitness;
    }
    
    LOG(8, 0, 0, "objective total population fitness = %f...\r\n", this->aggregate.value.fitness);
    
}

#pragma mark EA::FUNCTION::TEMPLATE: <*>cpd{}

// calculates an arbitrary population's cumulative probability distribution

template<typename genome> void objective<genome>::cpd() {
    
    LOG(8, 0, 0, "generating cpd\r\n");

    double cumulative_probability = 0.0;

    this->aggregate.value.cpd.clear();
    
    this->fitness();

    LOG(8, 0, 0, "sorting population descending fitness ...\r\n");

    std::sort(this->population.begin(), this->population.end(), compare_fitness<genome>);
    std::reverse(this->population.begin(), this->population.end());

    LOG(3, 0, 0, "objective %d cpd calculation total fitness = %f\r\n", this->id, this->aggregate.value.fitness);

    for(int i=0; i < this->population.size(); i++) {

        this->population[i].selection_distribution = this->population[i].fitness / this->aggregate.value.fitness;

        LOG(8, 0, 0, "calculated island %d solution %d fitness %f selection distribution = %f\r\n", this->id, i, this->population[i].fitness, this->population[i].selection_distribution);

        cumulative_probability += this->population[i].selection_distribution;

        LOG(8, 0, 0, "%d=%f|", i, cumulative_probability);

        this->aggregate.value.cpd.push_back(cumulative_probability);;

    }
    
    LOG(3, 0, 0, "(%d,%d,%d,%d,%d) end cpd() sum_fit = %f\r\n", mpi.id, this->id, this->run.id, this->run.cycle.id, this->run.cycle.eval.id, this->aggregate.value.fitness);

}

#pragma mark EA::INTERVAL::FUNCTION filter{}

// return only fully evaluated population members for intervals with
// high frequencies

template<typename genome> std::vector<genome> objective<genome>::filter(bool check_evals) {
    
    std::vector<genome> filtered;
    
    int max_rounds = this->run.max * this->run.cycle.max;
    
    if(this->population.size() == 0) { return filtered; }
    
    if(check_evals) {
        std::copy_if( this->population.begin(), this->population.end(), std::back_inserter(filtered), [&](const genome &item) { return item.fitness != 0.0 && this->run.cycle.eval.id >= max_rounds; });
    } else {
        std::copy_if( this->population.begin(), this->population.end(), std::back_inserter(filtered), [&](const genome &item) { return item.fitness != 0.0; });
    }
    
    std::sort(filtered.begin(), filtered.end(), compare_fitness<genome>);
    std::reverse(filtered.begin(), filtered.end());
    
    if(filtered.size() == 0) { filtered.push_back(this->population[0]); }
    
    return filtered;
    
}

#pragma mark EA::PARALLEL::INTERVAL reduce{}

// aggregate min|max|sum <field> statistics from all islands and return the average

template<typename genome> template<typename i> double objective<genome>::reduce(i &interval, mpi_local *locf, mpi_local *minf, mpi_local *maxf, mpi_local *sumf) {
    
    MPI_Reduce(locf, &interval.stat_recv, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);
    
    if(mpi.id == 0) { this->minmax(minf, interval.stat_recv, std::min<double>); }
    
    MPI_Reduce(locf, &interval.stat_recv, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);
    
    if(mpi.id == 0) { this->minmax(maxf, interval.stat_recv, std::max<double>); }
    
    MPI_Reduce(locf, sumf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    double avg = 0.0;
    
    if(mpi.id == 0) { avg = sumf->value / mpi.size; }

    return avg;
}

#pragma mark EA::META::OBJECTIVE::FUNCTION: <>measure()

// accept a target interval type to receive metrics from a source interval
// 
template<typename genome> template<typename sint, typename tint> void objective<genome>::measure(tint &target, sint &source, genome *current) {
    
    // aggregate time values from all islands

    target.stats.avg_t = this->reduce(target, &source.stats.local_t, &target.stats.min_t, &target.stats.max_t, &target.stats.sum_t);
    target.stats.avg_scatter_t = this->reduce(target, &source.stats.local_scatter_t, &target.stats.min_scatter_t, &target.stats.max_scatter_t, &target.stats.sum_scatter_t);
    target.stats.avg_gather_t = this->reduce(target, &source.stats.local_gather_t, &target.stats.min_gather_t, &target.stats.max_gather_t, &target.stats.sum_gather_t);
    target.stats.avg_migration_t = this->reduce(target, &source.stats.local_migration_t, &target.stats.min_migration_t, &target.stats.max_migration_t, &target.stats.sum_migration_t);
    
    if(mpi.id != 0) { return; }
    
    current->measure(target);
    
    std::vector<genome> valid = this->filter(false);

    if(valid.size() == 0) { return; }
    
    if(target.stats.best_fitness == 0.0 || valid[0].fitness > target.stats.best_fitness) {
        target.stats.best_fitness = valid[0].fitness;
        target.stats.best.push_back(valid[0].fitness);
    }
    
    if(this->run.cycle.stats.best_fitness == 0.0 || target.stats.best_fitness > this->run.cycle.stats.best_fitness) {
        this->run.cycle.stats.best_fitness = target.stats.best_fitness;
        this->run.cycle.stats.best.push_back(target.stats.best_fitness);
    }
    
    if(this->run.stats.best_fitness == 0.0 || this->run.cycle.stats.best_fitness > this->run.stats.best_fitness) {
        this->run.stats.best_fitness = this->run.cycle.stats.best_fitness;
        this->run.stats.best.push_back(this->run.cycle.stats.best_fitness);
    }
    
    auto result = minmax_element(valid.begin(), valid.end(),[](genome const &g1, genome const &g2) { return g1.fitness < g2.fitness;});
    
    target.stats.min_fitness = result.first->fitness;
    target.stats.max_fitness = result.second->fitness;
    target.stats.total_fitness = std::accumulate(valid.begin(), valid.end(), 0.0, [](double v, const genome& o){ return o.fitness + v; });
    target.stats.total_o1_fitness = std::accumulate(valid.begin(), valid.end(), 0.0, [](double v, const genome& o){ return o.fitness_multi.first + v; });
    target.stats.total_o2_fitness = std::accumulate(valid.begin(), valid.end(), 0.0, [](double v, const genome& o){ return o.fitness_multi.first + v; });
    target.stats.total_best_fitness = std::accumulate(target.stats.best.begin(), target.stats.best.end(), 0.0);
    target.stats.avg_fitness = target.stats.total_fitness / valid.size();
    target.stats.avg_o1_fitness = target.stats.total_o1_fitness / valid.size();
    target.stats.avg_o2_fitness = target.stats.total_o2_fitness / valid.size();
    
    target.stats.avg_best_fitness = target.stats.total_best_fitness / target.stats.best.size();
    
    current->fitness_multi.second = source.stats.avg_best_fitness;
    
    printf("");

}

template<typename genome> template<typename i> void objective<genome>::measure(i &interval, genome *current) {
     
    interval.stats.avg_t = this->reduce(interval, &interval.stats.local_t, &interval.stats.min_t, &interval.stats.max_t, &interval.stats.sum_t);
    interval.stats.avg_scatter_t = this->reduce(interval, &interval.stats.local_scatter_t, &interval.stats.min_scatter_t, &interval.stats.max_scatter_t, &interval.stats.sum_scatter_t);
    interval.stats.avg_gather_t = this->reduce(interval, &interval.stats.local_gather_t, &interval.stats.min_gather_t, &interval.stats.max_gather_t, &interval.stats.sum_gather_t);
    interval.stats.avg_migration_t = this->reduce(interval, &interval.stats.local_migration_t, &interval.stats.min_migration_t, &interval.stats.max_migration_t, &interval.stats.sum_migration_t);
    
    if(mpi.id != 0) { return; }
    
    if(current) {
        current->measure(interval);
    }
    
    std::vector<genome> valid = this->filter(false);

    if(valid.size() == 0) { return; }

    if(interval.stats.best_fitness == 0.0 || valid[0].fitness > interval.stats.best_fitness) {
        interval.stats.best_fitness = valid[0].fitness;
        interval.stats.best.push_back(valid[0].fitness);
    }
    
    if(this->run.cycle.stats.best_fitness == 0.0 || interval.stats.best_fitness > this->run.cycle.stats.best_fitness) {
        this->run.cycle.stats.best_fitness = interval.stats.best_fitness;
        this->run.cycle.stats.best.push_back(interval.stats.best_fitness);
    }
    
    if(this->run.stats.best_fitness == 0.0 || this->run.cycle.stats.best_fitness > this->run.stats.best_fitness) {
        this->run.stats.best_fitness = this->run.cycle.stats.best_fitness;
        this->run.stats.best.push_back(this->run.cycle.stats.best_fitness);
    }
    
    auto result = minmax_element(valid.begin(), valid.end(),[](genome const &g1, genome const &g2) { return g1.fitness < g2.fitness;});
    
    interval.stats.min_fitness = result.first->fitness;
    interval.stats.max_fitness = result.second->fitness;
    
    interval.stats.total_fitness = std::accumulate(valid.begin(), valid.end(), 0.0, [](double v, const genome& o){ return o.fitness + v; });
    double total_best_fitness = std::accumulate(interval.stats.best.begin(), interval.stats.best.end(), 0.0);
    
    interval.stats.avg_fitness = interval.stats.total_fitness / valid.size();
    interval.stats.avg_best_fitness = total_best_fitness / interval.stats.best.size();
    
    
}

#pragma mark EA::OBJECTIVE::INTERVAL statistics logging

template<typename genome> template<typename i> void objective<genome>::log_stats(i &interval) {
    
    if(interval.log_interval != 0 && interval.id % interval.log_interval == 0) {
        
        if(interval.log_stdout == true) {
    
            if(mpi.id != 0) { return; }
            
            LOG(2, 0, 0, "%04d,%06d,%08d,%11.6f,%11.6f,%11.6f,%11.6f,%11.6f,%11.6f,%11.6f(%d),%11.6f(%d),%11.6f\r\n",
                
                this->run.id,
                this->run.cycle.id,
                this->run.cycle.eval.id,
                interval.stats.avg_fitness,
                interval.stats.best_fitness,
                interval.stats.sum_scatter_t.value,
                interval.stats.sum_gather_t.value,
                interval.stats.sum_migration_t.value,
                interval.stats.sum_t.value,
                interval.stats.min_t.value,
                interval.stats.min_t.island,
                interval.stats.max_t.value,
                interval.stats.max_t.island,
                interval.stats.avg_t);
        
        }
        
        if(interval.log_fout == true) {

            if(mpi.id != 0) { return; }

            fprintf(interval.stats_out, "%d,%d,%d,%f,%f,%f,%f,%d,%f,%d,%f,%f,%f,%d,%f,%d,%f,%f,%f,%d,%f,%d,%f,%f,%f,%d,%f,%d,%f,%f,%f,%f\r\n",

                this->run.id,
                this->run.cycle.id,
                this->run.cycle.eval.id,
                interval.stats.avg_fitness,
                interval.stats.best_fitness,
                interval.stats.avg_best_fitness,
                interval.stats.min_scatter_t.value,
                interval.stats.min_scatter_t.island,
                interval.stats.max_scatter_t.value,
                interval.stats.max_scatter_t.island,
                interval.stats.sum_scatter_t.value,
                interval.stats.avg_scatter_t,
                interval.stats.min_gather_t.value,
                interval.stats.min_gather_t.island,
                interval.stats.max_gather_t.value,
                interval.stats.max_gather_t.island,
                interval.stats.sum_gather_t.value,
                interval.stats.avg_gather_t,
                interval.stats.min_migration_t.value,
                interval.stats.min_migration_t.island,
                interval.stats.max_migration_t.value,
                interval.stats.max_migration_t.island,
                interval.stats.sum_migration_t.value,
                interval.stats.avg_migration_t,
                interval.stats.min_t.value,
                interval.stats.min_t.island,
                interval.stats.max_t.value,
                interval.stats.max_t.island,
                interval.stats.sum_t.value,
                interval.stats.avg_t,
                interval.stats.start,
                interval.stats.local_t.value);

        }
        
        fflush(interval.stats_out);

    }
    
}

#pragma mark EA::OBJECTIVE::INTERVAL logging wrappers for standard (independent) ea

template<typename genome> template<typename i> void objective<genome>::log_begin(i &interval) {
    
    LOG(7, 0, 0, "\r\n--- (%d) BEGIN OBJECTIVE %s %d OF %d ---\r\n", mpi.id, interval.name, interval.id, interval.max);
    LOG(4, mpi.id, 0, "\r\n--- (%d) BEGIN OBJECTIVE %s %d OF %d ---\r\n", mpi.id, interval.name, interval.id, interval.max);
    
}

template<typename genome> template<typename i> void objective<genome>::log_end(i &interval) {
    
    LOG(7, 0, 0, "\r\n--- (%d) END OBJECTIVE %s %d OF %d ---\r\n", mpi.id, interval.name, interval.id, interval.max);
    LOG(4, mpi.id, 0, "\r\n--- (%d) END OBJECTIVE %s %d OF %d ---\r\n", mpi.id, interval.name, interval.id, interval.max);
    
    this->measure(interval, interval.local);
    
    if(interval.log_interval != 0 && interval.id % interval.log_interval == 0) {
        this->log_stats(interval);
    }
    
}

#pragma mark EA::OBJECTIVE::INTERVAL logging wrappers for dependent (meta) ea

template<typename genome> template<typename source, typename target> void objective<genome>::log_begin(target &interval, source &data_interval, genome *target_genome) {
    
}

template<typename genome> template<typename source, typename target> void objective<genome>::log_end(target &interval, source &data_interval, genome *target_genome) {
    
    this->measure(interval, data_interval, target_genome);
    
    if(interval.log_interval != 0 && interval.id % interval.log_interval == 0) {
        this->log_stats(interval);
    }
    
}

#pragma mark EA::OBJECTIVE::INTERVAL control wrappers for independent ea

template<typename genome> template<typename i> void objective<genome>::begin(i &interval, genome *current) {
    
    interval.begin(current);
    this->log_begin(interval);
    
}

template<typename genome> template<typename i> void objective<genome>::end(i &interval, genome *current) {
    
    interval.end(current);
    this->log_end(interval);
    
    if(interval.log_population_interval != 0 && interval.id%interval.log_population_interval == 0) {
        this->log_population(interval);
    }
    
    if(interval.log_genome_interval != 0 && interval.id%interval.log_genome_interval == 0) {
        current->log(interval);
    }
    
}

#pragma mark EA::OBJECTIVE::INTERVAL control wrappers for dependent (meta) ea

template<typename genome> template<typename source, typename target> void objective<genome>::begin(target &interval, source &data_interval, genome *target_genome) {

    interval.begin(target_genome);
    
    this->log_begin(interval, data_interval, target_genome);
    
}

template<typename genome> template<typename source, typename target> void objective<genome>::end(target &interval, source &data_interval, genome *target_genome) {
    
    interval.end(target_genome);
    
    this->log_end(interval, data_interval, target_genome);
    
    if(interval.log_population_interval != 0 && interval.id%interval.log_population_interval == 0) {
        this->log_population(interval);
    }
    
    if(interval.log_genome_interval != 0 && interval.id%interval.log_genome_interval == 0) {
        target_genome->log(interval);
    }
    
}

template<typename genome>
void objective<genome>::minmax(mpi_local *field, mpi_local result, double const &(*func)(double const&, double const&)) {
    *field = (*field).value == 0.0 ? result : func((*field).value, result.value) != (*field).value ? result : *field;
}

#endif /* objective_h */

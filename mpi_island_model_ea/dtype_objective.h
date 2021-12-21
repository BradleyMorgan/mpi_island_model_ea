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
    
    char name[128] = "";
    
    int id = 1;
    int mu = 0;
    int lambda = 0;
    int log_interval = 0;
    
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
    
    objective<genome>() : id(1), local(run.cycle.eval.local), run() {};
    
    template<typename i> void end(i &interval, genome *current = NULL);
    template<typename sint, typename tint> void end(tint &target, sint &source, genome *current);
    template<typename i> void begin(i &interval, genome *current = NULL);

    #pragma mark EA::FUNCTION::TEMPLATES: initialization logging
    
    template<typename g> std::vector<g> filter(bool check_evals);

    template<typename i> void log_begin(i &interval);
    template<typename i> void log_end(i &interval);
    template<typename sint, typename tint> void log_end(tint &target, sint &source, genome *current);
    template<typename i> void log_stats(i &interval);
    template<typename sint, typename tint> void log_stats(tint &target, sint &source, genome *current);
    
    template<typename i> void measure(i &interval, genome *g);
    template<typename sint, typename tint> void measure(tint &target, sint &source, genome *current);
    template<typename i> double reduce(i &interval, mpi_local *locf, mpi_local *minf, mpi_local *maxf, mpi_local *sumf);
    
    void minmax(mpi_local *field, mpi_local result, double const &(*func)(double const&, double const&));
    
    #pragma mark EA::FUNCTION::TEMPLATES: evolution cycle

    template<typename f> genome select(f function) { return function(*this); }
    template<typename f, typename v> genome crossover(v &variant, f function) { return function(*this, variant); }
    template<typename f, typename v> void populate(v &variant, f function) { function(variant); }
    template<typename f, typename v> void distribute(v &variant, f function) { function(variant); }
    template<typename f, typename v, typename m> void evolve(v &variant, m &meta, f function) { function(variant, meta); }
    template<typename f, typename v, typename m, typename g> void evolve(v &variant, m &meta, g &individual, f function) { function(variant, meta, individual); }
    
};

template<typename genome> struct obj_eval_stats;
template<typename genome> struct obj_cycle_stats;
template<typename genome> struct obj_run_stats;

template<typename genome> struct objective<genome>::evaluation_interval {
    
    int id = 0;
    int max = 0;
    int log_interval = 0;
    
    bool log_stdout = false;
    bool log_fout = true;
    
    FILE *stats_out;
    
    obj_eval_stats<genome> stats;
    
    genome *local;
    
    mpi_local stat_send;
    mpi_local stat_recv;
    
    void init() {
        this->id = 0;
        this->local = {};
        this->stat_send = {};
        this->stat_recv = {};
    }
    
    void begin(genome *current = NULL) {
        this->id++;
        LOG(9, 0, 0, "BEGIN INTERVAL EVAL %d\r\n", this->id);
        this->stats = {};
        this->stat_send.init();
        this->stat_recv.init();
        this->local = current;
        this->stats.start = MPI_Wtime();
        this->stats.begin_header();
    }

    void end(genome *current = NULL) {
        LOG(9, 0, 0, "END INTERVAL EVAL %d\r\n", this->id);
        this->stats.local_t.value = MPI_Wtime() - this->stats.start;
        this->stats.end_header();
    }
    
    evaluation_interval(): stats() {};

};

template<typename genome> struct objective<genome>::cycle_interval {

    int id = 0;
    int max = 0;
    int log_interval = 0;
    
    bool log_stdout = true;
    bool log_fout = true;
    
    FILE *stats_out;
    
    mpi_local stat_send;
    mpi_local stat_recv;
    
    obj_cycle_stats<genome> stats;
    
    evaluation_interval eval;
    
    genome *local = eval.local;
    
    void init() {
        this->id = 0;
        this->local = {};
        this->stat_send = {};
        this->stat_recv = {};
        this->eval.init();
    }
    
    void begin(genome *current = NULL) {
        LOG(9, 0, 0, "BEGIN INTERVAL CYCLE %d\r\n", this->id);
        this->stats = {};
        this->stat_send.init();
        this->stat_recv.init();
        this->local = current;
        if(this->local) { LOG(2, mpi.id, 0, "GENOME ID %d", this->local->numeric_id()); }
        this->stats.start = MPI_Wtime();
        this->stats.begin_header();
    }

    void end(genome *current = NULL) {
        LOG(9, 0, 0, "END INTERVAL CYCLE %d\r\n", this->id);
        this->stats.local_t.value = MPI_Wtime() - this->stats.start;
        this->stats.end_header();
    }
    
    cycle_interval(): stats(), eval() {};
    
};

template<typename genome> struct objective<genome>::run_interval {
    
    int id = 0;
    int max = 0;
    int log_interval = 0;
    
    bool log_stdout = true;
    bool log_fout = true;
    
    FILE *stats_out;
    
    obj_run_stats<genome> stats;
    
    cycle_interval cycle;
    
    genome *local;
    
    mpi_local stat_send;
    mpi_local stat_recv;
    
    void init() {
        this->stats = {};
        this->local = {};
        this->stat_send = {};
        this->stat_recv = {};
        this->cycle.init();
    }
    
    void begin(genome *current = NULL) {
        this->stats = {};
        this->cycle.stats = {};
        this->cycle.eval.stats = {};
        this->stat_send.init();
        this->stat_recv.init();
        this->local = current;
        this->stats.start = MPI_Wtime();
        if(this->local) { LOG(2, mpi.id, 0, "GENOME ID %d", this->local->numeric_id()); }
        LOG(3, 0, 0, "BEGIN INTERVAL RUN %d\r\n", this->id);
        this->stats.begin_header();
    }

    void end(genome *current = NULL) {
        LOG(3, 0, 0, "END INTERVAL RUN %d\r\n", this->id);
        this->stats.local_t.value = MPI_Wtime() - this->stats.start;
        this->stats.end_header();
    }
    
    run_interval(): stats(), cycle() {};
    
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

#pragma mark EA::OBJECTIVE::FUNCTION filter{}

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

template<typename genome> template<typename i> double objective<genome>::reduce(i &interval, mpi_local *locf, mpi_local *minf, mpi_local *maxf, mpi_local *sumf) {
    
    interval.stat_recv.init();
    MPI_Reduce(locf, &interval.stat_recv, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    if(mpi.id == 0) { this->minmax(minf, interval.stat_recv, std::min<double>); }
    
    interval.stat_recv.init();
    MPI_Reduce(locf, &interval.stat_recv, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    if(mpi.id == 0) { this->minmax(maxf, interval.stat_recv, std::max<double>); }
    
    interval.stat_recv.init();
    MPI_Reduce(locf, sumf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    double avg = 0.0;
    
    if(mpi.id == 0) { avg = sumf->value / mpi.size; }

    return avg;
}

#pragma mark EA::META::OBJECTIVE::FUNCTION: <>measure()

template<typename genome> template<typename sint, typename tint> void objective<genome>::measure(tint &target, sint &source, genome *current) {
    
    target.stats.avg_t = this->reduce(target, &source.stats.local_t, &target.stats.min_t, &target.stats.max_t, &target.stats.sum_t);
    target.stats.avg_scatter_t = this->reduce(target, &source.stats.local_scatter_t, &target.stats.min_scatter_t, &target.stats.max_scatter_t, &target.stats.sum_scatter_t);
    target.stats.avg_gather_t = this->reduce(target, &source.stats.local_gather_t, &target.stats.min_gather_t, &target.stats.max_gather_t, &target.stats.sum_gather_t);
    target.stats.avg_migration_t = this->reduce(target, &source.stats.local_migration_t, &target.stats.min_migration_t, &target.stats.max_migration_t, &target.stats.sum_migration_t);
    
    if(mpi.id != 0) { return; }
    
    target.local->fitness= target.stats.avg_t * -1;
    
    std::vector<genome> valid = this->filter(false);

    if(valid.size() == 0) { return; }
    
    target.stats.best_fitness = valid[0].fitness;
    target.stats.best.push_back(valid[0].fitness);

    if(valid[0].fitness > target.stats.best_fitness) {
        target.stats.best_fitness = valid[0].fitness;
        target.stats.best.push_back(valid[0].fitness);
    }
    
    auto result = minmax_element(valid.begin(), valid.end(),[](genome const &g1, genome const &g2) { return g1.fitness < g2.fitness;});
    
    target.stats.min_fitness = result.first->fitness;
    target.stats.max_fitness = result.second->fitness;
    target.stats.total_fitness = std::accumulate(valid.begin(), valid.end(), 0.0, [](double v, const genome& o){ return o.fitness + v; });
    
    double total_best_fitness = std::accumulate(target.stats.best.begin(), target.stats.best.end(), 0.0);
    
    target.stats.avg_fitness = target.stats.total_fitness / valid.size();
    target.stats.avg_best_fitness = total_best_fitness / target.stats.best.size();
    
    current->measure(target);
    
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
    
    interval.stats.best_fitness = valid[0].fitness;
    interval.stats.best.push_back(valid[0].fitness);

    if(valid[0].fitness > interval.stats.best_fitness) {
        interval.stats.best_fitness = valid[0].fitness;
        interval.stats.best.push_back(valid[0].fitness);
    }
    
    auto result = minmax_element(valid.begin(), valid.end(),[](genome const &g1, genome const &g2) { return g1.fitness < g2.fitness;});
    
    interval.stats.min_fitness = result.first->fitness;
    interval.stats.max_fitness = result.second->fitness;
    
    interval.stats.total_fitness = std::accumulate(valid.begin(), valid.end(), 0.0, [](double v, const genome& o){ return o.fitness + v; });
    double total_best_fitness = std::accumulate(interval.stats.best.begin(), interval.stats.best.end(), 0.0);
    
    interval.stats.avg_fitness = interval.stats.total_fitness / valid.size();
    interval.stats.avg_best_fitness = total_best_fitness / interval.stats.best.size();
    
}

template<typename genome> template<typename i> void objective<genome>::log_stats(i &interval) {
    
    if(interval.id % interval.log_interval == 0) {
        
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

template<typename genome> template<typename sint, typename tint> void objective<genome>::log_stats(tint &target, sint &source, genome *current) {
    
    
    
}

template<typename genome> template<typename i> void objective<genome>::log_begin(i &interval) {
    
    LOG(7, 0, 0, "\r\n--- (%d) END OBJECTIVE CYCLE %d OF %d ---\r\n", mpi.id, interval.id, interval.max);
    LOG(4, mpi.id, 0, "\r\n--- (%d) END OBJECTIVE CYCLE %d OF %d ---\r\n", mpi.id, interval.id, interval.max);
    
}

template<typename genome> template<typename i> void objective<genome>::log_end(i &interval) {
    
    LOG(7, 0, 0, "\r\n--- (%d) END OBJECTIVE CYCLE %d OF %d ---\r\n", mpi.id, interval.id, interval.max);
    LOG(4, mpi.id, 0, "\r\n--- (%d) END OBJECTIVE CYCLE %d OF %d ---\r\n", mpi.id, interval.id, interval.max);
    
    this->measure(interval, interval.local);
    this->log_stats(interval);
    
}

template<typename genome> template<typename sint, typename tint> void objective<genome>::log_end(tint &target, sint &source, genome *current) {
    
    this->measure(target, source, current);
    this->log_stats(target);
    
}

template<typename genome> template<typename i> void objective<genome>::begin(i &interval, genome *current) {
    
    interval.begin(current);
    
    this->log_begin(interval);
    
}

template<typename genome> template<typename i> void objective<genome>::end(i &interval, genome *current) {
    
    interval.end(current);
    
    this->log_end(interval);
    
}

template<typename genome> template<typename sint, typename tint> void objective<genome>::end(tint &target, sint &source, genome *current) {
    
    target.end(current);
    
    this->log_end(target, source, current);
    
}

template<typename genome>
void objective<genome>::minmax(mpi_local *field, mpi_local result, double const &(*func)(double const&, double const&)) {
    *field = (*field).value == 0.0 ? result : func((*field).value, result.value) != (*field).value ? result : *field;
}

#endif /* objective_h */

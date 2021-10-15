//
//  objective.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 8/18/21.
//  Copyright © 2021 Bradley Morgan. All rights reserved.
//

#ifndef objective_h
#define objective_h

#pragma mark DATATYPE: objective

template<typename genome> struct objective {
    
    int id = 0;
    int max_runs = 0;
    int max_evo_evals = 0;
    int max_fit_evals = 0;
    int mu = 0;
    int lambda = 0;
    int islands = 0;
    int workers = 0;
    
    double mutation_rate = 0.0;
    
    objective_run run;
    
    std::vector<genome> population = {};
    
    #pragma mark GENERIC: evolution functions
    
    template<typename f> genome select(f function) { return function(*this); }
    template<typename f, typename v> genome select(f function, v &variant) { return function(*this, variant); }
    
    template<typename f> genome crossover(f function) { return function(*this); }
    template<typename f, typename v> genome crossover(v &variant, f function) { return function(*this, variant); }
    
    template<typename f> void populate(f function) { function(*this); }
    template<typename f, typename v> void populate(v &variant, f function) { function(variant); }
    template<typename f, typename v, typename m> void populate(v &variant, m &meta, f function) { function(*this, variant, meta); }
    
    template<typename f> void distribute(f function) { function(*this); }
    template<typename f, typename v> void distribute(v &variant, f function) { function(variant); }
    
    template<typename f, typename e> void evaluate(e &ea, genome &individual, f function) { function(ea, individual); }
    template<typename f, typename m, typename v> void evaluate(m &meta, v &variant, genome &individual, f function) { function(meta, variant, individual); }
    
    template<typename f> void calculate(f function) { function(); }
    template<typename f, typename v> void calculate(f function, v &variant) { function(*this, variant); }
    
    template<typename f> void evolve(f function) { function(*this); }
    template<typename f, typename v> void evolve(v &variant, f function) { function(*this, variant); }
    template<typename f, typename v, typename m> void evolve(v &variant, m &meta, f function) { function(variant, meta); }
    template<typename f, typename v, typename m, typename g> void evolve(v &variant, m &meta, g &individual, f function) { function(variant, meta, individual); }
    
    template<typename e> void log_begin(objective_run &run, objective<genome> &o, e &ea);
    template<typename e> void log_end(objective_run &run, objective<genome> &o, e &ea);
    template<typename e> void log_begin(objective_eval &eval, objective<genome> &o, e &ea);
    template<typename e> void log_end(objective_eval &eval, objective<genome> &o, e &ea);
    
//    template<typename f, typename e> void begin(objective_run &run, e &ea) { run.begin(ea, *this); }
//    template<typename f, typename e> void end(objective_run &run, e &ea) { run.end(ea, *this); }
//
//    template<typename f, typename e> void begin(objective_eval &eval, e &ea) { run.eval.begin(ea, *this); }
//    template<typename e> void end(objective_eval &eval, e &ea) { run.eval.end(ea, *this); }
    
    void cpd();
    void fitness();
    
    #pragma mark DATATYPE: metrics
    
    struct metrics {
        
        struct values {
            
            std::vector<double> cpd = {};
            double fitness = 0.0;
            
            values() : cpd(0.0), fitness(0.0) {}
            
        };
        
        metrics(void) {}
        
        values value;
        
    };
    
    metrics aggregate;
    
};

// comparator for parent fitness values ...

//template<typename genome> bool compare_fitness(const genome &p1, const genome &p2) {
//    return p1.fitness < p2.fitness;
//}


template<typename genome> void objective<genome>::fitness() {
    
    this->aggregate.value.fitness = 0.0;
    
    for(typename std::vector<genome>::iterator it = this->population.begin(); it != this->population.end(); ++it) {
        LOG(10, 0, 0, "genome %d fitness = %f\r\n", this->id, it->fitness);
        this->aggregate.value.fitness += it->fitness;
    }
    
}

template<typename genome> void objective<genome>::cpd() {

    LOG(10, 0, 0, "generating cpd\r\n");

    double cumulative_probability = 0.0;

    this->aggregate.value.cpd.clear();

    LOG(6, 0, 0, "calculating total fitness ...\r\n");

    LOG(6, 0, 0, "sorting population descending fitness ...\r\n");

    std::sort(this->population.begin(), this->population.end(), compare_fitness<genome>);
    std::reverse(this->population.begin(), this->population.end());

    LOG(6, 0, 0, "objective %d cpd calculation total fitness = %f", this->id, this->aggregate.value.fitness);

    for(int i=0; i < this->population.size(); i++) {

        this->population[i].selection_distribution = this->population[i].fitness / this->aggregate.value.fitness;

        LOG(8, 0, 0, "calculated island %d solution %d fitness %f selection distribution = %f\r\n", this->id, i, this->population[i].fitness, this->population[i].selection_distribution);

        cumulative_probability += this->population[i].selection_distribution;

        LOG(8, 0, 0, "solution %d cumulative prob = %f\r\n", i, cumulative_probability);

        this->aggregate.value.cpd.push_back(cumulative_probability);

    }

}



template<typename e> void log_begin(objective_run &run, objective<topology> &obj, e &meta) {
    
    LOG(2, meta.variant.isle.id, 0, "\r\n --- BEGIN META RUN %d ---\r\n", meta.run.id);

    meta.run.eval.id = 1;
    meta.run.start = MPI_Wtime();
    meta.topologies.population.clear();
    meta.variant.isle.population.clear();
    meta.run.eval.stats.init();
    
}

template<typename e> void log_end(objective_run &run, objective<topology> &obj, e &meta) {
    
    if(meta.variant.isle.id != 0) { return; }

    LOG(2, meta.variant.isle.id, 0, "\r\n--- END META RUN %d ---\r\n", meta.run.id);

    double run_end = MPI_Wtime();

    meta.run.stats.run_duration = run_end - meta.run.start;

    fprintf(config::topo_run_stats_out, "average_topo_fitness, global_best_topo_id, global_best_topo_rounds, global_best_topo_channels, global_best_topo_round_fitness, global_best_topo_fitness1, local_best_topo_fitness, global_best_topo_fitness2, average_local_best_topo_fitness, average_global_best_topo_fitness, t_id, t_rounds, t_channels, t_fitness\r\n");

    std::fprintf(config::topo_run_stats_out, "%d,%f,%f,%f,%f,%f,%d", meta.run.id, meta.run.stats.run_duration, meta.run.eval.stats.average_local_best_topo_fitness, meta.run.eval.stats.average_global_best_topo_fitness, meta.run.eval.stats.global_best_topo_fitness, meta.run.eval.stats.total_migrate_time, meta.run.stats.total_channels);

    fflush(config::topo_run_stats_out);

}




template<typename e> void log_begin(objective_run &run, objective<solution> &obj, e &solver) {
    
    LOG(2, solver.variant.isle.id, 0, "   --- BEGIN SOLVER RUN %d ---\r\n", solver.run.id);

    solver.run.eval.id = 1;
    solver.run.start = MPI_Wtime();
    solver.solutions.population.clear();
    solver.variant.isle.population.clear();
    solver.variant.isle.population.resize(config::island_mu);
    solver.run.eval.stats.init();

    LOG(2, solver.variant.isle.id, 0, "\r\n%5s %6s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s", "r", "e", "avg_fit", "lbest_fit", "gbest_fit", "avg_lbest", "avg_gbest", "avg_scat_t", "avg_gathr_t", "avg_migr", "init_t", "eval_t");
    
}

template<typename e> void log_end(objective_run &run, objective<solution> &obj, e &solver) {

    if(solver.variant.isle.id != 0) { return; }

    double run_end = MPI_Wtime();

    solver.run.stats.run_duration = run_end - solver.run.start;

    LOG(2, solver.variant.isle.id, 0, "\r\n   --- END SOLVER RUN %d ---\r\n", solver.run.id);

    std::fprintf(config::run_stats_out, "%d,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d\r\n", solver.run.id, solver.run.eval.stats.global_best_fitness, solver.run.eval.stats.average_local_best_fitness, solver.run.eval.stats.average_global_best_fitness, solver.run.eval.stats.total_scatter_time, solver.run.eval.stats.total_gather_time, solver.run.eval.stats.total_migrate_time, solver.run.stats.run_duration, solver.run.stats.init_duration, solver.variant.islands, solver.variant.island_size);

    fflush(config::run_stats_out);
    
}

template<typename e> void log_begin(objective_eval &eval, objective<solution> &obj, e &solver) {
    //LOG(2, solver.variant.isle.id, 0, "\r\n  --- BEGIN SOLVER EVOLUTION %d (RUN %d) ---\r\n\r\n", solver.run.eval.id, solver.run.id);
}

template<typename e> void log_end(objective_eval &eval, objective<solution> &obj, e &solver) {
    //LOG(2, solver.variant.isle.id, 0, "\r\n  --- END SOLVER EVOLUTION %d (RUN %d) ---\r\n\r\n", solver.run.eval.id, solver.run.id);
}

template<typename e> void log_begin(objective_eval &eval, objective<topology> &obj, e &meta) {
    meta.run.eval.start = MPI_Wtime();
    LOG(2, meta.variant.isle.id, 0, "\r\n  --- BEGIN META GENOME %d (RUN %d) ---\r\n\r\n", meta.run.eval.id, meta.run.id);
}

template<typename e> void log_end(objective_eval &eval, objective<topology> &obj, e &meta) {
    meta.run.eval.stats.eval_duration = MPI_Wtime() - meta.run.eval.start;
    LOG(2, meta.variant.isle.id, 0, "\r\n  --- END META GENOME %d (RUN %d) duration = %f ---\r\n", meta.run.eval.id, meta.run.id, meta.run.eval.stats.eval_duration);
}





#endif /* objective_h */

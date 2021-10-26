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
    
    int id = 1;
    int max_runs = 0;
    int max_evo_cycles = 0;
    int max_fit_evals = 0;
    int mu = 0;
    int lambda = 0;
    int islands = 0;
    int workers = 0;
    
    double mutation_rate = 0.0;
    
    objective_run run;
    evolution_cycle cycle;
    
    std::vector<genome> population = {};
    
    
    #pragma mark DATATYPE: @metrics{}
    
    struct metrics {
        
        struct values {
            
            std::vector<double> cpd = {};
            double fitness = 0.0;
            
            values() : cpd(0.0), fitness(0.0) {}
            
        };
        
        metrics(void) {}
        
        values value;
        
    };
    
    #pragma mark EA::FUNCTION::TEMPLATES: metrics
    
    metrics aggregate;
    
    void cpd();
    void fitness();
    
    #pragma mark EA::FUNCTION::TEMPLATES: initialization
    
    objective<genome>() { this->id = 1; this->run.id = 1; this->run.eval.id = 1; }
    
    template<typename e> void begin(objective_run &run, e &ea);
    template<typename e> void begin(objective_eval &eval, e &ea);
    template<typename e> void begin(evolution_cycle &cycle, e &ea);
    template<typename e> void end(objective_run &run, e &ea);
    template<typename e> void end(objective_eval &eval, e &ea);
    template<typename e> void end(evolution_cycle &cycle, e &ea);
    template<typename e> void end(e &ea) { ea_end(ea, *this); }
    
    #pragma mark EA::FUNCTION::TEMPLATES: initialization logging
    
    template<typename e> void log_begin(objective_run &run, e &ea);
    template<typename e> void log_begin(objective_eval &eval, e &ea);
    template<typename e> void log_begin(evolution_cycle &cycle, e &ea);
    template<typename e> void log_end(objective_run &run, e &ea);
    template<typename e> void log_end(objective_eval &eval, e &ea);
    template<typename e> void log_end(evolution_cycle &cycle, e &ea);
    template<typename e, typename m, typename g> void log_stats(objective_eval &eval, e &solver, m &meta, g &individual);
    template<typename e, typename m, typename g> void log_stats(evolution_cycle &cycle, e &solver, m &meta, g &individual);
    
    #pragma mark EA::FUNCTION::TEMPLATES: evolution cycle

    template<typename f> genome select(f function) { return function(*this); }
    template<typename f, typename v> genome crossover(v &variant, f function) { return function(*this, variant); }
    template<typename f, typename v> void populate(v &variant, f function) { function(variant); }
    template<typename f, typename v> void distribute(v &variant, f function) { function(variant); }
    template<typename f, typename v, typename m> void evolve(v &variant, m &meta, f function) { function(variant, meta); }
    template<typename f, typename v, typename m, typename g> void evolve(v &variant, m &meta, g &individual, f function) { function(variant, meta, individual); }

    
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

    LOG(8, 0, 0, "objective %d cpd calculation total fitness = %f\r\n", this->id, this->aggregate.value.fitness);

    for(int i=0; i < this->population.size(); i++) {

        this->population[i].selection_distribution = this->population[i].fitness / this->aggregate.value.fitness;

        LOG(8, 0, 0, "calculated island %d solution %d fitness %f selection distribution = %f\r\n", this->id, i, this->population[i].fitness, this->population[i].selection_distribution);

        cumulative_probability += this->population[i].selection_distribution;

        LOG(8, 0, 0, "%d=%f|", i, cumulative_probability);

        this->aggregate.value.cpd.push_back(cumulative_probability);;

    }

}

#pragma mark EA::FUNCTION::TEMPLATE: <*>cpd{}

// calculates an arbitrary population's cumulative probability distribution

//template<typename genome, typename e> void objective<genome>::begin(objective_run &run, e &ea);

//template<typename e> void log_end(objective_run &run, e &meta) {
//    
//    if(meta.variant.isle.id != 0) { return; }
//
//    LOG(2, meta.variant.isle.id, 0, "\r\n--- END META RUN %d ---\r\n", meta.run.id);
//
//    double run_end = MPI_Wtime();
//
//    meta.run.eval.id = 1;
//    meta.run.stats.run_duration = run_end - meta.run.start;
//
//    fprintf(config::topo_run_stats_out, "average_topo_fitness, global_best_topo_id, global_best_topo_rounds, global_best_topo_channels, global_best_topo_round_fitness, global_best_topo_fitness1, local_best_topo_fitness, global_best_topo_fitness2, average_local_best_topo_fitness, average_global_best_topo_fitness, t_id, t_rounds, t_channels, t_fitness\r\n");
//
//    std::fprintf(config::topo_run_stats_out, "%d,%f,%f,%f,%f,%f,%d", meta.run.id, meta.run.stats.run_duration, meta.run.eval.stats.average_local_best_topo_fitness, meta.run.eval.stats.average_global_best_topo_fitness, meta.run.eval.stats.global_best_topo_fitness, meta.run.eval.stats.total_migrate_time, meta.run.stats.total_channels);
//
//    fflush(config::topo_run_stats_out);
//
//}

//template<typename e> void log_begin(objective_eval &eval, e &solver) {
//    LOG(5, this->variant.isle.id, 0, "\r\n\r\n  --- BEGIN SOLVER EVOLUTION CYCLE %d (RUN %d) ---\r\n\r\n", this->run.eval.id, solver.solutions.run.id);
//}
//
//template<typename genome> void objective<genome>::log_end(objective_eval &eval) {
//    LOG(5, this->variant.isle.id, 0, "\r\n  --- END SOLVER EVOLUTION CYCLE %d (RUN %d) ---\r\n\r\n", this->run.eval.id, solver.solutions.run.id);
//}
//
//template<typename e> void log_begin(objective_eval &eval, objective<topology> &obj, e &meta) {
//    LOG(2, meta.variant.isle.id, 0, "\r\n  --- BEGIN META GENOME %d (RUN %d) ---\r\n", meta.run.eval.id, meta.run.id);
//}
//
//template<typename e> void log_end(objective_eval &eval, objective<topology> &obj, e &meta) {
//    meta.run.eval.stats.eval_duration = MPI_Wtime() - meta.run.eval.start;
//    LOG(2, meta.variant.isle.id, 0, "\r\n  --- END META GENOME %d (RUN %d) duration = %f ---\r\n", meta.run.eval.id, meta.run.id, meta.run.eval.stats.eval_duration);
//}

#endif /* objective_h */

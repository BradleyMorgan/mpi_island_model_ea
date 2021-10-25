//
//  solution.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 4/12/21.
//  Copyright © 2021 Bradley Morgan. All rights reserved.
//

#ifndef solution_h
#define solution_h

struct offset {
    
    std::array<double, DIM> input;
    
    void generate(double min, double max, double step);
    
};

void offset::generate(double min, double max, double step) {
    
    std::array<double, DIM> offsets;
    std::vector<double> increments;
    
    for (double increment = min; increment <= max; increment += step) {
        increments.push_back(increment);
    }
    
    for(int i=0; i<DIM; i++) {
        offsets[i] = increments[rand()%increments.size()];
    }
    
    this->input = offsets;
    
}

std::array<double, DIM> generate_offsets(double min, double max, double step) {
    
    LOG(6, 0, 0, "generating offsets...\r\n");
    
    std::array<double, DIM> offsets;
    std::vector<double> increments;
    
    for (double increment = min; increment <= max; increment += step) {
        increments.push_back(increment);
    }
    
    for(int i=0; i<DIM; i++) {
        offsets[i] = increments[rand()%increments.size()];
    }
    
    return offsets;
    
}

template<> template<typename e> void objective<solution>::begin(objective_run &run, e &solver) {
    
    LOG(3, 0, 0, "BEGIN ISLAND %d SOLVER OBJECTIVE %d RUN %d -> ", solver.variant.isle.id, this->id, this->run.id);

    this->run.begin();

    this->population.clear();
    //this->population.resize(this->mu);
    
    solver.variant.isle.population.clear();
    solver.variant.isle.population.resize(config::island_mu);
    
    //solver.variant.isle.init();
    
    this->log_begin(run, solver);
    
}

template<> template<typename e> void objective<solution>::end(objective_run &run, e &solver) {
     
    LOG(3, 0, 0, "END ISLAND %d SOLVER OBJECTIVE %d RUN %d\r\n\r\n", solver.variant.isle.id, this->id, this->run.id);

    this->run.end();

    if(solver.variant.isle.id == 0) {
        solver.offsets = generate_offsets(-2.5, 2.5, .5);
    }
    
    MPI_Bcast(&solver.offsets, DIM, MPI_DOUBLE, 0, solver.variant.tcomm);
    
    this->log_end(run, solver);
    
}

//template<typename e> void log_begin(objective_run &run, objective<solution> &obj, solver &solver) {
//
//    LOG(2, solver.variant.isle.id, 0, "   --- BEGIN SOLVER RUN %d ---\r\n", solver.solutions.run.id);
//
//    solver.run.eval.id = 1;
//    solver.run.start = MPI_Wtime();
//    solver.solutions.population.clear();
//    solver.variant.isle.population.clear();
//    solver.variant.isle.population.resize(config::island_mu);
//    solver.run.eval.stats.init();
//
//    LOG(2, solver.variant.isle.id, 0, "\r\n%5s %6s %16s %16s %16s %16s %16s %16s %16s %16s %16s", "r", "e", "avg_fit", "lbest_fit", "gbest_fit", "avg_lbest", "avg_gbest", "avg_scat_t", "avg_gathr_t", "avg_migr", "eval_t");
//
//}

template<> template<typename e> void objective<solution>::begin(objective_eval &eval, e &solver) {
    
    LOG(3, 0, 0, "BEGIN ISLAND %d SOLVER OBJECTIVE %d EVAL %d -> ", solver.variant.isle.id, this->id, this->run.eval.id);

    this->run.eval.begin();
    
    this->log_begin(eval, solver);
    
}

template<> template<typename e> void objective<solution>::end(objective_eval &eval, e &solver) {
    
    LOG(3, 0, 0, "END ISLAND %d SOLVER OBJECTIVE %d EVAL %d -> ", solver.variant.isle.id, this->id, this->run.eval.id);

    this->run.eval.end();
    
    this->log_end(eval, solver);
    
}

template<> template<typename e> void objective<solution>::log_begin(objective_run &run, e &solver) {
    
    if(solver.variant.isle.id != 0 || this->run.id == 0) { return; }
    
    LOG(2, solver.variant.isle.id, 0, "%5s %6s %16s %16s %16s %16s %16s %16s %16s %16s %16s\r\n", "r", "e", "avg_fit", "lbest_fit", "gbest_fit", "avg_lbest", "avg_gbest", "avg_scat_t", "avg_gathr_t", "avg_migr", "eval_t");
    
}


template<> template<typename e> void objective<solution>::log_end(objective_run &run, e &solver) {

    if(solver.variant.isle.id != 0 || this->run.id == 0) { return; }

    std::fprintf(config::run_stats_out, "%d,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d\r\n", solver.solutions.run.id, solver.run.eval.stats.global_best_fitness, solver.run.eval.stats.average_local_best_fitness, solver.run.eval.stats.average_global_best_fitness, solver.run.eval.stats.total_scatter_time, solver.run.eval.stats.total_gather_time, solver.run.eval.stats.total_migrate_time, solver.run.stats.run_duration, solver.run.stats.init_duration, solver.variant.islands, solver.variant.island_size);

    fflush(config::run_stats_out);

}

//
//template<typename e> void log_end(objective_run &run, objective<solution> &obj, e &solver) {
//
//    if(solver.variant.isle.id != 0) { return; }
//
//    double run_end = MPI_Wtime();
//
//    solver.run.stats.run_duration = run_end - solver.run.start;
//
//    LOG(2, solver.variant.isle.id, 0, "\r\n\r\n   --- END SOLVER RUN %d ---\r\n", solver.solutions.run.id);
//
//    std::fprintf(config::run_stats_out, "%d,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d\r\n", solver.solutions.run.id, solver.run.eval.stats.global_best_fitness, solver.run.eval.stats.average_local_best_fitness, solver.run.eval.stats.average_global_best_fitness, solver.run.eval.stats.total_scatter_time, solver.run.eval.stats.total_gather_time, solver.run.eval.stats.total_migrate_time, solver.run.stats.run_duration, solver.run.stats.init_duration, solver.variant.islands, solver.variant.island_size);
//
//    fflush(config::run_stats_out);
//
//}

template<> template<typename e> void objective<solution>::log_begin(objective_eval &eval, e &solver) {
    
    if(solver.variant.isle.id != 0 || this->run.eval.id == 0) { return; }
    
    LOG(6, 0, 0, "\r\n --- BEGIN SOLVER EVOLUTION CYCLE %d ---\r\n", this->run.eval.id);
    
}

template<> template<typename e> void objective<solution>::log_end(objective_eval &eval, e &solver) {
    
    if(solver.variant.isle.id != 0 || this->run.eval.id == 0) { return; }
    
    LOG(6, 0, 0, "\r\n --- END SOLVER EVOLUTION CYCLE %d ---\r\n", this->run.eval.id);
    
    
}

template<> template<typename e, typename m, typename g> void objective<solution>::log_stats(objective_eval &eval, e &solver, m &meta, g &genome) {
    
    if(solver.variant.isle.id != 0 || solver.solutions.run.eval.id == 0) { return; }
    
    log_fn_eval_stats(solver, meta, genome);
    
}

//template<> template<typename e> void objective<solution>::end(e &solver) {
//    
//    ea_end(solver, *this);
//    
//}

#endif /* solution_h */

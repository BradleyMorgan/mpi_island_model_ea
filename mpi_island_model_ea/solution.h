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


/*----------------------------------------------------------------------------*/

#pragma mark EA::FUNCTION::TEMPLATES:INITIALIZATION: <solution>

// specialization of generic @objective{} methods for performing pre and post
// interval initialization and cleanup


#pragma mark EA::OBJECTIVE::RUN <solution>::begin()

// call at the start of an experimental run, usually after incrementing the run
// count near the beginning of a loop
//
// a new run marks the beginning of an experiment, so perform appropriate
// initializations for restarting the algorithm (e.g. populations and per-run
// metrics)


template<> template<typename e> void
objective<solution>::begin(objective_run &run, e &solver) {
    
    LOG(3, 0, 0, "BEGIN ISLAND %d SOLVER OBJECTIVE %d RUN %d -> ", \
        solver.variant.isle.id, this->id, this->run.id);

    this->run.begin();

    this->population.clear();
    
    solver.variant.isle.population.clear();
    solver.variant.isle.population.resize(config::island_mu);
    
    this->log_begin(run, solver);
    
}

#pragma mark EA::OBJECTIVE::CYCLE <solution>::begin()

// call at start of evolution cycle (generation)
//
// evolution cycle performs parent selection, crossover, child evaluation (if
// needed) and survival selection perform pre-evolution operations
//

template<> template<typename e>
void objective<solution>::begin(evolution_cycle &cycle, e &solver) {
    
    LOG(3, 0, 0, "BEGIN ISLAND %d SOLVER OBJECTIVE %d EVOLUTION CYCLE %d\r\n", \
        solver.variant.isle.id, this->id, this->cycle.id);
    
    this->run.eval.begin();

}

#pragma mark EA::OBJECTIVE::EVAL <solution>::begin()

// call before evaluating an individual, usually after incrementing
// the evaluation count, near the beginning of a loop
//
// an evaluation encompasses all tasks necessary to evaluate individual fitness
// use eval.begin() to gather statistics, perform measurements, logging, etc.

template<> template<typename e>
void objective<solution>::begin(objective_eval &eval, e &solver) {
    
    LOG(3, 0, 0, "BEGIN ISLAND %d SOLVER OBJECTIVE %d EVAL %d -> ", \
        solver.variant.isle.id, this->id, this->run.eval.id);

    this->run.eval.begin();
    
    this->log_begin(eval, solver);
    
}

/*----------------------------------------------------------------------------*/

#pragma mark EA::OBJECTIVE::RUN <solution>::end()

// call at the end of an experimental run, usually before incrementing the
// run count, near the end of a loop a new run marks the beginning of an
// experiment, so perform appropriate initializations for restarting the
// algorithm

template<> template<typename e>
void objective<solution>::end(objective_run &run, e &solver) {
     
    LOG(3, 0, 0, "END ISLAND %d SOLVER OBJECTIVE %d RUN %d\r\n\r\n", \
        solver.variant.isle.id, this->id, this->run.id);
    
    this->log_end(run, solver);
    
    this->run.end();

    if(solver.variant.isle.id == 0) {
        solver.offsets = generate_offsets(-2.5, 2.5, .5);
    }
    
    MPI_Bcast(&solver.offsets, DIM, MPI_DOUBLE, 0, solver.variant.tcomm);
    
}

#pragma mark EA::OBJECTIVE::CYCLE <solution>::end()

// call at end of evolution cycle (generation)
//
// evolution cycle performs parent selection, crossover, child evaluation (if
// needed) and survival selection perform post-evolution operations

template<> template<typename e>
void objective<solution>::end(evolution_cycle &cycle, e &solver) {
    
    LOG(3, 0, 0, "END ISLAND %d SOLVER OBJECTIVE %d EVOLUTION CYCLE %d -> ",
        solver.variant.isle.id, this->id, this->cycle.id);
    
    this->run.eval.end();

}

#pragma mark EA::OBJECTIVE::EVAL <solution>::end()

// call before evaluating an individual, usually after incrementing the
// evaluation count, near the beginning of a loop an evaluation encompasses all
// tasks necessary to evaluate individual fitness use eval.begin() to gather
// statistics, perform measurements, logging, etc.

template<> template<typename e>
void objective<solution>::end(objective_eval &eval, e &solver) {
    
    LOG(3, 0, 0, "END ISLAND %d SOLVER OBJECTIVE %d EVAL %d -> ",
        solver.variant.isle.id, this->id, this->run.eval.id);

    this->run.eval.end();
    
    this->log_end(eval, solver);
    
}

/*----------------------------------------------------------------------------*/

#pragma mark EA::FUNCTION::TEMPLATES:LOGGING: <solution>

// specialization of generic @objective{} methods for file and standard output


#pragma mark EA::OBJECTIVE::RUN <solution>::log_begin()

// log anything of interest at the start of an experimental run

template<> template<typename e>
void objective<solution>::log_begin(objective_run &run, e &solver) {
    
    if(solver.variant.isle.id != 0 || this->run.id == 0) { return; }
    
    LOG(2, solver.variant.isle.id, 0, "%5s %6s %16s %16s %16s %16s %16s %16s"
        "%16s %16s %16s\r\n", "r", "e", "avg_fit", "lbest_fit", "gbest_fit",
        "avg_lbest", "avg_gbest", "avg_scat_t", "avg_gathr_t", "avg_migr",
        "eval_t");
    
}

#pragma mark EA::OBJECTIVE::CYCLE <solution>::log_begin()

// log anything of interest at the start of an evolution cycle

template<> template<typename e>
void objective<solution>::log_begin(evolution_cycle &cycle, e &solver) {
    
    if(solver.variant.isle.id != 0 || this->run.id == 0) { return; }
    
    LOG(3, 0, 0, "LOG SOLVER OBJECTIVE %d EVOLUTION CYCLE %d BEGIN\r\n",
        this->id, this->run.id);
    
}

#pragma mark EA::OBJECTIVE::EVAL <solution>::log_begin()

// log anything of interest before evaluating an individual

template<> template<typename e>
void objective<solution>::log_begin(objective_eval &eval, e &solver) {
    
    if(solver.variant.isle.id != 0 || this->run.eval.id == 0) { return; }
    
    LOG(6, 0, 0, "\r\n --- BEGIN SOLVER EVOLUTION CYCLE %d ---\r\n",
        this->run.eval.id);
    
}

#pragma mark EA::OBJECTIVE::RUN <solution>::log_end()

// log anything of interest at the end of an experimental run

template<> template<typename e>
void objective<solution>::log_end(objective_run &run, e &solver) {

    if(solver.variant.isle.id != 0 || this->run.id == 0) { return; }

    std::fprintf(config::run_stats_out, "%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%d,"
        "%d\r\n",
                 solver.solutions.run.id,
                 solver.solutions.cycle.id,
                 solver.solutions.run.eval.id,
                 solver.run.eval.stats.global_best_fitness,
                 solver.run.eval.stats.average_local_best_fitness,
                 solver.run.eval.stats.average_global_best_fitness,
                 solver.run.eval.stats.total_scatter_time,
                 solver.run.eval.stats.total_gather_time,
                 solver.run.eval.stats.total_migrate_time,
                 solver.solutions.run.stats.run_duration,
                 solver.init_duration,
                 solver.variant.islands,
                 solver.variant.island_size);

    fflush(config::run_stats_out);

}

#pragma mark EA::OBJECTIVE::CYCLE <solution>::log_end()

// log anything of interest at the end of an evolution cycle

template<> template<typename e>
void objective<solution>::log_end(evolution_cycle &cycle, e &solver) {
    
    if(solver.variant.isle.id != 0 || this->cycle.id == 0) { return; }
    
    LOG(3, 0, 0, "LOG META OBJECTIVE %d EVOLUTION CYCLE %d END\r\n",
        this->id, this->run.id);
    
}

#pragma mark EA::OBJECTIVE::EVAL <solution>::log_end()

// log anything of interest at the end of an evaluation

template<> template<typename e>
void objective<solution>::log_end(objective_eval &eval, e &solver) {
    
    if(solver.variant.isle.id != 0 || this->run.eval.id == 0) { return; }
    
    LOG(6, 0, 0, "\r\n --- END SOLVER EVOLUTION CYCLE %d ---\r\n",
        this->run.eval.id);
    
}

#pragma mark EA::OBJECTIVE::CYCLE <solution>::log_stats()

// custom evolution cycle logging

template<> template<typename e, typename m, typename g>
void objective<solution>::log_stats(evolution_cycle &cycle,e &solver, m &meta,
                                    g &genome) {
    
    if(solver.variant.isle.id != 0 || this->cycle.id == 0) { return; }
    
    log_fn_cycle_stats(solver, meta, genome);
    
}

#pragma mark EA::OBJECTIVE::EVAL <solution>::log_stats()

// catchall function for custom logging

template<> template<typename e, typename m, typename g>
void objective<solution>::log_stats(objective_eval &eval, e &solver, m &meta,
                                    g &genome) {
    
    if(solver.variant.isle.id != 0 || this->run.eval.id == 0) { return; }
    
    log_fn_eval_stats(solver, meta, genome);
    
}



//template<> template<typename e> void objective<solution>::end(e &solver) {
//    
//    ea_end(solver, *this);
//    
//}


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





#endif /* solution_h */

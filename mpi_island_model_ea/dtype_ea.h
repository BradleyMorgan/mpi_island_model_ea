//
//  dtype_ea.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 10/13/21.
//  Copyright © 2021 Bradley Morgan. All rights reserved.
//

#ifndef dtype_ea_h
#define dtype_ea_h

#include "solution.h"

struct ea_stats: time_stats, parallel_stats, fitness_stats {
    
    void init() {
        
        LOG(5, 0, 0, "\r\n--- INIT EA STATS ---\r\n");
        
        this->start = 0.0;
        this->duration = 0.0;
        this->init_duration = 0.0;
        
    }
    
};

#pragma mark DATATYPE: @ea_eval{}

// track eval stats, etc. per island

struct ea_eval {
  
    int id = 0;
    int max = 0;
    int log_interval = 0;
    
    eval_stats stats;
  
    ea_eval(): id(0) {};
    
    void begin() {
        LOG(3, 0, 0, "\r\n--- BEGIN EA EVAL %d OF %d ---\r\n", this->id, this->max);
        this->stats.start = MPI_Wtime();
    }

    void end() {
        double eval_end = MPI_Wtime();
        this->stats.duration = eval_end - this->stats.start;
        LOG(3, 0, 0, "\r\n--- END EA EVAL %d OF %d ---\r\n", this->id, this->max);
    }

};

struct ea_cycle {
    
    int id = 0;
    int max = 0;
    int log_interval = 0;
    
    ea_eval eval;
    
    cycle_stats stats;
    
    void begin() {
        LOG(5, 0, 0, "\r\n--- BEGIN EA CYCLE %d OF %d ---\r\n", this->id, this->max);
        this->stats.duration = 0.0;
        this->stats.start = MPI_Wtime();
    }
    
    void end() {
        LOG(5, 0, 0, "\r\n--- END EA CYCLE %d OF %d ---\r\n", this->id, this->max);
        double cycle_end = MPI_Wtime();
        this->stats.duration = cycle_end - this->stats.start;
    }
    
    
};

#pragma mark EA::DATATYPE: @ea_run{}

// track run stats, etc. per island

struct ea_run {
    
    int id = 0;
    int max = 0;
    int log_interval = 0;
    
    ea_cycle cycle;
    run_stats stats;
  
    ea_run(): id(0) {};
    
    void begin() {
        this->id++;
        LOG(4, 0, 0, "\r\n--- BEGIN EA RUN %d OF %d ---\r\n", this->id, this->max);
        this->stats.duration = 0.0;
        this->stats.start = MPI_Wtime();
    }
    
    void end() {
        LOG(4, 0, 0, "\r\n--- END EA RUN %d OF %d ---\r\n", this->id, this->max);
        double run_end = MPI_Wtime();
        this->stats.duration = run_end - this->stats.start;
    }
    
};

// ----- start mpi derived datatypes...

template<typename m> void parallel_types(m &model) {
    
    MPI_Datatype ptype;

    MPI_Type_contiguous(128, MPI_CHAR, &ptype);
    MPI_Type_commit(&ptype);

    MPI_Datatype sol_types[11] = { MPI_CHAR, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, ptype };

    MPI_Aint sol_offsets[11] = { 0,                                  // id                       + 1 lluint
        sizeof(char)*64,                                             // input                    + 10 double
        sizeof(double)*DIM + sizeof(char)*64,                        // fitness                  + 11 double
        sizeof(double)*(DIM+1) + sizeof(char)*64,                    // selection_distribution   + 12 double
        sizeof(double)*(DIM+2) + sizeof(char)*64,                    // group                    + 13 double
        sizeof(double)*(DIM+3) + sizeof(char)*64,                    // source                   + 1 int
        sizeof(double)*(DIM+3) + sizeof(char)*64 + (sizeof(int)),    // locale                   + 1 int
        sizeof(double)*(DIM+3) + sizeof(char)*64 + (sizeof(int)*2),  // migrations               + 1 int
        sizeof(double)*(DIM+3) + sizeof(char)*64 + (sizeof(int)*3),  // selected                 + 1 int
        sizeof(double)*(DIM+3) + sizeof(char)*64 + (sizeof(int)*4),  // survival                 + 1 int
        sizeof(double)*(DIM+3) + sizeof(char)*64 + (sizeof(int)*5)   // parents                  + 2 int
    };

    int sol_lengths[11] = { 64, DIM, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

    MPI_Type_create_struct(11, sol_lengths, sol_offsets, sol_types, &model.solution_type);
    MPI_Type_commit(&model.solution_type);

    MPI_Datatype visa_types[4] = { MPI_INT, MPI_INT, MPI_INT, MPI_CHAR };
    MPI_Aint visa_offsets[4] = { 0, sizeof(int), sizeof(int)*2, sizeof(int)*3 };

    int visa_lengths[4] = { 1, 1, 1, 64 };

    MPI_Type_create_struct(4, visa_lengths, visa_offsets, visa_types, &model.visa_type);
    MPI_Type_commit(&model.visa_type);

}

#pragma mark DATATYPE: @ea_variant{}

// describes the properties needed for parallel (island model) EA

struct ea_model {

    unsigned long start;

    char name[128];
    
    int island_size;
    int root = 0;

    // isle holds context specific population and migration data
    // from the perspective of a specific island instance (e.g. mpi rank)
    //
    // @island.h::island{} datatype
    //
    
    island isle;
    
    // TODO: maintaining a list of island_ids doesn't seem necessary
    
    MPI_Datatype solution_type;
    MPI_Datatype visa_type;
    
    MPI_Comm tcomm;
    
    void init() {

        this->start = MPI_Wtime();

        // load configuration items ...

        config::load("config.txt", mpi.size, mpi.id);

        this->tcomm = MPI_COMM_WORLD;
        this->island_size = config::mu_sub;
        
        parallel_types(*this);
        
        this->isle.init();
        
    }
    
};

#pragma mark EA::DATATYPE: @ea{}

// wrapper type to serve as the trunk for ea properties
//
//  +----heirarchy---+----context-----+-----termination-----+---scope---+
//  | EA𝑛            | ea𝑛 logic      | all objective runs  | O₁₋𝑛,r,c,e |
//  |  ↳RUN𝑛         | ea𝑛 execution  | ea runtime          | O₁₋𝑛,c,e   |
//  |   ↳CYCLE𝑛      | ea𝑛 evolution  | max int, fitness    | O₁₋𝑛,e     |
//  |    ↳EVAL𝑛      | ea𝑛 population | fitness result      | O₁₋𝑛       |
//  |  ↳OBJECTIVE𝑛   | o𝑛 logic       | obj termination     | O𝑛,r,c,e   |
//  |    ↳RUN𝑛       | o𝑛 execution   | obj runtime         | O𝑛,c,e     |
//  |      ↳CYCLE𝑛   | o𝑛 evolution   | max int, fitness    | O𝑛,e       |
//  |        ↳EVAL𝑛  | O𝑛 population  | fit calc [1..n|max] | O𝑛         |
//  +----------------+----------------+---------------------+-----------+
//

template<class variant>
struct ea {
    
    char name[128] = "";
    
    ea_run run;
    ea_model model;
    ea_stats stats;
    
    // TODO: generalize collection type to store multiple objectives
    
    ea();
    
    void begin() {
        
        LOG(4, 0, 0, "\r\n--- BEGIN %s EA ---\r\n", this->name);
        
        this->model.isle.population.clear();
        this->model.isle.population.resize(config::mu_sub);
        this->stats.start = MPI_Wtime();
        
    }
    
    void end() {
        
        LOG(4, 0, 0, "\r\n--- END %s EA ---\r\n", this->name);
        
        double ea_end = MPI_Wtime();
        
        this->stats.local_t.value = ea_end - this->stats.start;
        this->stats.local_t.island = mpi.id;
        
        MPI_Reduce(&this->stats.local_t, &this->stats.min_t, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, this->model.tcomm);
        MPI_Reduce(&this->stats.local_t, &this->stats.max_t, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, this->model.tcomm);
        MPI_Reduce(&this->stats.local_t, &this->stats.sum_t, 1, MPI_DOUBLE, MPI_SUM, 0, this->model.tcomm);
        
        this->stats.avg_t = this->stats.sum_t / mpi.size;
        
    }
    
    #pragma mark EA::DATATYPE::FUNCTION::TEMPLATES
    
    template<typename f, typename o> void populate(objective<o> &obj, f function) { obj.populate(*this, function); }
    template<typename f, typename o, typename g> void evaluate(objective<o> &obj, g &individual, f function) { obj.evaluate(*this, individual, function); }
    template<typename f, typename o, typename m> void evolve(objective<o> &obj, m &meta, f function) { obj.evolve(*this, meta, function); }
    template<typename f, typename o, typename m, typename g> void evolve(objective<o> &obj, m &meta, g &individual, f function) { obj.evolve(*this, meta, individual, function); }
    
    template<typename e> void begin(ea<e> &target);
    template<typename genome> void begin(objective<genome> &obj, genome *current = NULL);
    template<typename genome> void end(objective<genome> &obj, genome *current = NULL);
    template<typename genome, typename i> void begin(objective<genome> &obj, i &interval = objective<genome>::run, genome *current = NULL);
    template<typename genome, typename i> void end(objective<genome> &obj, i &interval = objective<genome>::run, genome *current = NULL);
    
};

template<typename variant> template<typename genome> void ea<variant>::begin(objective<genome> &obj, genome *current) {
    
    this->begin();
    
    obj.init();
    
    this->begin(obj, obj.run, current);
    
}

template<typename variant> template<typename genome> void ea<variant>::end(objective<genome> &obj, genome *current) {
    
    this->end(obj, obj.run, obj.run.local);
    
}

template<typename variant> template<typename genome, typename i> void ea<variant>::begin(objective<genome> &obj, i &interval, genome *current) {

    obj.begin(interval, current);

    LOG(3, mpi.id, 0, "(%d,%d,%d,%d,%d) ea::%s::begin(%d,%d,%lu)\r\n", mpi.id, obj.id, obj.run.id, obj.run.cycle.id, obj.run.cycle.eval.id, this->name, obj.id, interval.id, sizeof(current));
    
}

template<typename variant> template<typename genome, typename i> void ea<variant>::end(objective<genome> &obj, i &interval, genome *current) {
    
    obj.end(interval, interval.local);
    
}

#endif /* dtype_ea_h */

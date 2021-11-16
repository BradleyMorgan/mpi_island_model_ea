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

// ----- start mpi derived datatypes...

template<typename e> void parallel_types(e &pea) {
    
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

    MPI_Type_create_struct(11, sol_lengths, sol_offsets, sol_types, &pea.variant.solution_type);
    MPI_Type_commit(&pea.variant.solution_type);

    MPI_Datatype visa_types[4] = { MPI_INT, MPI_INT, MPI_INT, MPI_CHAR };
    MPI_Aint visa_offsets[4] = { 0, sizeof(int), sizeof(int)*2, sizeof(int)*3 };

    int visa_lengths[4] = { 1, 1, 1, 64 };

    MPI_Type_create_struct(4, visa_lengths, visa_offsets, visa_types, &pea.variant.visa_type);
    MPI_Type_commit(&pea.variant.visa_type);

}

#pragma mark DATATYPE: @ea_variant{}

// describes the properties needed for parallel (island model) EA

struct ea_variant {

    unsigned long start;

    int islands;
    int island_size;
    int root = 0;

    // isle holds context specific population and migration data
    // from the perspective of a specific island instance (e.g. mpi rank)
    //
    // @island.h::island{} datatype
    //
    
    island isle;
    
    // TODO: maintaining a list of island_ids doesn't seem necessary
    
    std::vector<int> island_ids;
    
    MPI_Datatype solution_type;
    MPI_Datatype visa_type;
    
    MPI_Comm tcomm;
    
    template<typename e> void init(e &pea) {

        pea.variant.start = MPI_Wtime();
        
        // initialize MPI environment ...
        
        MPI_Init(NULL, NULL);
        MPI_Comm_size(MPI_COMM_WORLD, &pea.variant.islands);
        MPI_Comm_rank(MPI_COMM_WORLD, &pea.variant.isle.id);

        // load configuration items ...

        config::load("config.txt", pea.variant.islands, pea.variant.isle.id);

        for(int i=0; i < pea.variant.islands; i++) { pea.variant.island_ids.push_back(i); }

        pea.variant.tcomm = MPI_COMM_WORLD;
        pea.variant.island_size = config::island_mu;
        
        parallel_types(pea);
        
        pea.variant.isle.init();
        
    }
    
};

#pragma mark EA::DATATYPE: @ea{}

// wrapper type to serve as the trunk for ea properties

struct ea {
    
    ea_run run;
    ea_variant variant;
    
    double start = 0.0;
    double duration = 0.0;
    double init_start = 0.0;
    double init_duration = 0.0;
    
    // TODO: generalize collection type to store multiple objectives
    
    ea() {
        
//        this->start = MPI_Wtime();
//        this->init_start = MPI_Wtime();
//        this->run.id = 0;
//        this->run.eval.id = 0;
//        this->run.stats.init();
//        this->run.eval.stats.init();
        
    };
    
    #pragma mark EA::DATATYPE::FUNCTION::TEMPLATES
    
    template<typename o> void begin(ea_run &run, objective<o> &obj) { run.begin(); obj.begin(obj.run, *this); }
    template<typename o> void begin(ea_eval &eval, objective<o> &obj) { eval.begin(); obj.begin(obj.run.eval, *this); }
    
    //template<typename o> void begin(ea_run &run, objective<o> &obj, o &genome) { run.begin(); obj.run.begin(); obj.begin(obj.run, *this); }
    //template<typename o> void begin(ea_eval &eval, objective<o> &obj, o &genome) { eval.begin(); obj.begin(obj.run.eval, *this); }
    
    template<typename o> void end(ea_run &run, objective<o> &obj) { run.end(); obj.end(obj.run, *this); }
    template<typename o> void end(ea_eval &eval, objective<o> &obj) { eval.end(); obj.end(obj.run.eval, *this); }
    
//    template<typename o> void end(ea_run &run, objective<o> &obj, o &genome) { run.end(); }
//    template<typename o> void end(ea_eval &eval, objective<o> &obj, o &genome) { eval.end(); }
    
//    template<typename o, typename g> void log(ea_run &run, objective<o> &obj, ea &meta, g &genome) { log(run, obj, *this, meta, genome); }
//    template<typename o, typename g> void log(ea_eval &eval, objective<o> &obj, ea &meta, g &genome) { log(eval, obj, *this, meta, genome); }
//    
//    template<typename o, typename g> void log(objective_run &run, objective<o> &obj, ea &meta, g &genome) { obj.log(run, obj, *this, meta, genome); }
//    template<typename o, typename g> void log(objective_eval &eval, objective<o> &obj, ea &meta, g &genome) { obj.log(eval, obj, *this, meta, genome); }
    
    template<typename f, typename o> void populate(objective<o> &obj, f function) { obj.populate(*this, function); }
    
    template<typename f, typename o, typename g> void evaluate(objective<o> &obj, g &individual, f function) { obj.evaluate(*this, individual, function); }
    
    template<typename f, typename o, typename m> void evolve(objective<o> &obj, m &meta, f function) { obj.evolve(*this, meta, function); }
    template<typename f, typename o, typename m, typename g> void evolve(objective<o> &obj, m &meta, g &individual, f function) { obj.evolve(*this, meta, individual, function); }
    
    template<typename o> void end(objective<o> &obj) { obj.end(*this); }
    template<typename o> void end(ea &ea, objective<o> &obj) { ea_end(*this, obj); }
    
};

#endif /* dtype_ea_h */

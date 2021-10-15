//
//  main.cpp
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 2/26/20.
//  Copyright © 2020 Bradley Morgan. All rights reserved.
//

// EA -> Objective -> Population -> Genome

#include <mpi.h>
#include <algorithm>
#include <sys/time.h>
#include <ctime>
#include <string.h>
#include "config.h"
#include "utility.h"
#include "dtype_heirarchy.h"
#include "dtype_solution.h"
#include "dtype_island.h"
#include "dtype_topology.h"
#include "dtype_stats.h"
#include "objective.h"
#include "dtype_ea.h"
#include "ea.h"

ea solver_ea() {
    
    ea solver;
    
    solver.init(solver_init);
    
    return solver;
}

ea meta_ea(ea &variant) {
    
    ea meta;
    
    meta.init(meta_init, variant);
    
    return meta;
    
}

void solver_begin(ea &meta, ea &solver, topology &t) {

    for(solver.run.id = 1; solver.run.id <= solver.solutions.max_runs; solver.run.id++) {
    
        solver.begin(solver.run, solver.solutions);
        
        solver.populate(solver.solutions, solution_populate);
        solver.distribute(solver.solutions, solution_scatter);
        
        t.apply(solver.variant.isle, t);

        for(solver.run.eval.id = 1; solver.run.eval.id <= solver.solutions.max_evo_evals; solver.run.eval.id++) {

            solver.begin(solver.run.eval, solver.solutions);
            
            solver.evolve(solver.solutions, meta, t, solutions_evolve);
            
            solver.end(solver.run.eval, solver.solutions);

        }

        t.fitness = (t.total_migration_time / t.rounds) * -1;
        
        solver.end(solver.run, solver.solutions);
        
    }

}

void meta_begin(ea &meta, ea &solver) {
    
    for(meta.run.id = 1; meta.run.id <= meta.topologies.max_runs; meta.run.id++) {
        
        meta.begin(meta.run, meta.topologies);
        
        meta.populate(meta.topologies, topologies_populate);
        
        for(int i=0; i<meta.topologies.mu; i++) {
            
            meta.begin(meta.run.eval, meta.topologies);
            
            solver_begin(meta, solver, meta.topologies.population[i]);
            
            meta.end(meta.run.eval, meta.topologies);

        }

        for(meta.run.eval.id = 1; meta.run.eval.id <= meta.topologies.max_evo_evals; meta.run.eval.id++) {

            solver.evolve(meta.topologies, meta, topology_evolve);

        }
        
        meta.end(meta.run, meta.topologies);
        
    }
    
}

//meta_ea()
//
//    topology run [1 .. meta_ea::runs]
//
//        generate topology population
//
//        for each topology
//
//            until meta::fitness_termination
//
//                fitness = solver_ea(topology)
//
//        until meta::evo_termination
//
//            create O offspring
//
//            for each offspring
//
//                until meta::fitness_termination
//
//                    fitness = solver_ea(topology)
//
//            select survivors
//
//
//solver_ea(topology)
//
//    solver run [1 .. solver_ea::runs]
//
//        generate solution population
//
//        apply topology
//
//        until max_solution_evo_termination [1000]
//
//            create O offspring
//            evolve S solution generations [1..1000]
//
//        return T fitness


#pragma mark FUNCTION: main()

// core initialization, overarching EA logic:
// this implementation is a parallelized (island model) multi-objective ea
// evolving solutions for a minimization problem (offset rastrigin, objective 1)
// and a population of communication topologies, measured in duration (objective 2)
//
// note: functional dependency and potential fitness influence of objective 1 as it is
// applied to the communication (migration) pattern employed by the island model objective 2
//
// note: only objective (1) employs island model (parallel) evolution
//

int main(int argc, const char * argv[]) {
    
    ea solver = solver_ea();
    ea meta = meta_ea(solver);
    
    if(config::ea_mode > 0) {
        meta_begin(meta, solver);
    } else {
        benchmark_topology(meta);
        solver_begin(meta, solver, meta.topologies.population[0]);
    }
    
    solver.end(solver.solutions);
    meta.end(meta.topologies);
    
    
}

//
//  main.cpp
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 2/26/20.
//  Copyright © 2020 Bradley Morgan. All rights reserved.
//

#include <mpi.h>
#include <algorithm>
#include "config.h"
#include "utility.h"
#include "dtype_heirarchy.h"
#include "dtype_solution.h"
#include "dtype_island.h"
#include "dtype_topology.h"
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

    for(solver.solutions.run = 1; solver.solutions.run <= solver.solutions.max_runs; solver.solutions.run++) {
        
        solver.populate(solver.solutions, solution_populate);
        solver.distribute(solver.solutions, solution_scatter);
        
        t.apply(solver.variant.isle, t);

        for(solver.solutions.eval = 1; solver.solutions.eval <= solver.solutions.max_evo_evals; solver.solutions.eval++) {

            solver.evolve(meta, solver.solutions, solutions_evolve, t);

        }
//
//        t.fitness = (t.total_migration_time / t.rounds) * -1;
        
    }

}

void meta_begin(ea &meta, ea &solver) {
    
    for(meta.topologies.run = 1; meta.topologies.run <= meta.topologies.max_runs; meta.topologies.run++) {
        
        meta.run_init(meta_run_init);
        
        meta.populate(meta.topologies, topologies_populate);
        
        for(int i=0; i<meta.topologies.mu; i++) {
            
            printf("TOPOLOGY %d\r\n", meta.topologies.population[i].id);
            solver_begin(meta, solver, meta.topologies.population[i]);

        }

//        while(meta.topologies.eval <= meta.topologies.max_evo_evals) {
//
//            solver.solutions.evolve(solver, topology_evolve);
//
//        }
        
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
    
    meta_begin(meta, solver);
    
    
    
    // perform mpi and global variable initialization ...
    
//    ea multi = ea_init();
//    multi.topologies.run = 1;
//
//    // begin experimental runs ...
//
//    for(multi.topologies.run = 1; multi.topologies.run <= multi.topologies.max_runs; multi.topologies.run++) {  // n (param) experimental runs
//
//        for(multi.solutions.run = 1; multi.solutions.run <= multi.solutions.max_runs; multi.solutions.run++) {  // n (param) experimental runs
//
//            // clear any current populations and reset counters ...
//
//            multi.run_init();
//
//            // generate random objective (solution) population with fitness calculated
//            // and distribute (mu/n) subpopulations to each island
//
//            multi.populate(multi.solutions, solution_populate);
//            multi.distribute(multi.solutions, solution_scatter);
//
//            if(config::ea_mode > 0) {  // multi-objective, evolve topologies
//
//                if(multi.topologies.run == 1) {
//                    multi.populate(multi.topologies, topologies_populate);
//                }
//
//                // evaluate initial topology population by applying each topology and using it for n solution evals
//
//                for(int i=0; i<multi.topologies.mu; i++) {
//
//                    multi.evaluate(multi.topologies, multi.topologies.population[i], topology_evaluate);
//
//                    if(multi.solutions.eval >= multi.solutions.max_evo_evals) { break; }
//
//                }
//
//                if(multi.solutions.eval >= multi.solutions.max_evo_evals) { continue; }
//
//                // perform the remaining solution evolution cycles indirectly through topology evolution
//                // topology evolution also triggers evolutionary cycles of the solution population in order
//                // to determine topology fitness
//
//                while(multi.solutions.eval <= multi.solutions.max_evo_evals) {
//                    multi.evolve(multi.solutions, topology_evolve);
//                }
//
//            } else {  // benchmark, use static ring topology
//
//                // in the benchmark case, we only need to apply the topology once ...
//
//                if(multi.solutions.run == 1) {
//                    benchmark_topology(multi);
//                    multi.topologies.population[0].apply(multi.variant.isle, multi.topologies.population[0]);
//                }
//
//                // perform solver evolution ...
//
//                while(multi.solutions.eval <= multi.solutions.max_evo_evals) {
//                    multi.coevolve(multi.solutions, solutions_evolve, multi.topologies.population[0]);
//                }
//
//            }
//
//            // perform any value resets and write the run statistics to file output ...
//
//            multi.run_end();
//
//        } // ea 1 run end
//
//    } // ea 2 run_end
//
//    multi.ea_end();
//
//
////    logic test ...
//
////    int o1run = 1;
////    int o2run = 1;
////
////    while(o1run <= config::objective_1_runs && o2run <= config::objective_2_runs) {
////
////        // objective 1 (o1): topology
////        // generate initial topology population
////
////        // objective 2 (o2): solver
////        // generate initial solution population
////
////        // T1:run=1
////
////        int t = 0; // topology index
////
////        int o1_fit_eval = 1;
////        int o1eval = 1;
////        int o2eval = 1;
////
////        // (o1) apply topology t
////
////        while(t <= config::topo_mu) {
////
////            // (o1) calculate t.fitness as migration time from (o2)
////            // (o2) evolution cycle
////
////            // T2:t=7
////
////            for(int i=0; i<=config::lambda; i++) {
////                // (o2) parent selection
////                // (o2) crossover
////                // (o2) mutation
////            }
////
////            // (o2) survival selection
////            // (o2) migration
////            // (o1) fitness += migration time
////
////            o1_fit_eval++;
////            o2eval++;
////
////            if(o1_fit_eval%config::objective_1_max_fit_evals == 0) {
////                t++;
////                // (o1) apply topology t
////            }
////
////            if(o2eval >= config::objective_2_max_evo_evals) {
////                o2run++;
////                // (o2) init run
////            }
////
////            // T3:o1eval=100
////            // T4:o2eval=200
////
////        }
////
////        // (o1) initial population fitness calculated
////
////        while (o1eval <= config::objective_1_max_evo_evals) {
////
////            // (o1) evolution cycle
////
////            for(int i=0; i<=config::topo_lambda; i++) {
////                // (o1) parent selection
////                // (o1) crossover
////                // (o1) mutation
////            }
////
////            // (o1) calculate child fitness -> (o2)
////            // (o1) survival selection
////
////            o1eval++;
////
////        }
////
////        o2run++;
////
////    }
    
}

//
//  ea.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 3/30/20.
//  Copyright © 2020 Bradley Morgan. All rights reserved.
//

#ifndef ea_h
#define ea_h

#include <sys/time.h>
#include "solution.h"
#include "island.h"
#include "topology.h"
#include "objective.h"
#include "stats.h"

template<typename genome> genome parent(objective<genome> &o) {
    
    genome p;

    o.cpd();
    
    // implementation uses the single armed roulette wheel approach to select
    // an individual from the population
    
    int i = 1;
    
    // random double uniformly distributed between 0 and 1
    
    double r = ((double)rand()/(double)RAND_MAX);
    
    // spin the wheel
    
    while (o.aggregate.value.cpd[i] < r ) { i++; }
    
    p = o.population[i];
    
    return p;
        
}

// pick a parent from the population using the provided method ...

solution select_parent(island &isle) {

    // implementation uses the single armed roulette wheel approach to select
    // an individual from the population

    int i = 1;

    // random double uniformly distributed between 0 and 1

    LOG(8, 0, 0, "random double uniformly distributed between 0 and 1\r\n");

    double r = ((double)rand()/(double)RAND_MAX);

    // spin the wheel

    LOG(10, 0, 0, "spin the wheel ...\r\n");

    while (isle.cpd[i] < r ) { LOG(10, 0, 0, "island %d, cpd[%d] = %f, r = %f\r\n", isle.id, i, isle.cpd[i], r); i++; }

    LOG(8, 0, 0, "island %d selected individual %d cpd[%d] = %f = ...\r\n", isle.id, i, i, isle.cpd[i]);
    
    isle.population[i].selected++;
    
    LOG(5, 0, 0, "sol %s selected=%d", isle.population[i].id, isle.population[i].selected);
    
    solution s = isle.population[i];
    
    return s;

}


void select_survivors(island &isle, std::vector<solution> &children, int island_mu) {

    isle.calculator.average_fitness(isle);
    
    LOG(6, 0, 0, "island %d survival, population before: %lu, fitness before: %f, best fit: %f\r\n", isle.id, isle.population.size(), isle.average_fitness, isle.population[0].fitness);
    
    // truncation: add new children to the population, and then kill the weakest
    
    isle.population.insert(isle.population.end(), children.begin(), children.end());
    
    std::sort(isle.population.begin(), isle.population.end(), compare_fitness<solution>);
    std::reverse(isle.population.begin(), isle.population.end());
    
    isle.population.erase(isle.population.begin()+island_mu, isle.population.end());
    
    for_each(isle.population.begin(), isle.population.end(), [](solution &g) { g.survival++; });
    
    isle.calculator.average_fitness(isle);
    
    LOG(6, 0, 0, "island %d survival, population after: %lu, fitness after: %f, best fit: %f\r\n", isle.id, isle.population.size(), isle.average_fitness, isle.population[0].fitness);
    
}


void mutate(solution &mutant) {

    for(int i=0; i<DIM; i++){
        
        if(rand()/(RAND_MAX+1.0) < config::ea_1_mutation_rate) {
     
            mutant.input[i] = drand(-5.12, 5.12);
            
        }
        
    }
    
}

std::vector<solution> crossover(ea &multi) {
      
    multi.variant.isle.calculator.cpd(multi.variant.isle);
    
    std::vector<solution> children;
    
    for(int i = 0; i < config::island_lambda; i++) {
        
        LOG(5, 0, 0, "island %d creating objective<solution> child %d ...\r\n", multi.variant.isle.id, i);
        
        solution p1 = select_parent(multi.variant.isle);
        LOG(6, 0, 0, "island %d selected p1<solution> cpd = %f, fitness = %f ...\r\n", multi.variant.isle.id, p1.selection_distribution, p1.fitness);
        solution p2 = select_parent(multi.variant.isle);
        LOG(6, 0, 0, "island %d selected p2<solution> cpd = %f, fitness = %f...\r\n", multi.variant.isle.id, p2.selection_distribution, p2.fitness);
        
        solution child;
        
        strcpy(child.id, uniqid(sinstances++));
        
        for(int j=0; j<DIM; j++) {
            if(rand()%2 == 1) {
                int gene = rand()%DIM;
                child.input[j] = p1.input[gene];
                child.group += p1.input[gene];
                LOG(7, 0, 0, "assigning gene %d = %f from p1<solution> to child %d (%d) ...\r\n", j, child.input[j], i, gene);
            } else {
                int gene = rand()%DIM;
                child.input[j] = p2.input[gene];
                child.group += p2.input[gene];
                LOG(7, 0, 0, "assigning gene %d = %f from p2<solution> to child %d (%d) ...\r\n", j, child.input[j], i, gene);
            }
        }
        
        if(rand()/(RAND_MAX+1.0) < multi.solutions.mutation_rate) {
            
            LOG(6, 0, 0, "island %d mutating child<solution> %d ...\r\n", multi.variant.isle.id, i);
            
            mutate(child);
            
        }
        
        child.fitness = offset_rastrigin(child.input, multi.offsets);
        
        if(child.fitness > (p1.fitness / 4)) {
            LOG(4, 0, 0, "LOW island %d %f<->%f child<solution> %d fitness %f\r\n", multi.variant.isle.id, p1.fitness, p2.fitness, i, child.fitness);
        } else {
            LOG(8, 0, 0, "child<solution> %d fitness %f\r\n", i, child.fitness);
        }
        
        child.source = multi.variant.isle.id;
        child.locale = multi.variant.isle.id;
        child.migrations = 0;
        
        strcpy(child.parents[0], p1.id);
        strcpy(child.parents[1], p2.id);
        
        LOG(5, 0, 0, "child id %s p1 %s p2 %s\r\n", child.id, child.parents[0], child.parents[1]);
        
        children.push_back(child);
        
    }
    
    LOG(6, 0, 0, "rank %d returning %lu child solutions from crossover\r\n", multi.variant.isle.id, children.size());
    
    return children;
    
}

#pragma mark FUNCTION: solution_populate()

// returns a vector of a randomly generated offset rastrigin solution population                        |
// accepts an array of floating point numbers of n-dimensions @offsets[param:0]                         |
// used to calculate the solution fitness.                                                              |

void solution_populate(ea &solver) {
    
    LOG(6, 0, 0, "rank %d of %d entered solution_populate\r\n", solver.variant.isle.id, solver.variant.islands);
    
    if(solver.variant.isle.id != 0) {
        LOG(6, 0, 0, "rank %d leaving solution_populate\r\n", solver.variant.isle.id);
        return;
    }
    
    // the "solution" datatype represents a single offset rastrigin solution as
    // an array of size DIM = @dim[config.txt:12] holding the solution's randomly
    // generated gene values.
    
    LOG(4, solver.variant.isle.id, 0, "initializing objective (solution) population ...\r\n");
    LOG(6, 0, 0, "island %d (root) initializing mu=%d solutions ...\r\n", solver.variant.isle.id, solver.solutions.mu);

    for(int i=0; i<solver.solutions.mu; i++) {
        
        LOG(6, 0, 0, "creating solution %d ...\r\n", i);
        
        solution p;
        
        //p.id = uniqid(sinstances++);
        //strcpy(p.id, uniqid(sinstances++));
        
        p.source = solver.variant.isle.id;
        p.migrations++;
        
        //p.parents[0] = 0;
        //p.parents[1] = 0;
        
        strcpy(p.parents[0], "0");
        strcpy(p.parents[1], "0");
        
        LOG(6, 0, 0, "rank %d assigning solution values\r\n", solver.variant.isle.id);
        
        for (int j = 0; j < DIM; j++) {
            p.input[j] = drand(-5.12, 5.12); // rastrigin says: x[i] ∈ [-5.12,5.12]
            p.group += p.input[j];
            LOG(10, solver.variant.isle.id, 0, "rank %d solution %d group = %f\r\n", solver.variant.isle.id, i, p.group);
        }
        
        LOG(10, solver.variant.isle.id, 0, "rank %d solution %d FINAL id = %f\r\n", solver.variant.isle.id, i, p.group);
                
        LOG(6, 0, 0, "rank %d assigning solution fitness\r\n", solver.variant.isle.id);
        
        p.fitness = offset_rastrigin(p.input, solver.offsets);
        
        LOG(6, 0, 0, "rank %d assigned solution fitness = %f\r\n", solver.variant.isle.id, p.fitness);
        
        LOG(6, 0, 0, "rank %d adding solution %d with fitness %f to population, current size = %lu\r\n", solver.variant.isle.id, i, p.fitness, solver.solutions.population.size());
        
        solver.solutions.aggregate.value.fitness += p.fitness;

        visa v(solver.run.eval.id, p.source, p.source, p.id);

        //multi.variant.isle.visas.push_back(v);

        solver.solutions.population.push_back(p);
    
        LOG(6, 0, 0, "island %d (root) initialized solution %d with fitness %f ...\r\n", solver.variant.isle.id, i, p.fitness);
        
    }
    
    LOG(4, solver.variant.isle.id, 0, "initialized objective (solution) population: total fitness = %f\r\n", solver.solutions.aggregate.value.fitness);
    LOG(6, 0, 0, "%lu solutions initialized ...\r\n", solver.solutions.population.size());

//    std::sort(multi.solutions.population.begin(), multi.solutions.population.end(), compare_fitness);
//    std::reverse(multi.solutions.population.begin(), multi.solutions.population.end());
    
    if(solver.run.eval.stats.global_best_fitness == 0.0) {
        solver.run.eval.stats.global_best_fitness = solver.solutions.population.data()[0].fitness;
    }
    
    
    LOG(6, 0, 0, "rank %d leaving solution_populate\r\n", solver.variant.isle.id);
    
    // we have our initial primary population with the calculated fitnesses
    
}

#pragma mark FUNCTION: solution_scatter()

// separate the single full population from the root process to subpopulations across all processes ...

void solution_scatter(ea &solver) {
    
    LOG(6, 0, 0, "rank %d entered solution_scatter\r\n", solver.variant.isle.id);

    LOG(6, 0, 0, "rank %d of %d resized local island population to %lu\r\n", solver.variant.isle.id, solver.variant.islands, solver.variant.isle.population.size());
    
    if(solver.variant.isle.id == 0) {
        LOG(4, solver.variant.isle.id, 0, "rank 0 scattering population root size = %lu mem 0 = %f...\r\n", solver.solutions.population.size(), solver.solutions.population[0].fitness);
    } else {
          LOG(4, 0, 0, "rank %d instantiating scatter with = %lu subpopulation size ...\r\n", solver.variant.isle.id, solver.variant.isle.population.size());
    }

    double scatter_start = MPI_Wtime();
    
    MPI_Scatter(&solver.solutions.population[0], config::island_mu, solver.variant.solution_type, &solver.variant.isle.population[0], config::island_mu, solver.variant.solution_type, 0, solver.variant.tcomm);
    
    LOG(4, solver.variant.isle.id, 0, "rank %d scatter return ...\r\n", solver.variant.isle.id);
    
    double scatter_end = MPI_Wtime();
    double scatter_time = scatter_end - scatter_start;

    LOG(4, 0, 0, "population scattered, island %d population %lu mem 0 = %f...\r\n", solver.variant.isle.id, solver.variant.isle.population.size(), solver.variant.isle.population[0].fitness);

    solver.variant.isle.calculator.average_fitness(solver.variant.isle);

    LOG(4, 0, 0, "island %d scattered population size %lu, average fitness: %f\r\n", solver.variant.isle.id, solver.variant.isle.population.size(), solver.variant.isle.average_fitness);

    MPI_Reduce(&scatter_time, &solver.run.eval.stats.total_scatter_time, 1, MPI_DOUBLE, MPI_SUM, 0, solver.variant.tcomm);

    solver.run.eval.stats.total_scatter_time += scatter_time;

    LOG(6, 0, 0, "rank %d leaving solution_scatter\r\n", solver.variant.isle.id);
    
}

#pragma mark FUNCTION: solution_eval()

void solutions_evolve(ea &solver, ea &meta, topology &t) {
    
    LOG(5, solver.variant.isle.id, 0, "beginning objective<solution> evaluation %d, topology %d\r\n", solver.run.eval.id, t.id);
    
    LOG(8, solver.variant.isle.id, 0, "calculating island %d cpd, topology %d, eval %d\r\n", solver.variant.isle.id, t.id, solver.run.eval.id);
    
    solver.solutions.fitness();
    solver.solutions.cpd();
    
    LOG(6, solver.variant.isle.id, 0, "performing crossover island %d, topology %d, eval %d\r\n", solver.variant.isle.id, t.id, solver.run.eval.id);

    // crossover function performs parent selection iteratively, creating n child solutions with fitness calculated ...
    
    std::vector<solution> children = crossover(solver);

    LOG(4, solver.variant.isle.id, 0, "island %d created %lu child solutions, population size = %lu ... selecting %d survivors  ...\r\n", solver.variant.isle.id, children.size(), solver.variant.isle.population.size(), solver.variant.island_size);

    // select the best solutions and remove n children ...
    
    select_survivors(solver.variant.isle, children, solver.variant.island_size);

    LOG(6, solver.variant.isle.id, 0, "island migrations ...\r\n");
    
    // perform migration of best solutions using the currently applied topology ...
    
    if(solver.run.eval.id%config::migration_interval == 0) {
        
        double migrate_start = MPI_Wtime();
        
        // issue migration imports and exports ...
        
        island::migration::send(solver.variant.isle, solver.variant.solution_type, solver.run.eval.id);
        island::migration::receive(solver.variant.isle, solver.variant.solution_type, solver.run.eval.id);
        
        double migrate_end = MPI_Wtime();
        double migrate_time = migrate_end - migrate_start;
        
        LOG(8, 0, 0, "migrate start = %3.10f, migrate end = %3.10f, migrate time = %3.10f\r\n", migrate_start, migrate_end, migrate_time);
        
        solver.run.eval.stats.total_migrate_time += migrate_time;
     
        // aggregate the migration time from each island into the topology fitness ...
        
        MPI_Reduce(&migrate_time, &t.round_fitness, 1, MPI_DOUBLE, MPI_SUM, 0, solver.variant.tcomm);
        
    }

    LOG(4, solver.variant.isle.id, 0, "rank %d (root) topology %d, eval %d, round %d, tfitness = %f\r\n", solver.variant.isle.id, t.id, solver.run.eval.id, t.rounds, t.fitness);
    
    LOG(6, 0, 0, "rank %d topology = %d rounds = %d eval = %d\r\n", solver.variant.isle.id, t.id, t.rounds, meta.run.eval.id);

    if(solver.variant.isle.id == 0) {
        
        t.rounds++;
        t.total_migration_time += t.round_fitness;
        t.round_fitness = t.round_fitness * -1;
        t.fitness = (t.total_migration_time / t.rounds) * -1;
        
        meta.run.eval.stats.topo_migrate_time += t.total_migration_time;
        
    }
    
    LOG(8, solver.variant.isle.id, 0, "rounds=%d, total_time=%013.10f, round_fit=%013.10f, fit=%013.10f\r\n", t.rounds, t.total_migration_time, t.round_fitness, t.fitness);
    
    // output for various intervals ...
    
//    if(multi.solutions.eval%10 == 0) {
//        LOG(6, 0, 0, "island %d average fit %f\r\n", multi.meta.isle.id, multi.meta.isle.average_fitness);
//    }
    
    if(solver.run.eval.id%config::ea_1_log_interval == 0) {
    
        LOG(6, 0, 0, "gathering population, subpopulation %d size %lu avg fitness %f...\r\n", solver.variant.isle.id, solver.variant.isle.population.size(), solver.variant.isle.average_fitness);
        
        double gather_start = MPI_Wtime();
        
        // gather island subpopulations back into the aggregate population on rank 0 ...
        
        MPI_Gather(&solver.variant.isle.population[0], solver.variant.island_size, solver.variant.solution_type, &solver.solutions.population[0], solver.variant.island_size, solver.variant.solution_type, 0, solver.variant.tcomm);
        
        double gather_end = MPI_Wtime();
        double gather_time = gather_end - gather_start;
                    
        LOG(4, solver.variant.isle.id, 0, "population gathered by rank %d, aggregate size %lu, subpopulation size = %lu ...\r\n", solver.variant.isle.id, solver.solutions.population.size(), solver.variant.isle.population.size());
                    
        solver.run.eval.stats.total_gather_time += gather_time;
        
    }
    
    if(solver.run.eval.id%config::ea_1_log_interval == 0 && solver.variant.isle.id == 0) {
        
        LOG(4, 0, 0, "population size %lu, member = %2.10f\r\n", solver.variant.isle.population.size(), solver.variant.isle.population[0].fitness);
    
        if(solver.variant.isle.id == 0) {
            
            log_fn_eval_stats(solver, meta, t);
            
        }
        
    }
    
    if(solver.run.eval.id%config::ea_1_max_fit_evals == 0) {
        
        LOG(4, 0, 0, "population size %lu, member = %2.10f\r\n", solver.variant.isle.population.size(), solver.variant.isle.population[0].fitness);
        
        if(solver.variant.isle.id == 0) {
            log_fn_topology_stats(solver, meta, t);
        }
        
    }
    
    if(solver.run.eval.id%config::ea_1_population_log_interval == 0) {
        
        log_pop_stats(solver, solver.variant.isle, solver.variant.visa_type);
        
    }

}

void solver_init(ea &solver) {
    
    solver.run.id = 1;
    solver.run.eval.id = 1;
    solver.run.stats.init();
    solver.run.eval.stats.init();
    
    solver.start = MPI_Wtime();
    solver.init_start = MPI_Wtime();
    solver.variant.start = MPI_Wtime();
    
    // initialize MPI environment ...

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &solver.variant.islands);
    MPI_Comm_rank(MPI_COMM_WORLD, &solver.variant.isle.id);

    // load configuration items ...

    config::load("config.txt", solver.variant.islands, solver.variant.isle.id);
    
    for(int i=0; i<solver.variant.islands; i++) { solver.variant.island_ids.push_back(i); }

    // ----- start mpi derived datatypes...

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
    
    MPI_Type_create_struct(11, sol_lengths, sol_offsets, sol_types, &solver.variant.solution_type);
    MPI_Type_commit(&solver.variant.solution_type);
    
    MPI_Datatype visa_types[4] = { MPI_INT, MPI_INT, MPI_INT, MPI_CHAR };
    MPI_Aint visa_offsets[4] = { 0, sizeof(int), sizeof(int)*2, sizeof(int)*3 };
    
    int visa_lengths[4] = { 1, 1, 1, 64 };
    
    MPI_Type_create_struct(4, visa_lengths, visa_offsets, visa_types, &solver.variant.visa_type);
    MPI_Type_commit(&solver.variant.visa_type);
    
    // ----- end derived mpi datatypes
    
    solver.variant.tcomm = MPI_COMM_WORLD;
    solver.variant.island_size = config::island_mu;
    solver.variant.isle.init();
    
    solver.solutions.mu = config::ea_1_mu;
    solver.solutions.lambda = config::ea_1_lambda;
    solver.solutions.max_runs = config::ea_1_runs;
    solver.solutions.mutation_rate = config::ea_1_mutation_rate;
    solver.solutions.max_evo_evals = config::ea_1_max_evo_evals;
    solver.solutions.max_fit_evals = config::ea_1_max_fit_evals;
    
    if(solver.variant.isle.id == 0) {
        solver.offsets = generate_offsets(-2.5, 2.5, .5);
    }
    
    MPI_Bcast(&solver.offsets, DIM, MPI_DOUBLE, 0, solver.variant.tcomm);
    
    // collect the time consumed by all islands in this initialization ...

    // TODO: the init time gather segfaults on higher core count runs, so may need to debug at some point, but currently the init duration is somewhat insigificant
    // MPI_Gather(&local_init_duration, 1, MPI_DOUBLE, &multi.run.stats.init_duration, 1, MPI_DOUBLE, 0, multi.meta.isle.tcomm);
    
    LOG(2, solver.variant.isle.id, 0, "world size: %d\r\n", solver.variant.islands);
    LOG(2, solver.variant.isle.id, 0, "mu mode: %d\r\n", config::mu_mode);
    LOG(2, solver.variant.isle.id, 0, "global mu %s %d\r\n", config::mu_msg, config::ea_1_mu);
    LOG(2, solver.variant.isle.id, 0, "island mu %s %d\r\n", config::subpop_msg, config::island_mu);
    LOG(2, solver.variant.isle.id, 0, "island lambda %s %d\r\n", config::lambda_msg, stoi(config::items["island_lambda"]));
    
    solver.init_duration += (MPI_Wtime() - solver.init_start);
    
}

void meta_init(ea &meta, ea &solver) {
    
    meta.variant = solver.variant;
    
    meta.start = MPI_Wtime();
    meta.init_start = MPI_Wtime();
    meta.variant.start = MPI_Wtime();
    
    meta.run.id = 1;
    meta.run.eval.id = 1;
    meta.run.stats.init();
    meta.run.eval.stats.init();
    
    meta.topologies.mu = config::ea_2_mu;
    meta.topologies.lambda = config::ea_2_lambda;
    meta.topologies.max_runs = config::ea_2_runs;
    meta.topologies.mutation_rate = config::ea_2_mutation_rate;
    meta.topologies.max_evo_evals = config::ea_2_max_evo_evals;
    meta.topologies.max_fit_evals = config::ea_2_max_fit_evals;
    meta.topologies.islands = meta.variant.islands;
    
    meta.init_duration += (MPI_Wtime() - meta.init_start);
    
}






#pragma mark FUNCTION: topology_crossover()

std::vector<topology> topology_crossover(ea &meta) {
   
    std::vector<topology> children;
    
    if(meta.variant.isle.id != 0) { children.resize(meta.topologies.lambda); return children; }
        
    meta.topologies.cpd();
    
    for(int n = 0; n < meta.topologies.lambda; n++) { // loop lambda

        // create child skeleton ...
        
        topology child;
        
        child.world_size = meta.variant.islands;
        child.fitness = 0.0;
        child.channel_count = 0;
        child.round_fitness = 0.0;
        child.selection_distribution = 0.0;
        child.channels.resize(meta.variant.islands);
        child.channels.clear();
        
        LOG(3, meta.variant.isle.id, 0, "creating topo kids\r\n");

        topology t1 = meta.topologies.select(parent<topology>);
        topology t2 = meta.topologies.select(parent<topology>);

        // calculate an adjacency matrix for each parent's associated topology for use in
        // generating child topology ...

        LOG(3, meta.variant.isle.id, 0, "parents<topology> t1=%2.10f,t2=%2.10f ...\r\n", t1.fitness, t2.fitness);

        std::vector<std::vector<int>> m1 = topology::create::matrix(t1);
        std::vector<std::vector<int>> m2 = topology::create::matrix(t2);

        // recombine the parent adjacency matrices, initialize ...

        LOG(3, meta.variant.isle.id, 0, "recombining topology %d <-> %d ...\r\n", t1.id, t2.id);

        std::vector<std::vector<int>> child_matrix;
        child_matrix.resize(meta.variant.islands);

        int comm_count = 0;
        int rec_count[meta.variant.islands];
        int snd_count[meta.variant.islands];

        for(int i=0; i<meta.variant.islands; i++) {
            rec_count[i] = 0;
            snd_count[i] = 0;
        }

        LOG(3, meta.variant.isle.id, 0, "child<topology> %d initialized \r\n", child.id);

        // iterate row->column for each x,y element in the child matrix, and for each
        // gene and randomly choose a parent from which to assign the value ...

        LOG(3, meta.variant.isle.id, 0, "performing child<topology> %d matrix crossover ...\r\n", child.id);

        while(comm_count == 0) { // failsafe to prevent empty matrix

            for(int i=0; i<m1.size(); i++) {  // child matrix row

                child_matrix[i].resize(meta.variant.islands);

                for(int j=0; j<m2.size(); j++) { // child matrix column

                    if(rec_count[j] >= config::migration_cap) {
                        LOG(6, 0, 0, "migration cap limit reached for process %d\r\n", i);
                        continue;
                    }

                    if(snd_count[i] >= config::send_cap) {
                        LOG(6, 0, 0, "send cap limit reached for process %d\r\n", i);
                        continue;
                    }

                    // we don't want an island sending migrants to itself
                    // also consider sparsity parameter ...
                    
                    if(i != j) {

                        if(rand()%2 == 1) { // coin flip, heads take the row index value from parent 1
                            LOG(10, 0, 0, "assigning child<topology> m1[%d][%d] -> %d\r\n", i, j, m1[i][j]);
                            
                            if(m1[i][j] == 1 && rand()/(RAND_MAX+1.0) > config::sparsity) {
                                child_matrix[i][j] = m1[i][j];
                                rec_count[j]++;
                                snd_count[i]++;
                            } else {
                                child_matrix[i][j] = m1[i][j];
                            }
                            
                        } else { // tails, take it from parent 2
                            LOG(10, 0, 0, "assigning child<topology> m2[%d][%d] -> %d\r\n", i, j, m2[i][j]);

                            if(m2[i][j] == 1 && rand()/(RAND_MAX+1.0) > config::sparsity) {
                                child_matrix[i][j] = m2[i][j];
                                rec_count[j]++;
                                snd_count[i]++;
                            } else {
                                child_matrix[i][j] = m2[i][j];
                            }
                        }

                        if(child_matrix[i][j] == 1) { comm_count++; }

                    } else {

                        child_matrix[i][j] = 0;

                    }

                } // child matrix column

                LOG(10, 0, 0, "------\r\n");

            }  // child matrix row

        }  // end failsafe
        
        // mutation, interate through the matrix and flip the bit with probability m,
        // also considering any sparsity constraints ...
        
        for(int i=0; i<child_matrix.size(); i++) { // matrix row

            for(int j=0; j<child_matrix[i].size(); j++) {  // matrix col

                if(rec_count[j] >= config::migration_cap) {
                    LOG(6, 0, 0, "migration cap limit reached for process %d\r\n", i);
                    continue;
                }

                if(snd_count[i] >= config::send_cap) {
                    LOG(6, 0, 0, "send cap limit reached for process %d\r\n", i);
                    continue;
                }

                if(rand()/(RAND_MAX+1.0) < meta.topologies.mutation_rate) {

                    LOG(3, meta.variant.isle.id, 0, "mutating child<topology> %d ...\r\n", child.id);

                    if(child_matrix[i][j] == 0 && rand()/(RAND_MAX+1.0) < config::sparsity) {

                        child_matrix[i][j] = 1;
                        rec_count[j]++;
                        snd_count[i]++;
                        comm_count++;

                    }

                    if(child_matrix[i][j] == 1 && rand()/(RAND_MAX+1.0) > config::sparsity && comm_count > 1) {

                        child_matrix[i][j] = 0;
                        comm_count--;

                    }

                }

            } // matrix col

        } // matrix row
        
        // convert the adjaceny matrix to a sender and receiver arrays for use with MPI send\recv ...
    
        topology::create::channels(child, child_matrix);

        LOG(3, meta.variant.isle.id, 0, "child<topology> %d born!\r\n", child.id);

        children.push_back(child);

        LOG(3, 0, 0, "child<topology> %d created by rank %d from matrix with %lu senders and %lu receivers\r\n", child.id, meta.variant.isle.id, child.channels[n].senders.size(), child.channels[n].receivers.size());

    } // loop lambda
    
    return children;
    
}

#pragma mark FUNCTION: topology_evolve()

void topology_evolve(ea &solver, ea &meta) {
    
    std::vector<topology> children;
    
    if(meta.variant.isle.id == 0) {
        
        // the root island holds the only valid, authoritative topology population,
        // we only need to perform parent selection and recombination on that process ...
        
        children = topology_crossover(meta);
        
        LOG(6, 0, 0, "island %d topology survival, population before: %lu, fitness before: %f, best fit: %f\r\n", meta.variant.isle.id, meta.topologies.population.size(), meta.topologies.aggregate.value.fitness, meta.topologies.population[0].fitness);
            
    } else {
        
        // ensure that the secondary islands have memory allocated to avoid potential segfault ....
        
        children.resize(meta.topologies.lambda);
        
    }

    LOG(6, 0, 0, "island %d <topology> survival, population after: %lu, fitness after: %f, best fit: %f\r\n", meta.variant.isle.id, meta.topologies.population.size(), meta.topologies.aggregate.value.fitness, meta.topologies.population[0].fitness);

    // in order to determine a topoology fitness, perform evolutionary cycles of the solver population,
    // with all islands participating using the send\recv channels as defined in the topology
    // and then assign the aggregate time spent in migration as topology fitness ...
    
    for (std::vector<topology>::iterator it = children.begin(); it != children.end(); ++it) {

        LOG(4, meta.variant.isle.id, 0, "rank %d (root) applying child<topology>[%d] for eval\r\n", meta.variant.isle.id, it->id);

        // the apply function distributes a topology's send\recv channels to all islands ...

        it->apply(meta.variant.isle, *it);

        // once the send\recv channels have been established, we perform a user defined number of evolutionary cycles
        // of the solution population ...

        for(meta.run.eval.id = 1; meta.run.eval.id <= meta.topologies.max_fit_evals; meta.run.eval.id++) {

            LOG(5, 0, 0, "begin eval %d on child topology %d\r\n", meta.run.eval.id, it->id);

            solver.evolve(solver.solutions, meta, *it, solutions_evolve);
            
        }

        meta.topologies.population.push_back(*it);

    }
    
    // children evaluated for fitness, perform survival selection ...
    
    if(meta.variant.isle.id == 0) {
    
        // similar to crossover, we only need to perform survival selection on the root island,
        // so only truncate the the worst individuals in the root island population ...
        
        std::sort(meta.topologies.population.begin(), meta.topologies.population.end(), compare_fitness<topology>);
        std::reverse(meta.topologies.population.begin(), meta.topologies.population.end());

        meta.topologies.population.erase(meta.topologies.population.begin()+(meta.topologies.mu-1), meta.topologies.population.end());

        meta.topologies.fitness();
        
    }
   
    
}

#pragma mark FUNCTION: topology_eval()

// apply a provided topology to the communication pattern to be
// used by objective 1 for island migrations ...
//
// ea metadata modified to reflect current topology, evaluate
// iterate as determined by ea parameter @config.txt[topo_evals]
// perform n evolutionary cycles of objective (1)
//

//void topology_evaluate(ea &meta, ea &solver, topology &t) {
//
//    LOG(4, solver.variant.isle.id, 0, "applying topology %d to %d islands ...\r\n", t.id, solver.variant.islands);
//
//    t.fitness = 0.0;
//    t.round_fitness = 0.0;
//
//    // distribute the send\recv channels from the passed topology to all islands ...
//
//    t.apply(solver.variant.isle, t);
//
//    LOG(8, 0, 0, "evaluating topology %d\r\n", t.id);
//
//    // evaluate the supplied topology using the solution (o1) evolution cycle
//    // which will assign a fitness value based on the time spent during migration
//
//    // this function should be called iteratively for each individual in the
//    // topology (o2) population, so it is assumed that the following loop will
//    // execute µ * o2.max_fitness_evals (for example 256 * 100 = 25,600) times
//
//    // the (o1) population may assert a maxmimum iteration constraint as well, so
//    //
//
//    for(meta.run.eval.id = 1; meta.run.eval.id%(meta.topologies.max_fit_evals+1) != 0; meta.run.eval.id++) {
//
//        LOG(5, 0, 0, "begin meta ea run %d solution objective on topology %d\r\n", solver.run.eval.id, t.id);
//
//        solver.evolve(solver.solutions, meta, t, solutions_evolve);
//
//        fflush(config::topo_stats_out);
//
//    }
//
//    LOG(8, solver.variant.isle.id, 0, "evaluated topology id=%d rounds=%d fitness=%3.10f\r\n", t.id, t.rounds, t.fitness);
//
//}

#pragma mark FUNCTION: ea_init()

// collect ea properties based on the runtime environment and config parameters.
// successful completion results in a data heirarchy describing the ea components
// needed to implement a parallelized (island model) multi-objective evolutionary algorithm.

//ea ea_init() {
//
//    // each process ("island") calling this function will reference its own unique instance
//    // with some values determined from the process ("island") id ...
//
//    ea multi;
//
//    multi.variant.start = std::clock();
//
//    // initialize MPI environment ...
//
//    MPI_Init(NULL, NULL);
//    MPI_Comm_size(MPI_COMM_WORLD, &multi.variant.islands);
//    MPI_Comm_rank(MPI_COMM_WORLD, &multi.variant.isle.id);
//
//    // load configuration items ...
//
//    config::load("config.txt", multi.variant.islands, multi.variant.isle.id);
//
//    for(int i=0; i<multi.variant.islands; i++) { multi.variant.island_ids.push_back(i); }
//
//    // ----- start mpi derived datatypes...
//
//    MPI_Datatype ptype;
//
//    MPI_Type_contiguous(128, MPI_CHAR, &ptype);
//    MPI_Type_commit(&ptype);
//
//    MPI_Datatype sol_types[11] = { MPI_CHAR, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, ptype };
//
//    MPI_Aint sol_offsets[11] = { 0,                                  // id                       + 1 lluint
//        sizeof(char)*64,                                             // input                    + 10 double
//        sizeof(double)*DIM + sizeof(char)*64,                        // fitness                  + 11 double
//        sizeof(double)*(DIM+1) + sizeof(char)*64,                    // selection_distribution   + 12 double
//        sizeof(double)*(DIM+2) + sizeof(char)*64,                    // group                    + 13 double
//        sizeof(double)*(DIM+3) + sizeof(char)*64,                    // source                   + 1 int
//        sizeof(double)*(DIM+3) + sizeof(char)*64 + (sizeof(int)),    // locale                   + 1 int
//        sizeof(double)*(DIM+3) + sizeof(char)*64 + (sizeof(int)*2),  // migrations               + 1 int
//        sizeof(double)*(DIM+3) + sizeof(char)*64 + (sizeof(int)*3),  // selected                 + 1 int
//        sizeof(double)*(DIM+3) + sizeof(char)*64 + (sizeof(int)*4),  // survival                 + 1 int
//        sizeof(double)*(DIM+3) + sizeof(char)*64 + (sizeof(int)*5)   // parents                  + 2 int
//    };
//
//    int sol_lengths[11] = { 64, DIM, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
//
//    MPI_Type_create_struct(11, sol_lengths, sol_offsets, sol_types, &multi.variant.solution_type);
//    MPI_Type_commit(&multi.variant.solution_type);
//
//    MPI_Datatype visa_types[4] = { MPI_INT, MPI_INT, MPI_INT, MPI_CHAR };
//    MPI_Aint visa_offsets[4] = { 0, sizeof(int), sizeof(int)*2, sizeof(int)*3 };
//
//    int visa_lengths[4] = { 1, 1, 1, 64 };
//
//    MPI_Type_create_struct(4, visa_lengths, visa_offsets, visa_types, &multi.variant.visa_type);
//    MPI_Type_commit(&multi.variant.visa_type);
//
//    // ----- end derived mpi datatypes
//
//    multi.eval = eval_init(0);
//    multi.eval.stats.init();
//
//    multi.variant.tcomm = MPI_COMM_WORLD;
//    multi.variant.island_size = config::island_mu;
//    multi.variant.isle.init();
//
//    multi.solutions.mu = config::ea_1_mu;
//    multi.solutions.lambda = config::ea_1_lambda;
//    multi.solutions.max_runs = config::ea_1_runs;
//    multi.solutions.mutation_rate = config::ea_1_mutation_rate;
//    multi.solutions.max_evo_evals = config::ea_1_max_evo_evals;
//    multi.solutions.max_fit_evals = config::ea_1_max_fit_evals;
//    multi.solutions.eval.id = 1;
//
//    multi.topologies.mu = config::ea_2_mu;
//    multi.topologies.lambda = config::ea_2_lambda;
//    multi.topologies.max_runs = config::ea_2_runs;
//    multi.topologies.mutation_rate = config::ea_2_mutation_rate;
//    multi.topologies.max_evo_evals = config::ea_2_max_evo_evals;
//    multi.topologies.max_fit_evals = config::ea_2_max_fit_evals;
//    multi.topologies.islands = multi.variant.islands;
//    multi.topologies.eval.id = 1;
//
//    multi.run.stats.init_duration = ( std::clock() - multi.variant.start ) / (double) CLOCKS_PER_SEC;
//
//    if(multi.variant.isle.id == 0) {
//        multi.offsets = generate_offsets(-2.5, 2.5, .5);
//    }
//
//    MPI_Bcast(&multi.offsets, DIM, MPI_DOUBLE, 0, multi.variant.tcomm);
//
//    // collect the time consumed by all islands in this initialization ...
//
//    // TODO: the init time gather segfaults on higher core count runs, so may need to debug at some point, but currently the init duration is somewhat insigificant
//    // MPI_Gather(&local_init_duration, 1, MPI_DOUBLE, &multi.run.stats.init_duration, 1, MPI_DOUBLE, 0, multi.meta.isle.tcomm);
//
//    LOG(2, multi.variant.isle.id, 0, "world size: %d\r\n", multi.variant.islands);
//    LOG(2, multi.variant.isle.id, 0, "mu mode: %d\r\n", config::mu_mode);
//    LOG(2, multi.variant.isle.id, 0, "global mu %s %d\r\n", config::mu_msg, config::ea_1_mu);
//    LOG(2, multi.variant.isle.id, 0, "island mu %s %d\r\n", config::subpop_msg, config::island_mu);
//    LOG(2, multi.variant.isle.id, 0, "island lambda %s %d\r\n", config::lambda_msg, stoi(config::items["island_lambda"]));
//
//    return multi;
//
//}

#pragma mark FUNCTION: topology_populate()

// returns a collection of randomly generated adjaceny matrices, representing                           |
// an island (communication) topology.  accepts a reference to a list of island                         |
// (process) identifiers from @ea{@meta{@island_ids[param:0}} to use as indices.                        |

void topologies_populate(ea &meta) {
    
    meta.topologies.population.clear();
    
    LOG(4, meta.variant.isle.id, 0, "rank %d (root) initializing topology population ... \r\n", meta.variant.isle.id);

    for(int i=0; i<meta.topologies.mu; i++) {
        
        topology t;
        
        t.id = i;
        t.rounds = 0;
        t.world_size = meta.variant.islands;
        t.fitness = 0.0;
        t.channel_count = 0;
        t.round_fitness = 0.0;
        t.selection_distribution = 0.0;
        t.channels = {};
        
        if(meta.variant.isle.id != 0) {
            LOG(8, 0, 0, "island %d (leaf) adding topology stub %d, returning ... \r\n", meta.variant.isle.id, i);
            meta.topologies.population.push_back(t);
        } else {
            LOG(8, 0, 0, "island %d (root) topology %d\r\n", meta.variant.isle.id, i);
            topology::create::dynamic(t);
            LOG(6, 0, 0, "dynamic topology %d created\r\n", i);
            meta.topologies.population.push_back(t);
            meta.run.stats.total_channels += t.channel_count;
        }
        
    }
    
    // we have our initial topology population, *NOT* evaluated for fitness
    
    LOG(4, meta.variant.isle.id, 0, "initialized objective (topology) population size %lu, \r\n", meta.topologies.population.size());
    
}

void benchmark_topology(ea &meta) {
    
    meta.topologies.population.clear();
    
    topology t;
    
    t.rounds = 0;
    t.world_size = meta.variant.islands;
    t.fitness = 0.0;
    t.round_fitness = 0.0;
    t.selection_distribution = 0.0;
    t.channels = {};
    t.channels.resize(meta.variant.islands);
    
    for(int i=0; i<meta.variant.islands; i++) {
        
        t.channels[i].senders = {};
        t.channels[i].receivers = {};
        
        int next = i+1 < meta.variant.islands ? i+1 : 0;
        int prev = i-1 < 0 ? (int)meta.variant.islands-1 : i-1;

        t.channels[i].senders.push_back(prev);
        t.channels[i].receivers.push_back(next);
        
    }
    
    t.channel_count = t.world_size * 2;

    meta.topologies.population.push_back(t);

}


#endif /* ea_h */

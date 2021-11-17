//
//  solver.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 10/17/21.
//  Copyright © 2021 Bradley Morgan. All rights reserved.
//

#ifndef solver_h
#define solver_h


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

    while (isle.metrics.value.cpd[i] < r ) { LOG(10, 0, 0, "island %d, cpd[%d] = %f, r = %f\r\n", isle.id, i, isle.metrics.value.cpd[i], r); i++; }

    LOG(8, 0, 0, "island %d selected individual %d cpd[%d] = %f = ...\r\n", isle.id, i, i, isle.metrics.value.cpd[i]);
    
    isle.population[i].selected++;
    
    LOG(5, 0, 0, "sol %s selected=%d", isle.population[i].id, isle.population[i].selected);
    
    solution s = isle.population[i];
    
    return s;

}


void select_survivors(island &isle, std::vector<solution> &children, int island_mu) {

    isle.average_fitness();
    
    LOG(6, 0, 0, "EA SELECTION: island %d pre-survival population size: %lu, average fitness: %f, best fit: %f\r\n", isle.id, isle.population.size(), isle.metrics.value.average_fitness, isle.population[0].fitness);
    
    // truncation: add new children to the population, and then kill the weakest
    
    isle.population.insert(isle.population.end(), children.begin(), children.end());
    
    LOG(4, isle.id, 0, "EA SELECTION: island %d added %lu child solutions, population size = %lu, begin selecting %d survivors\r\n", isle.id, children.size(), isle.population.size(), island_mu);
    
    std::sort(isle.population.begin(), isle.population.end(), compare_fitness<solution>);
    std::reverse(isle.population.begin(), isle.population.end());
    
    LOG(4, isle.id, 0, "EA SELECTION: island %d fitness individual[0] = %f, individual[mu] = %f\r\n", isle.id, isle.population[0].fitness, isle.population[isle.population.size()-1].fitness);
    
    isle.population.erase(isle.population.begin()+island_mu, isle.population.end());
    
    for_each(isle.population.begin(), isle.population.end(), [](solution &g) { g.survival++; });
    
    isle.average_fitness();
    
    LOG(4, isle.id, 0, "EA SELECTION: island %d survival population size: %lu, average fitness: %f, best fit: %f\r\n", isle.id, isle.population.size(), isle.metrics.value.average_fitness, isle.population[0].fitness);
    
    LOG(5, 0, 0, "EA SELECTION: island %d survival population size: %lu, average fitness: %f, best fit: %f\r\n", isle.id, isle.population.size(), isle.metrics.value.average_fitness, isle.population[0].fitness);
    
}


void mutate(solution &mutant) {

    for(int i=0; i<DIM; i++){
        
        if(rand()/(RAND_MAX+1.0) < config::ea_1_mutation_rate) {
     
            mutant.input[i] = drand(-5.12, 5.12);
            
        }
        
    }
    
}

std::vector<solution> crossover(ea_solver &solver) {
      
    solver.variant.isle.cpd();
    
    std::vector<solution> children;
    
    for(int i = 0; i < config::island_lambda; i++) {
        
        LOG(5, 0, 0, "island %d creating objective<solution> child %d ...\r\n",
            solver.variant.isle.id, i);
        
        solution p1 = select_parent(solver.variant.isle);
        LOG(6, 0, 0, "island %d selected p1<solution> cpd = %f, fitness = %f ...\r\n", solver.variant.isle.id, p1.selection_distribution, p1.fitness);
        
        solution p2 = select_parent(solver.variant.isle);
        LOG(6, 0, 0, "island %d selected p2<solution> cpd = %f, fitness = %f...\r\n", solver.variant.isle.id, p2.selection_distribution, p2.fitness);
        
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
        
        if(rand()/(RAND_MAX+1.0) < solver.solutions.mutation_rate) {
            
            LOG(6, 0, 0, "island %d mutating child<solution> %d ...\r\n", solver.variant.isle.id, i);
            
            mutate(child);
            
        }
        
        solver.solutions.begin(solver.solutions.run.eval, solver);
        
        child.fitness = offset_rastrigin(child.input, solver.offsets);
        
        solver.solutions.end(solver.solutions.run.eval, solver);
        
        if(child.fitness > (p1.fitness / 4)) {
            LOG(3, 0, 0, "LOW island %d %f<->%f child<solution> %d fitness %f\r\n", solver.variant.isle.id, p1.fitness, p2.fitness, i, child.fitness);
        } else {
            LOG(8, 0, 0, "child<solution> %d fitness %f\r\n", i, child.fitness);
        }
        
        child.source = solver.variant.isle.id;
        child.locale = solver.variant.isle.id;
        child.migrations = 0;
        
        strcpy(child.parents[0], p1.id);
        strcpy(child.parents[1], p2.id);
        
        LOG(5, 0, 0, "child id %s p1 %s p2 %s\r\n", child.id, child.parents[0], child.parents[1]);
        
        children.push_back(child);
        
    }
    
    LOG(6, 0, 0, "rank %d returning %lu child solutions from crossover\r\n", solver.variant.isle.id, children.size());
    
    return children;
    
}

#pragma mark FUNCTION: solution_populate()

// returns a vector of a randomly generated offset rastrigin solution population                        |
// accepts an array of floating point numbers of n-dimensions @offsets[param:0]                         |
// used to calculate the solution fitness.                                                              |

void solution_populate(ea_solver &solver) {
    
    LOG(6, 0, 0, "rank %d of %d entered solution_populate\r\n", solver.variant.isle.id, solver.variant.islands);

    if(solver.variant.isle.id != 0) {
        LOG(6, 0, 0, "rank %d leaving solution_populate\r\n", solver.variant.isle.id);
        return;
    }
    
    // the "solution" datatype represents a single offset rastrigin solution as
    // an array of size DIM = @dim[config.txt:12] holding the solution's randomly
    // generated gene values.
    
    LOG(4, solver.variant.isle.id, 0, "SOLVER: initializing objective %d population, %d individuals ...\r\n", solver.solutions.id, solver.solutions.mu);
    LOG(6, 0, 0, "island %d initializing mu=%d solutions ...\r\n", solver.variant.isle.id, solver.solutions.mu);

    for(int i=0; i < solver.solutions.mu; i++) {
        
        LOG(6, 0, 0, "creating solution %d ...\r\n", i);
        
        solution p;
        
        p.source = solver.variant.isle.id;
        p.migrations++;
        
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
        
        solver.solutions.begin(solver.solutions.run.eval, solver);
        
        p.fitness = offset_rastrigin(p.input, solver.offsets);
        
        solver.solutions.end(solver.solutions.run.eval, solver);
        
        solver.solutions.run.eval.end();
        
        LOG(6, 0, 0, "rank %d assigned solution fitness = %f\r\n", solver.variant.isle.id, p.fitness);
        
        LOG(6, 0, 0, "rank %d adding solution %d with fitness %f to population, current size = %lu\r\n", solver.variant.isle.id, i, p.fitness, solver.solutions.population.size());
        
        solver.solutions.aggregate.value.fitness += p.fitness;

        //visa v(solver.solutions.run.eval.id, p.source, p.source, p.id);

        //multi.variant.isle.visas.push_back(v);

        solver.solutions.population.push_back(p);
    
        LOG(6, 0, 0, "island %d (root) initialized solution %d with fitness %f ...\r\n", solver.variant.isle.id, i, p.fitness);
        
    }
    
    LOG(4, solver.variant.isle.id, 0, "SOLVER initialized objective %d population with size %lu, total fitness = %f\r\n", solver.solutions.id, solver.solutions.population.size(), solver.solutions.aggregate.value.fitness);
    LOG(6, 0, 0, "SOLVER ISLAND %d initialized objective %d local population with size %lu ...\r\n", solver.variant.isle.id, solver.solutions.id, solver.solutions.population.size());

    if(solver.solutions.run.eval.stats.global_best_fitness == 0.0) {
        solver.solutions.run.eval.stats.global_best_fitness = solver.solutions.population.data()[0].fitness;
    }
    
    LOG(6, 0, 0, "rank %d leaving solution_populate\r\n", solver.variant.isle.id);
    
    // we have our initial primary population with the calculated fitnesses
    
}

#pragma mark FUNCTION: solution_scatter()

// separate the single full population from the root process to subpopulations across all processes ...

void solution_scatter(ea_solver &solver) {
    
    LOG(6, 0, 0, "rank %d entered solution_scatter\r\n", solver.variant.isle.id);

//    LOG(6, 0, 0, "rank %d of %d resized local island population to %lu\r\n", solver.variant.isle.id, solver.variant.islands, solver.variant.isle.population.size());
    
//    if(solver.variant.isle.id == 0) {
//        LOG(4, solver.variant.isle.id, 0, "rank 0 scattering population root size = %lu mem 0 = %f...\r\n", solver.solutions.population.size(), solver.solutions.population[0].fitness);
//    } else {
//          LOG(4, 0, 0, "rank %d instantiating scatter with = %lu subpopulation size ...\r\n", solver.variant.isle.id, solver.variant.isle.population.size());
//    }

    double scatter_start = MPI_Wtime();
    
    LOG(6, 0, 0, "\r\nISLAND %d of %d SCATTER INIT: distributing %d of %d solutions\r\n", solver.variant.isle.id, solver.variant.islands, config::island_mu, solver.solutions.mu);
    
    MPI_Scatter(&solver.solutions.population[0], config::island_mu, solver.variant.solution_type, &solver.variant.isle.population[0], config::island_mu, solver.variant.solution_type, 0, solver.variant.tcomm);
    
    double scatter_end = MPI_Wtime();
    double scatter_time = scatter_end - scatter_start;

    solver.variant.isle.average_fitness();

    LOG(4, 0, 0, "\r\nISLAND %d of %d SCATTER END: received %lu solutions average fitness %f\r\n", solver.variant.isle.id, solver.variant.islands, solver.variant.isle.population.size(), solver.variant.isle.metrics.value.average_fitness);

    MPI_Reduce(&scatter_time, &solver.solutions.run.eval.stats.total_scatter_time, 1, MPI_DOUBLE, MPI_SUM, 0, solver.variant.tcomm);

    solver.solutions.run.eval.stats.total_scatter_time += scatter_time;

    LOG(6, 0, 0, "RETURN: rank %d leaving solution_scatter\r\n", solver.variant.isle.id);
    
}

#pragma mark FUNCTION: solution_eval()

void solutions_evolve(ea_solver &solver, ea_meta &meta, topology &t) {
    
    LOG(6, 0, 0, "ISLAND %d ENTER objective<solution> evolution cycle %d, topology %d, eval %d\r\n", solver.variant.isle.id, solver.solutions.run.eval.id, t.id, solver.solutions.run.eval.id);
    LOG(8, solver.variant.isle.id, 0, "calculating island %d cpd, topology %d, eval %d\r\n", solver.variant.isle.id, t.id, solver.solutions.run.eval.id);
    
    solver.solutions.fitness();
    solver.solutions.cpd();
    
    LOG(6, solver.variant.isle.id, 0, "performing crossover island %d, topology %d, eval %d\r\n", solver.variant.isle.id, t.id, solver.solutions.run.eval.id);

    // crossover function performs parent selection iteratively, creating n child solutions with fitness calculated ...
    
    std::vector<solution> children = crossover(solver);

    // select the best solutions and remove n children ...
    
    select_survivors(solver.variant.isle, children, solver.variant.island_size);

    LOG(6, solver.variant.isle.id, 0, "island migrations ...\r\n");
    
    // perform migration of best solutions using the currently applied topology ...
    
    if(solver.solutions.cycle.id%config::migration_interval == 0) {
        
        double migrate_start = MPI_Wtime();
        
        // issue migration imports and exports ...
                
        island::migration::send(solver.variant.isle, solver.variant.solution_type, solver.solutions.run.eval.id);
        island::migration::receive(solver.variant.isle, solver.variant.solution_type, solver.solutions.run.eval.id);
        
        double migrate_end = MPI_Wtime();
        double migrate_time = migrate_end - migrate_start;
        
        LOG(8, 0, 0, "migrate start = %3.10f, migrate end = %3.10f, migrate time = %3.10f\r\n", migrate_start, migrate_end, migrate_time);
        
        solver.solutions.run.eval.stats.total_migrate_time += migrate_time;
     
        // aggregate the migration time from each island into the topology fitness ...
        
        MPI_Reduce(&migrate_time, &t.round_fitness, 1, MPI_DOUBLE, MPI_SUM, 0, solver.variant.tcomm);
        
    }

    LOG(4, solver.variant.isle.id, 0, "ISLAND %d topology %d, eval %d, round %d, tfitness = %f\r\n", solver.variant.isle.id, t.id, solver.solutions.run.eval.id, t.rounds, t.fitness);
    
    LOG(7, 0, 0, "rank %d topology = %d rounds = %d eval = %d\r\n", solver.variant.isle.id, t.id, t.rounds, meta.topologies.run.eval.id);

    if(solver.variant.isle.id == 0) {
        
        t.rounds++;
        t.total_migration_time += t.round_fitness;
        t.round_fitness = t.round_fitness * -1;
        //t.fitness = (t.total_migration_time / t.rounds) * -1;
        //t.fitness += (MPI_Wtime() - meta.topologies.run.eval.stats.eval_start) * -1;
        meta.topologies.run.eval.stats.topo_migrate_time += t.total_migration_time;
        
    }
    
    LOG(8, solver.variant.isle.id, 0, "rounds=%d, total_time=%013.10f, round_fit=%013.10f, fit=%013.10f\r\n", t.rounds, t.total_migration_time, t.round_fitness, t.fitness);
    
    // output for various intervals ...
    
    if(solver.solutions.cycle.id%config::ea_1_log_interval == 0) {
    
        LOG(6, 0, 0, "\r\nISLAND %d of %d GATHER INIT: %d solutions from %d islands = (%d * %d) = mu = %d\r\n", solver.variant.isle.id, solver.variant.islands, config::island_mu, solver.variant.islands, config::island_mu, solver.variant.islands, solver.solutions.mu);
        
        double gather_start = MPI_Wtime();
        
        // gather island subpopulations back into the aggregate population on rank 0 ...
        
        MPI_Gather(&solver.variant.isle.population[0], solver.variant.island_size, solver.variant.solution_type, &solver.solutions.population[0], solver.variant.island_size, solver.variant.solution_type, 0, solver.variant.tcomm);
        
        double gather_end = MPI_Wtime();
        double gather_time = gather_end - gather_start;
        
        solver.solutions.run.eval.stats.total_gather_time += gather_time;
        
        LOG(4, 0, 0, "\r\nISLAND %d of %d GATHER END: %lu solutions from %d islands = (%d / %d) = island_mu = %d\r\n", solver.variant.isle.id, solver.variant.islands, solver.solutions.population.size(), solver.variant.islands, solver.solutions.mu, solver.variant.islands, config::island_mu);
        
    }
    
    if(solver.solutions.cycle.id%config::ea_1_log_interval == 0 && solver.variant.isle.id == 0) {
        
        LOG(5, 0, 0, "population size %lu, member = %2.10f\r\n", solver.variant.isle.population.size(), solver.variant.isle.population[0].fitness);
    
        if(solver.variant.isle.id == 0) {
            
            //log_fn_eval_stats(solver, meta, t);
            solver.solutions.log_stats(solver.solutions.cycle, solver, meta, t);
            
        }
        
    }
    
    if(solver.solutions.cycle.id%config::ea_2_log_interval == 0) {
        
        LOG(5, 0, 0, "population size %lu, member = %2.10f\r\n", solver.variant.isle.population.size(), solver.variant.isle.population[0].fitness);
        
        if(solver.variant.isle.id == 0) {
            //log_fn_topology_stats(solver, meta, t);
            meta.topologies.log_stats(meta.topologies.run.eval, solver, meta, t);
        }
        
    }
    
    if(solver.solutions.cycle.id%config::ea_1_population_log_interval == 0) {
        
        log_pop_stats(solver, solver.variant.isle, solver.variant.visa_type);
        
    }

}

#pragma mark EA::FUNCTION::SOLVER solver_begin()

// TODO: check eval critera ... how many evals do we need?
// TODO: number of runs is for statistical certainty ---

void solver_begin(ea_meta &meta, ea_solver &solver, topology &t, int runs = config::ea_1_runs, int cycles = config::ea_1_max_evo_cycles) {

    LOG(6, 0, 0, "BEGIN ISLAND %d objective<solutions> EVOLUTION (objective<topology> %d) AT SOLVER[%d,%d] META[%d,%d]\r\n", solver.variant.isle.id, t.id, solver.solutions.run.id, solver.solutions.run.eval.id, meta.topologies.run.id, meta.topologies.run.eval.id);
    
    // 𝑆𝑟𝑚𝑎𝑥 * 𝑆𝑒𝑚𝑎𝑥 nested iterations ...
    
    solver.start = MPI_Wtime();
    
    for(solver.solutions.run.id = 1; solver.solutions.run.id <= runs; solver.solutions.run.id++) {
     
        solver.solutions.begin(solver.solutions.run, solver);
        
        meta.variant = solver.variant;
        
        solver.solutions.populate(solver, solution_populate);
        solver.solutions.distribute(solver, solution_scatter);
        
        t.apply(solver.variant.isle, t);
        
        for(solver.solutions.cycle.id = 1; solver.solutions.cycle.id <= cycles; solver.solutions.cycle.id++) {

            solver.solutions.begin(solver.solutions.cycle, solver);
            
            solver.solutions.evolve(solver, meta, t, solutions_evolve);
            
            solver.solutions.end(solver.solutions.cycle, solver);
                
            t.fitness -= solver.solutions.cycle.duration;
            
        }

        solver.solutions.end(solver.solutions.run, solver);
        
    }
    
    solver.duration = MPI_Wtime() - solver.start;
    
    LOG(6, 0, 0, "END ISLAND %d objective<solutions> EVOLUTION (objective<topology> %d) AT SOLVER[%d,%d] META[%d,%d]\r\n", solver.variant.isle.id, t.id, solver.solutions.run.id, solver.solutions.cycle.id, meta.topologies.run.id, meta.topologies.run.eval.id);
    
}


#endif /* solver_h */

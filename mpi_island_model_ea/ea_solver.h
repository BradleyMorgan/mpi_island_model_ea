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

    while (isle.metrics.value.cpd[i] < r ) { LOG(10, 0, 0, "island %d, cpd[%d] = %f, r = %f\r\n", mpi.id, i, isle.metrics.value.cpd[i], r); i++; }

    LOG(8, 0, 0, "island %d selected individual %d cpd[%d] = %f = ...\r\n", mpi.id, i, i, isle.metrics.value.cpd[i]);
    
    isle.population[i].selected++;
    
    LOG(5, 0, 0, "sol %s selected=%d", isle.population[i].id, isle.population[i].selected);
    
    solution s = isle.population[i];
    
    return s;

}


void select_survivors(island &isle, std::vector<solution> &children, int island_mu) {

    isle.average_fitness();
    
    LOG(6, 0, 0, "EA SELECTION: island %d pre-survival population size: %lu, average fitness: %f, best fit: %f\r\n", mpi.id, isle.population.size(), isle.metrics.value.average_fitness, isle.population[0].fitness);
    
    // truncation: add new children to the population, and then kill the weakest
    
    isle.population.insert(isle.population.end(), children.begin(), children.end());
    
    LOG(4, mpi.id, 0, "EA SELECTION: island %d added %lu child solutions, population size = %lu, begin selecting %d survivors\r\n", mpi.id, children.size(), isle.population.size(), island_mu);
    
    std::sort(isle.population.begin(), isle.population.end(), compare_fitness<solution>);
    std::reverse(isle.population.begin(), isle.population.end());
    
    LOG(4, mpi.id, 0, "EA SELECTION: island %d fitness individual[0] = %f, individual[mu] = %f\r\n", mpi.id, isle.population[0].fitness, isle.population[isle.population.size()-1].fitness);
    
    isle.population.erase(isle.population.begin()+island_mu, isle.population.end());
    
    for_each(isle.population.begin(), isle.population.end(), [](solution &g) { g.survival++; });
    
    isle.average_fitness();
    
    LOG(4, mpi.id, 0, "EA SELECTION: island %d survival population size: %lu, average fitness: %f, best fit: %f\r\n", mpi.id, isle.population.size(), isle.metrics.value.average_fitness, isle.population[0].fitness);
    
    LOG(5, 0, 0, "EA SELECTION: island %d survival population size: %lu, average fitness: %f, best fit: %f\r\n", mpi.id, isle.population.size(), isle.metrics.value.average_fitness, isle.population[0].fitness);
    
}


void mutate(solution &mutant) {

    for(int i=0; i<DIM; i++){
        
        if(rand()/(RAND_MAX+1.0) < config::ea_1_mutation_rate) {
     
            mutant.input[i] = drand(-5.12, 5.12);
            
        }
        
    }
    
}

std::vector<solution> crossover(ea_solver &solver) {
      
    solver.model.isle.cpd();
    
    std::vector<solution> children;
    
    for(int i = 0; i < solver.model.island_lambda; i++) {
        
        LOG(5, 0, 0, "island %d creating objective<solution> child %d ...\r\n",
            mpi.id, i);
        
        solution p1 = select_parent(solver.model.isle);
        LOG(6, 0, 0, "island %d selected p1<solution> cpd = %f, fitness = %f ...\r\n", mpi.id, p1.selection_distribution, p1.fitness);
        
        solution p2 = select_parent(solver.model.isle);
        LOG(6, 0, 0, "island %d selected p2<solution> cpd = %f, fitness = %f...\r\n", mpi.id, p2.selection_distribution, p2.fitness);
        
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
            
            LOG(6, 0, 0, "island %d mutating child<solution> %d ...\r\n", mpi.id, i);
            
            mutate(child);
            
        }
        
        
        // child evaluation
        
        solver.ea::begin(solver.solutions, solver.solutions.run.cycle.eval, &child);
        
        child.fitness = offset_rastrigin(child.input, solver.offsets);
        
        solver.ea::end(solver.solutions, solver.solutions.run.cycle.eval, &child);
        
        // end child evaluation
        
        
        if(child.fitness > (p1.fitness / 4)) {
            LOG(4, 0, 0, "LOW island %d %f<->%f child<solution> %d fitness %f\r\n", mpi.id, p1.fitness, p2.fitness, i, child.fitness);
        } else {
            LOG(8, 0, 0, "child<solution> %d fitness %f\r\n", i, child.fitness);
        }
        
        child.source = mpi.id;
        child.locale = mpi.id;
        
        child.migrations = 0;
        
        strcpy(child.parents[0], p1.id);
        strcpy(child.parents[1], p2.id);
        
        LOG(5, 0, 0, "child id %s p1 %s p2 %s\r\n", child.id, child.parents[0], child.parents[1]);
        
        children.push_back(child);
        
        solver.solutions.run.cycle.eval.local = &child;
        
    }
    
    LOG(6, 0, 0, "rank %d returning %lu child solutions from crossover\r\n", mpi.id, children.size());
    
    return children;
    
}

#pragma mark FUNCTION: solution_populate()

// returns a vector of a randomly generated offset rastrigin solution population                        |
// accepts an array of floating point numbers of n-dimensions @offsets[param:0]                         |
// used to calculate the solution fitness.                                                              |

void solution_populate(ea_solver &solver) {
    
    LOG(6, 0, 0, "\r\nSOLVER RUN %d POPULATE BEGIN: rank %d of %d\r\n", solver.solutions.run.id, mpi.id, mpi.size);

    if(mpi.id != 0) {
        LOG(6, 0, 0, "\r\nSOLVER RUN %d POPULATE RETURN: rank %d leaving solution_populate\r\n", solver.solutions.run.id, mpi.id);
        return;
    }
    
    // the "solution" datatype represents a single offset rastrigin solution as
    // an array of size DIM = @dim[config.txt:12] holding the solution's randomly
    // generated gene values.
    
    LOG(4, mpi.id, 0, "\r\nSOLVER RUN %d POPULATE STEP 1: island %d generate mu=%d genomes ...\r\n", mpi.id, mpi.id, solver.solutions.mu);

    for(int i=0; i < solver.solutions.mu; i++) {
        
        LOG(6, 0, 0, "\r\nSOLVER RUN %d POPULATE STEP 1a: creating solution %d ...\r\n", solver.solutions.run.id, i);
        
        solution p;
        
        //p.source = mpi.id;
        //p.migrations = 0;
        
        strcpy(p.parents[0], "0");
        strcpy(p.parents[1], "0");
        
        LOG(6, 0, 0, "\r\nSOLVER RUN %d POPULATE STEP 1b: rank %d solution %s source=%d migrations=%d parents=%s,%s", solver.solutions.run.id, mpi.id, p.id, p.source, p.migrations, p.parents[0], p.parents[1]);
        
        LOG(10, mpi.id, 0, "\r\nSOLVER RUN %d POPULATE STEP 1c: island %d solution %s genome: [", solver.solutions.run.id, mpi.id, p.id);
        
        for (int j = 0; j < DIM; j++) {
            p.input[j] = drand(-5.12, 5.12); // rastrigin says: x[i] ∈ [-5.12,5.12]
            p.group += p.input[j];
            LOG(10, mpi.id, 0, "%f", p.input[j]);
        }
        
        LOG(10, mpi.id, 0, "]\r\n");
        
        //solver.ea::begin(solver.solutions, solver.solutions.run.cycle.eval, &p);
        
        p.fitness = offset_rastrigin(p.input, solver.offsets);
        
        //solver.ea::end(solver.solutions, solver.solutions.run.cycle.eval);
        
        LOG(6, 0, 0, "\r\nSOLVER RUN %d POPULATE STEP 1d: island %d adding solution %s (%d of %d) fitness=%f to population (size = %lu)\r\n", solver.solutions.run.id, mpi.id, p.id, i, solver.solutions.mu, p.fitness, solver.solutions.population.size());
        
        solver.solutions.aggregate.value.fitness += p.fitness;

        //visa v(solver.solutions.run.cycle.eval.id, p.source, p.source, p.id);

        //multi.variant.isle.visas.push_back(v);

        solver.solutions.population.push_back(p);
        
    }
    
    LOG(4, mpi.id, 0, "\r\nSOLVER RUN %d POPULATE:  initialized objective %d population with size %lu, total fitness = %f\r\n", solver.solutions.run.id, solver.solutions.id, solver.solutions.population.size(), solver.solutions.aggregate.value.fitness);
    
    LOG(6, 0, 0, "rank %d leaving solution_populate\r\n", mpi.id);
    
    // solver.solutions.log_population(solver.solutions.run);
    
    // we have our initial primary population with the calculated fitnesses
    
}

#pragma mark FUNCTION: solution_scatter()

// called by all islands
// separate the single full population from the root process to subpopulations across all processes ...

void solution_scatter(ea_solver &solver) {
    
    LOG(3, 0, 0, "rank %d entered solution_scatter\r\n", mpi.id);

    double scatter_start = MPI_Wtime();
    
    LOG(4, mpi.id, 0, "\r\nSOLVER RUN %d SCATTER BEGIN: island %d distributing %d solutions to each of %d islands \r\n", solver.solutions.run.id, mpi.id, solver.model.island_mu, mpi.size);
    
    MPI_Scatter(&solver.solutions.population[0], solver.model.island_mu, solver.model.solution_type, &solver.model.isle.population[0], solver.model.island_mu, solver.model.solution_type, 0, solver.model.tcomm);
    
    solver.solutions.run.stats.local_scatter_t.value = MPI_Wtime() - scatter_start;
    
    solver.model.isle.average_fitness();
    
    // assuming that distribution marks the beginning of a new run,
    // assign set the run local solution reference to the first individual
    
    solver.solutions.run.local = &solver.solutions.population[0];

    LOG(3, 0, 0, "\r\nSOLVER RUN %d SCATTER END: island %d received %lu solutions average fitness %f\r\n", mpi.id, mpi.size, solver.model.isle.population.size(), solver.model.isle.metrics.value.average_fitness);

    // solution distribution occurs at the beginning of solver run
    // log time metrics to the current run
        
    LOG(3, mpi.id, 0, "\r\n\r\nSOLVER RUN %d SCATTER RETURN: rank %d leaving solution_scatter cycle: %f %f %f %f", solver.solutions.run.id, mpi.id, solver.solutions.run.stats.min_scatter_t.value, solver.solutions.run.stats.max_scatter_t.value, solver.solutions.run.stats.sum_scatter_t.value, solver.solutions.run.stats.avg_scatter_t);
    
}

#pragma mark FUNCTION: solution_eval()

void solutions_evolve(ea_solver &solver, ea_meta &meta, topology &t) {
    
    LOG(7, 0, 0, "(%d) ENTER objective<solution> evolution cycle %d, topology %d, eval %d\r\n", mpi.id, solver.solutions.run.cycle.eval.id, t.id, solver.solutions.run.cycle.eval.id);
    LOG(8, mpi.id, 0, "calculating island %d cpd, topology %d, eval %d\r\n", mpi.id, t.id, solver.solutions.run.cycle.eval.id);
    
    // calculate population fitness and corresponding sorted probability distribution
    
    //solver.solutions.fitness();
    //solver.solutions.cpd();
    
    LOG(6, mpi.id, 0, "(%d) topology %d, eval %d\r\n", mpi.id, t.id, solver.solutions.run.cycle.eval.id);

    // crossover function performs parent selection iteratively, creating n child solutions with fitness calculated ...
    
    std::vector<solution> children = crossover(solver);

    // select the best solutions and remove n children ...
    
    select_survivors(solver.model.isle, children, solver.model.island_mu);

    LOG(6, mpi.id, 0, "island migrations ...\r\n");
    
    // perform migration of best solutions using the currently applied topology ...
    
    if(solver.solutions.run.cycle.id%config::migration_interval == 0) {
        
        double migrate_start = MPI_Wtime();
        
        // issue migration imports and exports ...

        solver.model.isle.stats.departures = 0;
        solver.model.isle.stats.arrivals = 0;
        
        int scount = 0;
        int rcount = 0;
        
        LOG(7, mpi.id, 0, "%d %d SEND [%d,%d,%d] -> chans=%lu ->\r\n", solver.model.isle.id, solver.solutions.run.cycle.id, t.stats.departures, t.stats.arrivals, t.stats.migrations, solver.model.isle.receivers.size());
        island::migration::send(solver.model.isle, solver.model.solution_type, solver.solutions.run.cycle.eval.id);
        LOG(7, mpi.id, 0, "%d %d SEND [%d,%d,%d]\r\n", solver.model.isle.id, solver.solutions.run.cycle.id, t.stats.departures, t.stats.arrivals, t.stats.migrations);
        
        LOG(7, mpi.id, 0, "%d %d RECV [%d,%d,%d] -> chans=%lu ->\r\n", solver.model.isle.id, solver.solutions.run.cycle.id, t.stats.departures, t.stats.arrivals, t.stats.migrations, solver.model.isle.receivers.size());
        island::migration::receive(solver.model.isle, solver.model.solution_type, solver.solutions.run.cycle.eval.id);
        LOG(7, mpi.id, 0, "%d %d RECV [%d,%d,%d]\r\n", solver.model.isle.id, solver.solutions.run.cycle.id, t.stats.departures, t.stats.arrivals, t.stats.migrations);
        
        LOG(9, 0, 0, "%d %d MIGR REDUCE SEND [%d,%d,%d]\r\n", solver.model.isle.id, solver.solutions.run.cycle.id, solver.model.isle.stats.departures, solver.model.isle.stats.arrivals, solver.model.isle.stats.migrations);
        
        MPI_Reduce(&solver.model.isle.stats.departures, &scount, 1, MPI_INT, MPI_SUM, 0, solver.model.isle.tcomm);
        MPI_Reduce(&solver.model.isle.stats.arrivals, &rcount, 1, MPI_INT, MPI_SUM, 0, solver.model.isle.tcomm);
        
        solver.model.isle.stats.migrations += scount + rcount;
        
        LOG(9, mpi.id, 0, "%d %d MIGR REDUCE RECV [%d,%d,%d]\r\n", solver.model.isle.id, solver.solutions.run.cycle.id, scount, rcount, solver.model.isle.stats.migrations);
        
        if(mpi.id == 0) {
            t.stats.departures += scount;
            t.stats.arrivals += rcount;
            t.stats.migrations += scount + rcount;
        }
        
        solver.solutions.run.cycle.stats.local_migration_t.value = MPI_Wtime() - migrate_start;
        solver.solutions.run.stats.local_migration_t.value += solver.solutions.run.cycle.stats.local_migration_t.value;
        
        LOG(8, 0, 0, "migrate start = %3.10f, migrate end = %3.10f, migrate time = %3.10f\r\n", migrate_start, MPI_Wtime() - migrate_start, solver.solutions.run.cycle.stat_send.value);
        
    }

    LOG(4, 0, 0, "ISLAND %d MIGRATION topology=%d migr[%d,%d,%d] schan=%lu rchan=%lu, fit=%f\r\n", mpi.id, t.id, t.stats.departures, t.stats.arrivals, t.stats.migrations, solver.model.isle.senders.size(), solver.model.isle.receivers.size(), t.fitness);

}

#pragma mark EA::FUNCTION::SOLVER solver_begin()

// TODO: check eval critera ... how many evals do we need?
// TODO: number of runs is for statistical certainty ---

void solver_begin(ea_meta &meta, ea_solver &solver, topology &t, int runs = config::ea_1_o1_max_runs, int cycles = config::ea_1_o1_max_cycles) {

    LOG(6, 0, 0, "BEGIN ISLAND %d objective<solutions> EVOLUTION (objective<topology> %d) AT SOLVER[%d,%d] META[%d,%d]\r\n", mpi.id, t.id, solver.solutions.run.id, solver.solutions.run.cycle.eval.id, meta.topologies.run.id, meta.topologies.run.cycle.eval.id);
    
    // 𝑆𝑟𝑚𝑎𝑥 * 𝑆𝑒𝑚𝑎𝑥 nested iterations ...
    
    t.apply(solver.model.isle, t);
    
    for(solver.solutions.run.id = 1; solver.solutions.run.id <= runs; ++solver.solutions.run.id) {
    
        if(mpi.id == 0) {
            solver.offsets = generate_offsets(-2.5, 2.5, .5);
        }
        
        MPI_Bcast(&solver.offsets, DIM, MPI_DOUBLE, 0, solver.model.tcomm);
        
        solver.ea::begin(solver.solutions);
        
        solver.solutions.populate(solver, solution_populate);
        solver.solutions.distribute(solver, solution_scatter);
        
        for(solver.solutions.run.cycle.id = 1; solver.solutions.run.cycle.id <= cycles; solver.solutions.run.cycle.id++) {
        
            solver.ea::begin(solver.solutions, solver.solutions.run.cycle, solver.solutions.run.local);
            
            solver.solutions.evolve(solver, meta, t, solutions_evolve);
            
            solver.ea::end(solver.solutions, solver.solutions.run.cycle, solver.solutions.run.cycle.local);
            
            LOG(5, mpi.id, 0, "**** TOPOLOGY %d FITNESS AT CYCLE %d: %f ****\r\n", t.id, solver.solutions.run.cycle.id, t.fitness);
            
        }
        
        solver.ea::end(solver.solutions);
        
        LOG(4, mpi.id, 0, "**** TOPOLOGY %d FITNESS AT RUN %d: %f ****\r\n", t.id, solver.solutions.run.id, t.fitness);
        
    }
    
    LOG(6, mpi.id, 0, "**** TOPOLOGY %d FITNESS AT EA END %d: %f ****\r\n", t.id, solver.solutions.run.id, t.fitness);
    
    LOG(5, 0, 0, "END ISLAND %d objective<solutions> EVOLUTION (objective<topology> %d) AT SOLVER[%d,%d] META[%d,%d]\r\n", mpi.id, t.id, solver.solutions.run.id, solver.solutions.run.cycle.id, meta.topologies.run.id, meta.topologies.run.cycle.eval.id);
    
}

#endif /* solver_h */

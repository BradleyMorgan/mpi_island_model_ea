//
//  ea.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 3/30/20.
//  Copyright © 2020 Bradley Morgan. All rights reserved.
//

#ifndef ea_h
#define ea_h

#include "utility.h"
#include "solution.h"
#include "island.h"
#include "topology.h"
#include "stats.h"

#pragma mark DATATYPE: @ea_eval{}

// track eval stats, etc. per island

struct ea_eval {
  
    int id;
    double start;
    estats stats;
};

#pragma mark DATATYPE: @ea_run{}

// track run stats, etc. per island

struct ea_run {
    
    int id;
    double start;
    rstats stats;
    ea_eval eval;
    
};

template<typename O> struct objective {
    
    int id;
    int runs;
    int evals;
    int run_id;
    int eval_id;
    int world_size;
    int mu;
    
    double total_fitness;
    
    std::vector<O> population;
    
    std::vector<double> cpd;
    
    struct calculate {
     
        static void cpd(objective &o);
        static void total_fitness(objective &o);
        
    };
    
};

template<typename C> void objective<C>::calculate::total_fitness(objective &o) {
    
    o.total_fitness = 0.0;
    
    for(int i=0; i<o.population.size(); i++) {
    
        o.total_fitness += o.population[i].fitness;
        
    }
    
}


#pragma mark FUNCTION: @objective<C>::calculate::cpd()

// calculate a population's cumulative probability distribution
// (and individual selection distribution) for fitness proportional selection 

template<typename C> void objective<C>::calculate::cpd(objective &o) {
    
    LOG(10, 0, 0, "generating cpd\r\n");
    
    double cumulative_probability = 0.0;
    
    o.cpd.clear();
        
    island::calculate::total_fitness(o);
    
    LOG(6, 0, 0, "island %d cpd calculation total fitness = %f", o.total_fitness)
    std::reverse(o.population.begin(), o.population.end());
    
    for(int i=0; i<o.population.size(); i++) {

        o.population[i].selection_distribution = (double)o.population[i].fitness / o.total_fitness;

        cumulative_probability += o.population[i].selection_distribution;
        o.cpd.push_back(cumulative_probability);
        
    }
    
}

#pragma mark DATATYPE: @ea_meta{}

// describes the properties needed for parallel (island model) EA

struct ea_meta {

    std::clock_t start;

    int islands;
    int island_size;
    int root = 0;

    // isle holds context specific population and migration data
    // from the perspective of *this* island (mpi rank)
    // see @island{} datatype
    
    island isle;
    
    // TODO: maintaining a list of island_ids doesn't seem necessary
    
    std::vector<int> island_ids;
    

    // track derived mpi datatypes used to transfer custom structs
    
    MPI_Datatype solution_type;
    MPI_Datatype topology_type;
    MPI_Datatype comm_type;
    
    MPI_Comm tcomm;
    
};


#pragma mark DATATYPE: @ea{}

// wrapper type to serve as the trunk for ea properties

struct ea {
    
    ea_meta meta;
    ea_run run;
    ea_eval eval;
    
    objective<topology> topologies;
    objective<solution> solutions;
    
    std::array<double, DIM> offsets;

    template<typename func, typename otype>
    void populate(func function, objective<otype> &o) {
        function(*this, o);
    };
    
    template<typename func, typename otype>
    void distribute(func function, objective<otype> &o) {
        function(*this, o);
    };
    
    template<typename func, typename otype>
    void evaluate(func function, objective<otype> &o, const int &index) {
        function(*this, o, index);
    };
    
    template<typename func, typename otype>
    void parent_selection(func function, objective<otype> &o) {
        function(this->meta.isle, o);
    };
    
};

solution pselect(island &isle, objective<solution> &o) {
    
    solution p;

    // implementation uses the single armed roulette wheel approach to select
    // an individual from the population
    
    int i = 1;
    
    // random double uniformly distributed between 0 and 1
    
    double r = ((double)rand()/(double)RAND_MAX);
    
    // spin the wheel
    
    while (isle.cpd[i] < r ) { i++; }
    
    p = isle.population[i];
    
    return p;
        
}

// pick a parent from the population using the provided method ...

solution select_parent(island &isle) {
    
    solution p;

    // implementation uses the single armed roulette wheel approach to select
    // an individual from the population
    
    int i = 1;
    
    // random double uniformly distributed between 0 and 1
    
    LOG(8, 0, 0, "random double uniformly distributed between 0 and 1\r\n");
    
    double r = ((double)rand()/(double)RAND_MAX);
    
    // spin the wheel
    
    LOG(10, 0, 0, "spin the wheel ...\r\n");
    
    while (isle.cpd[i] < r ) { LOG(8, 0, 0, "island %d, cpd[%d] = %f, r = %f\r\n", isle.id, i, isle.cpd[i], r); i++; }
    
    LOG(8, 0, 0, "island %d selected individual %d ...\r\n", isle.id, i);
    
    p = isle.population[i];
    
    return p;
        
}


void select_survivors(island &isle, std::vector<solution> &children, int island_mu) {

    LOG(6, 0, 0, "island %d survival, population before: %lu, fitness before: %f, best fit: %f\r\n", isle.id, isle.population.size(), isle.average_fitness, isle.population[0].fitness);
    
    // truncation: add new children to the population, and then kill the weakest
    
    isle.population.insert(isle.population.end(), children.begin(), children.end());
    
    std::sort(isle.population.begin(), isle.population.end(), compare_fitness);
    std::reverse(isle.population.begin(), isle.population.end());
    
    isle.population.erase(isle.population.begin()+island_mu, isle.population.end());
    
    island::calculate::average_fitness(isle);
    
    LOG(6, 0, 0, "island %d survival, population after: %lu, fitness after: %f, best fit: %f\r\n", isle.id, isle.population.size(), isle.average_fitness, isle.population[0].fitness);
    
}


void mutate(solution &mutant) {

    mutant.input[0] = drand(-5.12, 5.12);
    
}

std::vector<solution> crossover(island &isle, std::array<double, DIM> &offsets) {

    std::vector<solution> children;
    
    for(int i = 0; i < config::lambda; i++) {
        
        LOG(5, 0, 0, "island %d creating child %d ...\r\n", isle.id, i);
        
        solution p1 = select_parent(isle);
        LOG(6, 0, 0, "island %d selected p1 cpd = %f, fitness = %f ...\r\n", isle.id, p1.selection_distribution, p1.fitness);
        solution p2 = select_parent(isle);
        LOG(6, 0, 0, "island %d selected p2 cpd = %f, fitness = %f...\r\n", isle.id, p2.selection_distribution, p2.fitness);
        
        solution child;
        
        for(int j=0; j<DIM; j++) {
            if(rand()%2 == 1) {
                int test = rand()%DIM;
                child.input[j] = p1.input[test];
                LOG(4, 0, 0, "assigning gene %d = %f from parent 1 to child %d (%d) ...\r\n", j, child.input[j], i, test);
            } else {
                int test = rand()%DIM;
                child.input[j] = p2.input[test];
                LOG(4, 0, 0, "assigning gene %d = %f from parent 2 to child %d (%d) ...\r\n", j, child.input[j], i, test);
            }
        }
        
        if(rand()/(RAND_MAX+1.0) < config::mutation_rate) {
            
            LOG(6, 0, 0, "island %d mutating child %d ...\r\n", isle.id, i);
            
            mutate(child);
            
        }
        
        child.fitness = offset_rastrigin(child.input, offsets);
        
        if(child.fitness > (p1.fitness / 4)) {
           LOG(4, 0, 0, "LOW island %d %f<->%f child %d fitness %f\r\n", isle.id, p1.fitness, p2.fitness, i, child.fitness);
        } else {
            LOG(4, 0, 0, "child %d fitness %f\r\n", i, child.fitness);
        }
        
        children.push_back(child);
        
    }
    
    return children;
    
}

#pragma mark FUNCTION: topology_populate()

// returns a collection of randomly generated adjaceny matrices, representing                           |
// an island (communication) topology.  accepts a reference to a list of island                         |
// (process) identifiers from @ea{@meta{@island_ids[param:0}} to use as indices.                        |

void topologies_populate(ea &multi, objective<topology> &o) {
    
    LOG(4, multi.meta.isle.id, 0, "rank %d (root) initializing topology population ... \r\n", multi.meta.isle.id);
    
    for(int i=0; i<config::topo_mu; i++) {
        
        topology t;
        
        t.id = i;
        t.rounds = 0;
        t.world_size = multi.meta.islands;
        t.fitness = 0.0;
        t.channel_count = 0;
        t.round_fitness = 0.0;
        t.selection_distribution = 0.0;
        
//        if(multi.meta.isle.id != 0) {
//            LOG(8, 0, 0, "island %d (leaf) adding topology stub %d, returning ... \r\n", multi.meta.isle.id, i);
//            multi.topologies.population.push_back(t);
//            return;
//        }
        
        LOG(8, 0, 0, "island %d (root) topology %d\r\n", multi.meta.isle.id, i);
        
        topology::create::dynamic(t);
        
        LOG(6, 0, 0, "dynamic topology %d created\r\n", i);
        
        multi.topologies.population.push_back(t);
        
    }
    
    LOG(4, multi.meta.isle.id, 0, "initialized topology population size %lu, \r\n", multi.topologies.population.size());
    
    // we have our initial topology population, evaluated for fitness
    
}


#pragma mark FUNCTION: solution_populate()

// returns a vector of a randomly generated offset rastrigin solution population                        |
// accepts an array of floating point numbers of n-dimensions @offsets[param:0]                         |
// used to calculate the solution fitness.                                                              |

void solution_populate(ea &multi, objective<solution> &o) {
    
    if(multi.meta.isle.id != 0) { return; }
    
    // the "solution" datatype represents a single offset rastrigin solution as
    // an array of size DIM = @dim[config.txt:12] holding the solution's randomly
    // generated gene values.
    
    LOG(6, 0, 0, "initializing mu=%d solutions ...\r\n", config::mu);
    
    for(int i=0; i<config::mu; i++) {
        
        solution p;
        
        for (int j = 0; j < DIM; j++) {
            p.input[j] = drand(-5.12, 5.12); // rastrigin says: x[i] ∈ [-5.12,5.12]
        }
        
        p.fitness = offset_rastrigin(p.input, multi.offsets);
        
        o.population.push_back(p);
        
        o.total_fitness += p.fitness;
        
        LOG(6, 0, 0, "island %d (root) initialized solution %d with fitness %f ...\r\n", multi.meta.isle.id, i, p.fitness);
        
    }
    
    LOG(6, 0, 0, "%d solutions initialized ...\r\n", config::mu);
    
    std::sort(o.population.begin(), o.population.end(), compare_fitness);
    std::reverse(o.population.begin(), o.population.end());
    
    multi.eval.stats.global_best_fitness = o.population[0].fitness;
    
    // we have our initial primary population with the calculated fitnesses
    
}

#pragma mark FUNCTION: solution_scatter()

// separate the single full population from the root process to subpopulations across all processes ...

void solution_scatter(ea &multi, objective<solution> &o) {
    
    multi.meta.isle.population.clear();
    multi.meta.isle.population.resize(multi.meta.island_size);
    
    LOG(4, multi.meta.isle.id, 0, "scattering population root size = %lu mem 0 = %f...\r\n", multi.meta.isle.population.size(), o.population[0].fitness);
    
    double scatter_start = MPI_Wtime();
    
    // MPI_Scatter(&source, count, type, &target, count, type, source_island, comm)
    
    MPI_Scatter(&o.population[0], multi.meta.island_size, multi.meta.solution_type, &multi.meta.isle.population[0], multi.meta.island_size, multi.meta.solution_type, 0, multi.meta.tcomm);
    
    double scatter_end = MPI_Wtime();
    double scatter_time = scatter_end - scatter_start;
    
    LOG(4, 0, 0, "population scattered, island %d population %lu mem 0 = %f...\r\n", multi.meta.isle.id, multi.meta.isle.population.size(), multi.meta.isle.population[0].fitness);
    
    island::calculate::average_fitness(multi.meta.isle);
    
    LOG(4, 0, 0, "island %d population size %lu, average fitness: %f\r\n", multi.meta.isle.id, multi.meta.isle.population.size(), multi.meta.isle.average_fitness);
    
    multi.eval.stats.total_scatter_time += scatter_time;
    
}

#pragma mark FUNCTION: solution_eval()

void solutions_evolve(topology &t, ea &multi) {
   
    LOG(5, multi.meta.isle.id, 0, "beginning objective<solution> evaluation %d, topology %d\r\n", multi.solutions.eval_id, t.id);
    
    multi.eval.start = std::clock();

    LOG(8, multi.meta.isle.id, 0, "calculating island %d cpd, topology %d, eval %d\r\n", multi.meta.isle.id, t.id, multi.solutions.eval_id);
    
    island::calculate::cpd(multi.meta.isle);
    
    LOG(6, multi.meta.isle.id, 0, "performing crossover island %d, topology %d, eval %d\r\n", multi.meta.isle.id, t.id, multi.solutions.eval_id);

    std::vector<solution> children = crossover(multi.meta.isle, multi.offsets);

    LOG(5, multi.meta.isle.id, 0, "island %d created %lu child solutions, population size = %lu ... selecting %d survivors  ...\r\n", multi.meta.isle.id, children.size(), multi.meta.isle.population.size(), multi.meta.island_size);

    select_survivors(multi.meta.isle, children, multi.meta.island_size);

    LOG(6, multi.meta.isle.id, 0, "island migrations ...\r\n");
    
    if(multi.solutions.eval_id%1 == 0) {
        
        double migrate_start = MPI_Wtime();
        
        island::migration::send(multi.meta.isle, multi.meta.isle.tcomm);
        island::migration::receive(multi.meta.isle, multi.meta.isle.tcomm);
        
        double migrate_end = MPI_Wtime();
        double migrate_time = migrate_end - migrate_start;
        
        LOG(8, 0, 0, "migrate start = %3.10f, migrate end = %3.10f, migrate time = %3.10f\r\n", migrate_start, migrate_end, migrate_time);
        
        multi.eval.stats.total_migrate_time += migrate_time;
     
        MPI_Reduce(&migrate_time, &t.fitness, 1, MPI_DOUBLE, MPI_SUM, 0, multi.meta.tcomm);
        
    }

    LOG(6, multi.meta.isle.id, 0, "rank %d (root) topology %d, eval %d, round %d, fitness = %f\r\n", multi.meta.isle.id, t.id, multi.solutions.eval_id, t.rounds, t.fitness);
    
    LOG(6, 0, 0, "rank %d topology = %d rounds = %d eval = %d\r\n", multi.meta.isle.id, t.id, t.rounds, multi.topologies.eval_id);
    
    LOG(6, 0, 0, "gathering population, subpopulation %d size %lu avg fitness %f...\r\n", multi.meta.isle.id, multi.meta.isle.population.size(), multi.meta.isle.average_fitness);
    
    double gather_start = MPI_Wtime();
    
    // gather island subpopulations back into the aggregate population on rank 0 ...
    
    multi.solutions.population.clear();
    multi.solutions.population.resize(multi.meta.island_size);
    
    MPI_Gather(&multi.meta.isle.population[0], multi.meta.island_size, multi.meta.solution_type, &multi.solutions.population[0], multi.meta.island_size, multi.meta.solution_type, 0, multi.meta.tcomm);
    
    double gather_end = MPI_Wtime();
    double gather_time = gather_end - gather_start;
                
    LOG(6, multi.meta.isle.id, 0, "population gathered, size %lu ...\r\n", multi.solutions.population.size());
                
    multi.eval.stats.total_gather_time += gather_time;
    
    //if(t.fitness >= 0.0) {
    t.fitness = t.fitness * -1;
    t.round_fitness += t.fitness;
    t.rounds++;
    //} else {
    //    LOG(1, 0, 0, "FOUND NEGATIVE FITNESS: %2.10f in topology %d\r\n", t.fitness, t.id);
   // }
    
    if(multi.solutions.eval_id%10 == 0) {
        LOG(6, 0, 0, "island %d average fit %f\r\n", multi.meta.isle.id, multi.meta.isle.average_fitness);
    }
    
    if(multi.solutions.eval_id%100 == 0) {
                
        LOG(4, 0, 0, "population size %lu, member = %2.10f\r\n", multi.meta.isle.population.size(), multi.meta.isle.population[0].fitness);
        
        std::sort(multi.meta.isle.population.begin(), multi.meta.isle.population.end(), compare_fitness);
        std::reverse(multi.meta.isle.population.begin(), multi.meta.isle.population.end());
        
        for(int i=0; i<multi.meta.isle.population.size(); i++) {
            if(multi.meta.isle.population[i].fitness == -0.0) {
                LOG(4, 0, 0, "ZERO fitness found on island %d for individual %d\r\n", multi.meta.isle.id, i);
            }
        }
        //std::sort(multi.topologies.population.begin(), multi.topologies.population.end(), compare_topo_fitness);
        //std::reverse(multi.topologies.population.begin(), multi.topologies.population.end());
    
    }
    
    if(multi.meta.isle.id == 0 && multi.solutions.eval_id%multi.topologies.evals == 0) {
        
        for(int i=0; i<multi.solutions.population.size(); i++) {
            if(multi.solutions.population[i].fitness == -0.0) {
                LOG(4, 0, 0, "ZERO fitness found for individual %d\r\n", i);
            }
        }
        
        log_fn_eval_stats(multi.solutions.population, multi.topologies.population, multi.run.id, multi.solutions.eval_id, multi.eval.stats, multi.run.stats, t);
        
    }
        
}

#pragma mark FUNCTION: eval_init()

// start of new eval, initialize ...

ea_eval eval_init(int id) {
    
    ea_eval eval;
    
    eval.id = id;
    eval.stats.eval_start = std::clock();
    
    return eval;
    
}

#pragma mark FUNCTION: init_topology_run()

// start of new run, initialize ...

ea_run run_init(int id, ea &multi) {
    
    LOG(10, 0, 0, "island %d run %d initialization\r\n", multi.meta.isle.id, id);

    multi.run.id = id;
    multi.run.start = MPI_Wtime();
    
    if(multi.meta.isle.id != 0) {
        return multi.run;
    }
  
    multi.run.eval = eval_init(0);
    
    LOG(4, 0, 0, "initializing rastrigin population...\r\n");
    multi.solutions.population.clear();
    multi.solutions.population.resize(config::mu);
    
    LOG(4, 0, 0, "initializing topology population...\r\n");
    multi.topologies.population.clear();
    multi.topologies.population.resize(config::topo_mu);

    LOG(4, 0, 0, "initializing run variables...\r\n");
    multi.offsets = generate_offsets(-2.5, 2.5, .5);
        
    return multi.run;
    
}

void solution_evaluate(ea &multi, objective<solution> &o, const int &index) {
    
    LOG(8, 0, 0, "evaluating topology %d\r\n", multi.topologies.population[index].id);

    //for(o.eval_id = 1; o.eval_id <= o.evals; o.eval_id++) {
        
        solutions_evolve(multi.topologies.population[index], multi);
    
    //}
    
}

#pragma mark FUNCTION: topology_eval()

// apply a provided topology to the communication pattern to be
// used by objective 1 for island migrations ...
//
// ea metadata modified to reflect current topology, evaluate
// iterate as determined by ea parameter @config.txt[topo_evals]
// perform n evolutionary cycles of objective (1)
//

void topology_evaluate(ea &multi, objective<topology> &o, int index) {
    
    LOG(4, multi.meta.isle.id, 0, "applying topology %d to %d islands ... %d\r\n", o.population[index].id, multi.meta.islands, multi.meta.islands);
    
    o.population[index].fitness = 0.0;
    o.population[index].round_fitness = 0.0;
    
    o.population[index].apply(multi.meta.isle, o.population[index]);
    
    LOG(8, 0, 0, "evaluating topology %d\r\n", o.population[index].id);
    
    for(o.eval_id = 1; o.eval_id <= o.evals; o.eval_id++) {

        LOG(5, 0, 0, "begin eval %d solution objective on topology %d\r\n", o.eval_id, o.population[index].id);

        multi.solutions.eval_id++;
        
        multi.evaluate(solution_evaluate, multi.solutions, index);
        
    }

}

#pragma mark FUNCTION: ea_init()

// collect ea properties based on the runtime environment and config parameters.
// successful completion results in a data heirarchy describing the ea components
// needed to implement a parallelized (island model) multi-objective evolutionary algorithm.

ea ea_init() {
  
    // each process ("island") calling this function will reference its own unique instance
    // with some values determined from the process ("island") id ...
    
    ea multi;
    
    multi.meta.start = std::clock();
    
    // initialize MPI environment ...

    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &multi.meta.islands);
    MPI_Comm_rank(MPI_COMM_WORLD, &multi.meta.isle.id);

    // load configuration items ...

    config::load("config.txt", multi.meta.islands, multi.meta.isle.id);
    
    for(int i=0; i<multi.meta.islands; i++) { multi.meta.island_ids.push_back(i); }

    // MPI derived datatype for solution individual ...

    int sol_lengths[4] = { DIM, 1, 1, 1 };
    MPI_Aint sol_displacements[4] = { 0, sizeof(double)*DIM, sizeof(double)*(DIM+1) };
    MPI_Datatype sol_types[4] = { MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE };
    MPI_Type_create_struct(4, sol_lengths, sol_displacements, sol_types, &multi.meta.solution_type);
    MPI_Type_commit(&multi.meta.solution_type);
    
    multi.meta.island_size = config::mu / multi.meta.islands;
    
    
    //  attempted a derived mpi datatype for topology, but never got it to work
    
    //    int comm_lengths[4] = { 1, 1, 1, 1 };
    //    MPI_Aint comm_displacements[4] = { 0, sizeof(int), sizeof(int)+sizeof(double),
    //    (int)((sizeof(int)*2)+(sizeof(double)))*multi.meta.islands };
    //    MPI_Datatype comm_types[4] = { MPI_INT, MPI_DOUBLE, MPI_INT, MPI_INT };
    //    MPI_Type_create_struct(4, comm_lengths, comm_displacements, comm_types, &multi.meta.comm_type);
    //    MPI_Type_commit(&multi.meta.comm_type);
    //
    //    int topo_lengths[6] = { 1, 1, 1, 1, 1, multi.meta.islands};
    //    MPI_Aint topo_displacements[6] = { 0, sizeof(int), sizeof(int)*2, sizeof(double)+(sizeof(int)*2),
    //    (int)((sizeof(double)*2)+(sizeof(int)*2))*multi.meta.islands, (sizeof(double)*3)+(sizeof(int)*2) };
    //    MPI_Datatype topo_types[6] = { MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, multi.meta.comm_type};
    //    MPI_Type_create_struct(6, topo_lengths, topo_displacements, topo_types, &multi.meta.topology_type);
    //    MPI_Type_commit(&multi.meta.topology_type);
    
    //  communication duplication testing for implicit rearrangment, placeholder ...
    //  duplicate the default MPI_COMM_WORLD communicator so that it can be manipulated ...
    
    //  LOG(4, 0, 0, "duplicating MPI_COMM_WORLD ...\r\n");
    //  MPI_Comm_dup(MPI_COMM_WORLD, &multi.meta.tcomm);

    multi.meta.tcomm = MPI_COMM_WORLD;
    
    multi.meta.isle.init();
    multi.meta.isle.population.resize(multi.meta.islands);
    
    multi.solutions.mu = config::mu;
    multi.solutions.runs = config::runs;
    multi.solutions.evals = config::evals;
    multi.solutions.eval_id = 0;
    
    multi.topologies.mu = config::topo_mu;
    //multi.topologies.runs = config::topo_runs;
    multi.topologies.evals = config::topo_evals;
    multi.topologies.world_size = multi.meta.islands;
    multi.topologies.eval_id = 0;
    
    multi.run.stats.init_duration = ( std::clock() - multi.meta.start ) / (double) CLOCKS_PER_SEC;
    
    return multi;
    
}


#endif /* ea_h */

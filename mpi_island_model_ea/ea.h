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


struct comp_fitness {
  template <typename otype>
  bool operator() (const otype &a, const otype &b) const {
    return a.fitness < b.fitness;
  }
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
    
    LOG(6, 0, 0, "calculating total fitness ...\r\n");
    
    objective::calculate::total_fitness(o);
    
    LOG(6, 0, 0, "sorting population descending fitness ...\r\n");
    
    std::sort(o.population.begin(), o.population.end(), comp_fitness());
    std::reverse(o.population.begin(), o.population.end());
    
    LOG(6, 0, 0, "objective %d cpd calculation total fitness = %f", o.id, o.total_fitness)

    for(int i=0; i<o.population.size(); i++) {

        o.population[i].selection_distribution = (double)o.population[i].fitness / o.total_fitness;

        LOG(8, 0, 0, "calculated solution %d fitness %f selection distribution = %f\r\n", i, o.population[i].fitness, o.population[i].selection_distribution);
        
        cumulative_probability += o.population[i].selection_distribution;
        
        LOG(8, 0, 0, "solution %d cumulative prob = %f\r\n", i, cumulative_probability);
        
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
    
    MPI_Datatype solution_type;
    MPI_Datatype topology_type;
    MPI_Datatype comm_type;
    
    MPI_Datatype MPI_TOPOLOGY;
    MPI_Datatype MPI_CHANNEL;
    
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
    void evaluate(func function, objective<otype> &o, topology &t) {
        function(*this, o, t);
    };
    
    template<typename func, typename otype>
    otype selection(func function, objective<otype> &o, std::vector<otype> &population) {
        return function(*this, o, population);
    };
    
    template<typename func, typename otype>
    otype crossover(func function, objective<otype> &o, std::vector<otype> &population) {
        return function(*this, o, population);
    };
    
    void run_init() {
        
        this->solutions.population.clear();
        this->meta.isle.population.clear();
        
        this->eval.stats.init();
        this->solutions.eval_id = 0;
        this->topologies.eval_id = 0;
        
        this->run.start = MPI_Wtime();
        
    }
    
    void run_end() {
        
        if(this->meta.isle.id != 0) { return; }
        
        double run_end = MPI_Wtime();
        
        std::fprintf(config::run_stats_out, "%d,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%2.10f,%d,%d,%2.10f,%2.10f,%2.10f\r\n", run.id, this->eval.stats.global_best_fitness, this->eval.stats.average_local_best_fitness, this->eval.stats.average_global_best_fitness, this->eval.stats.total_scatter_time, this->eval.stats.total_gather_time, this->eval.stats.total_migrate_time, run_end - this->run.start, this->run.stats.init_duration, this->meta.islands, this->meta.island_size, this->eval.stats.global_best_topo_fitness, this->eval.stats.average_local_best_topo_fitness, this->eval.stats.average_global_best_topo_fitness);
        
        fflush(config::run_stats_out);
        
    }
    
};

template<typename otype>
otype parent(ea &multi, objective<otype> &o, std::vector<otype> &population) {
    
    otype p;

    // implementation uses the single armed roulette wheel approach to select
    // an individual from the population
    
    int i = 1;
    
    // random double uniformly distributed between 0 and 1
    
    double r = ((double)rand()/(double)RAND_MAX);
    
    // spin the wheel
    
    objective<otype>::calculate::cpd(o);
    
    while (o.cpd[i] < r ) { i++; }
    
    p = o.population[i];
    
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
        
        LOG(5, 0, 0, "island %d creating objective<solution> child %d ...\r\n", isle.id, i);
        
        solution p1 = select_parent(isle);
        LOG(6, 0, 0, "island %d selected p1<solution> cpd = %f, fitness = %f ...\r\n", isle.id, p1.selection_distribution, p1.fitness);
        solution p2 = select_parent(isle);
        LOG(6, 0, 0, "island %d selected p2<solution> cpd = %f, fitness = %f...\r\n", isle.id, p2.selection_distribution, p2.fitness);
        
        solution child;
        
        for(int j=0; j<DIM; j++) {
            if(rand()%2 == 1) {
                int test = rand()%DIM;
                child.input[j] = p1.input[test];
                LOG(6, 0, 0, "assigning gene %d = %f from p1<solution> to child %d (%d) ...\r\n", j, child.input[j], i, test);
            } else {
                int test = rand()%DIM;
                child.input[j] = p2.input[test];
                LOG(6, 0, 0, "assigning gene %d = %f from p2<solution> to child %d (%d) ...\r\n", j, child.input[j], i, test);
            }
        }
        
        if(rand()/(RAND_MAX+1.0) < config::mutation_rate) {
            
            LOG(6, 0, 0, "island %d mutating child<solution> %d ...\r\n", isle.id, i);
            
            mutate(child);
            
        }
        
        child.fitness = offset_rastrigin(child.input, offsets);
        
        if(child.fitness > (p1.fitness / 4)) {
           LOG(4, 0, 0, "LOW island %d %f<->%f child<solution> %d fitness %f\r\n", isle.id, p1.fitness, p2.fitness, i, child.fitness);
        } else {
            LOG(4, 0, 0, "child<solution> %d fitness %f\r\n", i, child.fitness);
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
    
    o.population.clear();
    
    LOG(4, multi.meta.isle.id, 0, "rank %d (root) initializing topology population ... \r\n", multi.meta.isle.id);
    
    //if(multi.meta.isle.id != 0) { return; }
    
    for(int i=0; i<config::topo_mu; i++) {
        
        topology t;
        
        t.id = i;
        t.rounds = 0;
        t.world_size = multi.meta.islands;
        t.fitness = 0.0;
        t.channel_count = 0;
        t.round_fitness = 0.0;
        t.selection_distribution = 0.0;
        t.channels = {};
        
        if(multi.meta.isle.id != 0) {
            LOG(8, 0, 0, "island %d (leaf) adding topology stub %d, returning ... \r\n", multi.meta.isle.id, i);
            multi.topologies.population.push_back(t);
        } else {
            LOG(8, 0, 0, "island %d (root) topology %d\r\n", multi.meta.isle.id, i);
            topology::create::dynamic(t);
            LOG(6, 0, 0, "dynamic topology %d created\r\n", i);
            multi.topologies.population.push_back(t);
        }
        
    }
    
    LOG(4, multi.meta.isle.id, 0, "initialized objective (topology) population size %lu, \r\n", multi.topologies.population.size());
    
    // we have our initial topology population, evaluated for fitness
    
}

void benchmark_topology(ea &multi) {
    
    multi.topologies.population.clear();
    
    topology t;

    t.id = 0;
    t.rounds = 0;
    t.world_size = multi.meta.islands;
    t.fitness = 0.0;
    t.round_fitness = 0.0;
    t.selection_distribution = 0.0;
    t.channels = {};
    t.channels.resize(multi.meta.islands);
    
    for(int i=0; i<multi.meta.islands; i++) {
        
        t.channels[i].senders = {};
        t.channels[i].receivers = {};
        
        int next = i+1 < multi.meta.islands ? i+1 : 0;
        int prev = i-1 < 0 ? (int)multi.meta.islands-1 : i-1;

        t.channels[i].senders.push_back(prev);
        t.channels[i].receivers.push_back(next);
     
        //multi.meta.isle.receivers[0] = next;
        //multi.meta.isle.senders[0] = prev;
        
    }
    
    t.channel_count = t.world_size * 2;
    
    multi.topologies.population.push_back(t);

}

#pragma mark FUNCTION: solution_populate()

// returns a vector of a randomly generated offset rastrigin solution population                        |
// accepts an array of floating point numbers of n-dimensions @offsets[param:0]                         |
// used to calculate the solution fitness.                                                              |

void solution_populate(ea &multi, objective<solution> &o) {
    
    LOG(4, multi.meta.isle.id, 0, "initializing objective (solution) population ...\r\n");
    LOG(6, 0, 0, "rank %d entered solution_populate\r\n", multi.meta.isle.id);
    
    if(multi.meta.isle.id != 0) { return; }
    
    // the "solution" datatype represents a single offset rastrigin solution as
    // an array of size DIM = @dim[config.txt:12] holding the solution's randomly
    // generated gene values.
    
    LOG(6, 0, 0, "island %d (root) initializing mu=%d solutions ...\r\n", multi.meta.isle.id, config::mu);
    
    o.population.resize(config::mu);
    
    for(int i=0; i<config::mu; i++) {
        
        LOG(6, 0, 0, "creating solution %d ...\r\n", i);
        
        solution p;
        
        LOG(6, 0, 0, "rank %d assigning solution values\r\n", multi.meta.isle.id);
        
        for (int j = 0; j < DIM; j++) {
            p.input[j] = drand(-5.12, 5.12); // rastrigin says: x[i] ∈ [-5.12,5.12]
        }
        
        LOG(6, 0, 0, "rank %d assigning solution fitness\r\n", multi.meta.isle.id);
        
        p.fitness = offset_rastrigin(p.input, multi.offsets);
        
        LOG(6, 0, 0, "rank %d assigned solution fitness = %f\r\n", multi.meta.isle.id, p.fitness);
        
        //o.population.push_back(p);
        
        LOG(6, 0, 0, "rank %d adding solution %d to population\r\n", multi.meta.isle.id, i);
        
        o.population[i] = p;
        
        o.total_fitness += p.fitness;
        
        LOG(6, 0, 0, "island %d (root) initialized solution %d with fitness %f ...\r\n", multi.meta.isle.id, i, p.fitness);
        
    }
    
    LOG(4, multi.meta.isle.id, 0, "initialized objective (solution) population: total fitness = %f\r\n", multi.solutions.total_fitness);
    LOG(6, 0, 0, "%d solutions initialized ...\r\n", config::mu);
    
    std::sort(o.population.begin(), o.population.end(), compare_fitness);
    std::reverse(o.population.begin(), o.population.end());
    
    multi.eval.stats.global_best_fitness = o.population[0].fitness;
    
    LOG(6, 0, 0, "rank %d leaving solution_populate\r\n", multi.meta.isle.id);
    
    // we have our initial primary population with the calculated fitnesses
    
}

#pragma mark FUNCTION: solution_scatter()

// separate the single full population from the root process to subpopulations across all processes ...

void solution_scatter(ea &multi, objective<solution> &o) {
    
    LOG(6, 0, 0, "rank %d entered solution_scatter\r\n", multi.meta.isle.id);
    
    multi.meta.isle.population.clear();
    multi.meta.isle.population.resize(multi.meta.island_size);
    
    LOG(4, multi.meta.isle.id, 0, "scattering population root size = %lu mem 0 = %f...\r\n", multi.meta.isle.population.size(), o.population[0].fitness);
    
    double scatter_start = MPI_Wtime();
    
    MPI_Scatter(&o.population[0], multi.meta.island_size, multi.meta.solution_type, &multi.meta.isle.population[0], multi.meta.island_size, multi.meta.solution_type, 0, multi.meta.tcomm);
    
    double scatter_end = MPI_Wtime();
    double scatter_time = scatter_end - scatter_start;
    
    LOG(4, 0, 0, "population scattered, island %d population %lu mem 0 = %f...\r\n", multi.meta.isle.id, multi.meta.isle.population.size(), multi.meta.isle.population[0].fitness);
    
    island::calculate::average_fitness(multi.meta.isle);
    
    LOG(4, 0, 0, "island %d population size %lu, average fitness: %f\r\n", multi.meta.isle.id, multi.meta.isle.population.size(), multi.meta.isle.average_fitness);
    
    multi.eval.stats.total_scatter_time += scatter_time;
    
    LOG(6, 0, 0, "rank %d leaving solution_scatter\r\n", multi.meta.isle.id);
    
}

#pragma mark FUNCTION: solution_eval()

void solutions_evolve(topology &t, ea &multi) {
   
    LOG(5, multi.meta.isle.id, 0, "beginning objective<solution> evaluation %d, topology %d\r\n", multi.solutions.eval_id, t.id);
    
    multi.eval.start = std::clock();

    LOG(8, multi.meta.isle.id, 0, "calculating island %d cpd, topology %d, eval %d\r\n", multi.meta.isle.id, t.id, multi.solutions.eval_id);
    
    island::calculate::cpd(multi.meta.isle);
    
    LOG(6, multi.meta.isle.id, 0, "performing crossover island %d, topology %d, eval %d\r\n", multi.meta.isle.id, t.id, multi.solutions.eval_id);

    // crossover function performs parent selection iteratively, creating n child solutions with fitness calculated ...
    
    std::vector<solution> children = crossover(multi.meta.isle, multi.offsets);

    LOG(5, multi.meta.isle.id, 0, "island %d created %lu child solutions, population size = %lu ... selecting %d survivors  ...\r\n", multi.meta.isle.id, children.size(), multi.meta.isle.population.size(), multi.meta.island_size);

    // select the best solutions and remove n children ...
    
    select_survivors(multi.meta.isle, children, multi.meta.island_size);

    LOG(6, multi.meta.isle.id, 0, "island migrations ...\r\n");
    
    // perform migration of best solutions using the currently applied topology ...
    
    if(multi.solutions.eval_id%1 == 0) {
        
        double migrate_start = MPI_Wtime();
        
        // issue migration imports and exports ...
        
        island::migration::send(multi.meta.isle, multi.meta.isle.tcomm);
        island::migration::receive(multi.meta.isle, multi.meta.isle.tcomm);
        
        double migrate_end = MPI_Wtime();
        double migrate_time = migrate_end - migrate_start;
        
        LOG(8, 0, 0, "migrate start = %3.10f, migrate end = %3.10f, migrate time = %3.10f\r\n", migrate_start, migrate_end, migrate_time);
        
        multi.eval.stats.total_migrate_time += migrate_time;
     
        // aggregate the migration time from each island into the topology fitness ...
        
        MPI_Reduce(&migrate_time, &t.round_fitness, 1, MPI_DOUBLE, MPI_SUM, 0, multi.meta.tcomm);
        
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

    if(multi.meta.isle.id == 0) {
    t.rounds++;
    t.total_migration_time += t.round_fitness;
    t.round_fitness = t.round_fitness * -1;
    t.fitness = (t.total_migration_time / t.rounds) * -1;
        
    }
    
    LOG(8, multi.meta.isle.id, 0, "rounds=%d, total_time=%013.10f, round_fit=%013.10f, fit=%013.10f\r\n", t.rounds, t.total_migration_time, t.round_fitness, t.fitness);
    
    // output for various intervals ...
    
    if(multi.solutions.eval_id%10 == 0) {
        LOG(6, 0, 0, "island %d average fit %f\r\n", multi.meta.isle.id, multi.meta.isle.average_fitness);
    }
    
    if(multi.solutions.eval_id%100 == 0) {
                
        LOG(4, 0, 0, "population size %lu, member = %2.10f\r\n", multi.meta.isle.population.size(), multi.meta.isle.population[0].fitness);
        
        std::sort(multi.meta.isle.population.begin(), multi.meta.isle.population.end(), compare_fitness);
        std::reverse(multi.meta.isle.population.begin(), multi.meta.isle.population.end());
        
        //for(int i=0; i<multi.meta.isle.population.size(); i++) {
        //    if(multi.meta.isle.population[i].fitness == -0.0) {
        //        LOG(4, 0, 0, "ZERO fitness found on island %d for individual %d\r\n", multi.meta.isle.id, i);
        //    }
        //}
        
//        std::sort(multi.topologies.population.begin(), multi.topologies.population.end(), compare_topo_fitness);
//        std::reverse(multi.topologies.population.begin(), multi.topologies.population.end());
    
        if(multi.meta.isle.id == 0) {
            log_fn_eval_stats(multi.solutions.population, multi.topologies.population, multi.run.id, multi.solutions.eval_id, multi.eval.stats, multi.run.stats, t);
        }
        
    }
    
    if(multi.meta.isle.id == 0 && multi.solutions.eval_id%multi.topologies.evals == 0) {
        
        //for(int i=0; i<multi.solutions.population.size(); i++) {
        //    if(multi.solutions.population[i].fitness == -0.0) {
        //        LOG(4, 0, 0, "ZERO fitness found for individual %d\r\n", i);
        //    }
        //}
        
        //log_fn_eval_stats(multi.solutions.population, multi.topologies.population, multi.run.id, multi.solutions.eval_id, multi.eval.stats, multi.run.stats, t);
        
        //return;
        
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

//ea_run run_init(int id, ea &multi) {
//
//    LOG(10, 0, 0, "island %d run %d initialization\r\n", multi.meta.isle.id, id);
//
//    multi.run.id = id;
//    multi.run.start = MPI_Wtime();
//
//    if(multi.meta.isle.id != 0) {
//        return multi.run;
//    }
//
//    multi.run.eval = eval_init(0);
//
//    LOG(4, 0, 0, "initializing rastrigin population...\r\n");
//    multi.solutions.population.clear();
//    multi.solutions.population.resize(config::mu);
//
//    LOG(4, 0, 0, "initializing topology population...\r\n");
//    multi.topologies.population.clear();
//    multi.topologies.population.resize(config::topo_mu);
//
//    LOG(4, 0, 0, "initializing run variables...\r\n");
//    multi.offsets = generate_offsets(-2.5, 2.5, .5);
//
//    return multi.run;
//
//}

//void run_init(ea &multi) {
//
//    multi.solutions.population.clear();
//    multi.meta.isle.population.clear();
//
//    multi.eval.stats.init();
//    multi.solutions.eval_id = 0;
//    multi.topologies.eval_id = 0;
//
//}

void solution_evaluate(ea &multi, objective<solution> &o, topology &t) {
    
    //LOG(2, 0, 0, "evaluating topology %d\r\n", multi.topologies.population[index].id);
    
    o.eval_id++;
    
    //if(multi.solutions.eval_id > multi.solutions.evals) { return; }
    
    solutions_evolve(t, multi);
    
}

#pragma mark FUNCTION: topology_crossover()

std::vector<topology> topology_crossover(ea &multi) {
    
    std::vector<topology> children;
    //children.resize(config::topo_lambda);
    
    //if(multi.meta.isle.id != 0) { return children; }
    
    objective<topology>::calculate::cpd(multi.topologies);
    
    topology t1;
    topology t2;
    
    for(int n = 0; n < config::topo_lambda; n++) { // loop lambda

            topology child;
        
            child.id = config::topo_mu+n * multi.solutions.eval_id;
            child.world_size = multi.meta.islands;
            child.fitness = 0.0;
            child.channel_count = 0;
            child.round_fitness = 0.0;
            child.selection_distribution = 0.0;

            child.channels.resize(multi.meta.islands);
            child.channels.clear();
            //child.channels[n].senders.clear();
            //child.channels[n].receivers.clear();

            LOG(3, multi.meta.isle.id, 0, "creating topo kids\r\n");

            topology t1 = multi.selection(parent<topology>, multi.topologies, multi.topologies.population);
            topology t2 = multi.selection(parent<topology>, multi.topologies, multi.topologies.population);

            //child.fitness = (t1.fitness + t2.fitness) / 2; // an estimated fitness, average of the parents

            // calculate an adjacency matrix for each parent's associated topology for use in
            // generating child topology ...

            LOG(3, multi.meta.isle.id, 0, "parents<topology> t1=%2.10f,t2=%2.10f ...\r\n", t1.fitness, t2.fitness);

            std::vector<std::vector<int>> m1 = topology::create::matrix(t1);
            std::vector<std::vector<int>> m2 = topology::create::matrix(t2);

            // recombine the parent adjacency matrices, initialize ...

            LOG(3, multi.meta.isle.id, 0, "recombining topology %d <-> %d ...\r\n", t1.id, t2.id);

            std::vector<std::vector<int>> child_matrix;
            child_matrix.resize(multi.meta.islands);

            int comm_count = 0;
            int rec_count[multi.meta.islands];
            int snd_count[multi.meta.islands];

            for(int i=0; i<multi.meta.islands; i++) {
                rec_count[i] = 0;
                snd_count[i] = 0;
            }

            LOG(3, multi.meta.isle.id, 0, "child<topology> %d initialized \r\n", child.id);

            // iterate row->column for each x,y element in the child matrix, and for each
            // gene and randomly choose a parent from which to assign the value ...

            LOG(3, multi.meta.isle.id, 0, "performing child<topology> %d matrix crossover ...\r\n", child.id);

            while(comm_count == 0) { // failsafe to prevent empty matrix

                for(int i=0; i<m1.size(); i++) {  // child matrix row

                    child_matrix[i].resize(multi.meta.islands);

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

                    if(rand()/(RAND_MAX+1.0) < config::mutation_rate) {

                        LOG(3, multi.meta.isle.id, 0, "mutating child<topology> %d ...\r\n", child.id);

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

            LOG(3, multi.meta.isle.id, 0, "child<topology> %d born!\r\n", child.id);

            children.push_back(child);

            LOG(3, 0, 0, "child<topology> %d created by rank %d from matrix with %lu senders and %lu receivers\r\n", child.id, multi.meta.isle.id, child.channels[n].senders.size(), child.channels[n].receivers.size());

    } // loop lambda
    
    return children;
    
}

#pragma mark FUNCTION: topology_evolve()

void topology_evolve(ea &multi) {
    
    if(multi.solutions.eval_id > multi.solutions.evals) { return; }
    
    std::vector<topology> children;
    
    if(multi.meta.isle.id == 0) {
        
        children = topology_crossover(multi);
        
        LOG(6, 0, 0, "island %d topology survival, population before: %lu, fitness before: %f, best fit: %f\r\n", multi.meta.isle.id, multi.topologies.population.size(), multi.topologies.total_fitness, multi.topologies.population[0].fitness);
        
        // truncation: add new children to the population, and then kill the weakest
        
        //multi.topologies.population.insert(multi.topologies.population.end(), children.begin(), children.end());
            
    } else {
        
        //multi.topologies.population.resize(multi.topologies.mu + config::topo_lambda);
        
        children.resize(config::topo_lambda);
        
    }

    LOG(6, 0, 0, "island %d <topology> survival, population after: %lu, fitness after: %f, best fit: %f\r\n", multi.meta.isle.id, multi.topologies.population.size(), multi.topologies.total_fitness, multi.topologies.population[0].fitness);

    for (std::vector<topology>::iterator it = children.begin(); it != children.end(); ++it) {
        
    //for(int i=0; i<children.size(); i++) {

        LOG(4, multi.meta.isle.id, 0, "rank %d (root) applying child<topology>[%d] for eval\r\n", multi.meta.isle.id, it->id);

        it->apply(multi.meta.isle, *it);

        for(multi.topologies.eval_id = 1; multi.topologies.eval_id <= multi.topologies.evals; multi.topologies.eval_id++) {

            LOG(5, 0, 0, "begin eval %d on child topology %d\r\n", multi.topologies.eval_id, it->id);

            //multi.solutions.eval_id++;

            multi.evaluate(solution_evaluate, multi.solutions, *it);

        }
        
        multi.topologies.population.push_back(*it);

    }

    if(multi.meta.isle.id == 0) {
    
        std::sort(multi.topologies.population.begin(), multi.topologies.population.end(), comp_fitness());
        std::reverse(multi.topologies.population.begin(), multi.topologies.population.end());

        multi.topologies.population.erase(multi.topologies.population.begin()+(multi.topologies.mu-1), multi.topologies.population.end());

        objective<topology>::calculate::total_fitness(multi.topologies);
        
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

void topology_evaluate(ea &multi, objective<topology> &o, topology &t) {
    
    LOG(4, multi.meta.isle.id, 0, "applying topology %d to %d islands ...\r\n", t.id, multi.meta.islands);
    
    t.fitness = 0.0;
    t.round_fitness = 0.0;
    
    // distribute the send\recv channels from the passed topology to all islands ...
    
    t.apply(multi.meta.isle, t);
    
    LOG(8, 0, 0, "evaluating topology %d\r\n", t.id);
    
    for(o.eval_id = 1; o.eval_id%(o.evals+1) != 0; o.eval_id++) {

        LOG(5, 0, 0, "begin eval %d solution objective on topology %d\r\n", o.eval_id, t.id);
        
        multi.evaluate(solution_evaluate, multi.solutions, t);
        
        
    }
    
    LOG(8, multi.meta.isle.id, 0, "evaluated topology id=%d rounds=%d fitness=%3.10f\r\n", t.id, t.rounds, t.fitness);

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
    
    // track derived mpi datatypes used to transfer custom structs
    
//    channel chan;
//
//    MPI_Datatype channel_types[4] = { MPI_INT, MPI_DOUBLE, MPI_INT, MPI_INT };
//    int channel_blocks[4] = { 2, (int)chan.senders.size(), (int)chan.receivers.size(), 1 };
//
//    MPI_Aint channel_offsets[4];
//
//    MPI_Get_address(&chan.id,&channel_offsets[0]);
//    MPI_Get_address(&chan.fitness,&channel_offsets[1]);
//    MPI_Get_address(&chan.senders[0],&channel_offsets[2]);
//    MPI_Get_address(&chan.receivers[0],&channel_offsets[3]);
//
//    MPI_Type_create_struct(5, channel_blocks, channel_offsets, channel_types, &multi.meta.MPI_CHANNEL);
//    MPI_Type_commit(&multi.meta.MPI_CHANNEL);
//
//    topology t;
//
//    MPI_Datatype topology_types[3] = { MPI_INT, MPI_DOUBLE, multi.meta.MPI_CHANNEL };
//    int topology_blocks[3] = { 4, 3, (int)t.channels.size() };
//
//    MPI_Aint topology_offsets[3];
//
//    MPI_Get_address(&t.id,&topology_offsets[0]);
//    MPI_Get_address(&t.fitness,&topology_offsets[1]);
//    MPI_Get_address(&t.channels[0],&topology_offsets[2]);
//
//    MPI_Type_create_struct(3, topology_blocks, topology_offsets, topology_types, &multi.meta.MPI_TOPOLOGY);
//    MPI_Type_commit(&multi.meta.MPI_TOPOLOGY);

    // ----- end derived mpi datatypes
    
    multi.meta.island_size = config::mu / multi.meta.islands;
    
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
    
    multi.topologies.population.resize(config::mu);
    
    double total_init_duration = 0.0;

    multi.run.stats.init_duration = ( std::clock() - multi.meta.start ) / (double) CLOCKS_PER_SEC;
    
    multi.offsets = generate_offsets(-2.5, 2.5, .5);
    
    //MPI_Gather(&multi.run.stats.init_duration, 1, MPI_DOUBLE, &multi.run.stats.init_duration, 1, MPI_DOUBLE, 0, multi.meta.isle.tcomm);

    
    return multi;
    
}


#endif /* ea_h */

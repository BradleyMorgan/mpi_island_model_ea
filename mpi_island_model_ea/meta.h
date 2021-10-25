//
//  meta.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 10/17/21.
//  Copyright © 2021 Bradley Morgan. All rights reserved.
//

#ifndef meta_h
#define meta_h

#pragma mark FUNCTION: topology_crossover()

// specialized recombination operator for the @topology{} datatype

std::vector<topology> topology_crossover(meta &meta) {
   
    std::vector<topology> children;
    
    // return an empty vector to any non-root rank as a
    // placeholder for future use
    
    if(meta.variant.isle.id != 0) {
        children.resize(meta.topologies.lambda);
        return children;
    }
        
    meta.topologies.cpd();
    
    for(int n = 0; n < meta.topologies.lambda; n++) { // loop lambda

        int test = meta.topologies.mu + ((n * meta.topologies.run.id) + meta.topologies.run.eval.id);
        
        printf("%d + (%d * %d + %d) = %d\r\n", meta.topologies.mu, n, meta.topologies.run.id, meta.topologies.run.eval.id, test);
       
        // create child skeleton ...
        
        topology child;
        
        child.id = test;
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

void topology_evolve(solver &solver, meta &meta) {
    
    LOG(2, solver.variant.isle.id, 0, "\r\n\r\n --- ISLAND %d META EVOLUTION [%d,%d] SOLVER [%d,%d]  --- \r\n\r\n", meta.variant.isle.id, meta.topologies.run.id, meta.topologies.run.eval.id, solver.solutions.run.id, solver.solutions.run.eval.id);
    
    LOG(6, 0, 0, "\r\n\r\n --- ISLAND %d META EVOLUTION [%d,%d] SOLVER [%d,%d]  --- \r\n\r\n", meta.variant.isle.id, meta.topologies.run.id, meta.topologies.run.eval.id, solver.solutions.run.id, solver.solutions.run.eval.id);
    
    std::vector<topology> children;
    
    if(meta.variant.isle.id == 0) {
        
        // the root island holds the only valid, authoritative topology population,
        // we only need to perform parent selection and recombination on that process ...
        
        children = topology_crossover(meta);
        
        LOG(6, 0, 0, "island %d topology survival, population before: %lu, fitness before: %f, best fit: %f\r\n", meta.variant.isle.id, meta.topologies.population.size(), meta.topologies.aggregate.value.fitness, meta.topologies.population[0].fitness);
            
    } else {
        
        // ensure that the secondary islands have memory allocated to avoid potential segfault ....
        
        children.resize(meta.topologies.lambda);
        
        //return;
        
    }

    LOG(6, 0, 0, "island %d <topology> survival, population after: %lu, fitness after: %f, best fit: %f\r\n", meta.variant.isle.id, meta.topologies.population.size(), meta.topologies.aggregate.value.fitness, meta.topologies.population[0].fitness);

    // in order to determine a topoology fitness, perform evolutionary cycles of the solver population,
    // with all islands participating using the send\recv channels as defined in the topology
    // and then assign the aggregate time spent in migration as topology fitness ...
    
    for (std::vector<topology>::iterator it = children.begin(); it != children.end(); ++it) {

        LOG(4, meta.variant.isle.id, 0, "ISLAND %d (root) applying child<topology>[%d] for eval\r\n", meta.variant.isle.id, it->id);

        // the apply function distributes a topology's send\recv channels to all islands ...

        it->apply(meta.variant.isle, *it);

        // once the send\recv channels have been established, we perform a user defined number of evolutionary cycles
        // of the solution population ...
        
        //for(meta.topologies.run.eval.id = 1; meta.topologies.run.eval.id <= meta.topologies.max_fit_evals; meta.topologies.run.eval.id++) {
            
            meta.topologies.begin(meta.topologies.run.eval, meta);
      
            solver_begin(meta, solver, *it);
                
            meta.topologies.end(meta.topologies.run.eval, meta);
            
        //}
        
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


#pragma mark FUNCTION: topology_populate()

// returns a collection of randomly generated adjaceny matrices, representing                           |
// an island (communication) topology.  accepts a reference to a list of island                         |
// (process) identifiers from @ea{@meta{@island_ids[param:0}} to use as indices.                        |

void topologies_populate(meta &meta) {
    
    meta.topologies.population.clear();
    
    LOG(6, meta.variant.isle.id, 0, "ISLAND %d OBJECTIVE %d (root) initializing topology population ... \r\n", meta.variant.isle.id, meta.topologies.id);

    for(int i=1; i <= meta.topologies.mu; i++) {
        
        topology t;
        t.id = i;
        
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
    
    LOG(2, meta.variant.isle.id, 0, "initialized objective (topology) population size %lu, \r\n", meta.topologies.population.size());
    
}

void benchmark_topology(meta &meta) {
    
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

void meta_begin(meta &meta, solver &solver) {
    
    // number gens and num of offspring
    
    // solver evals = total toplogy evals * solver_runs * (solver_mu + num_solver_gens * solver_lambda) = mu + lambda * gens
    // total topo evals = topology mu + (topology generations * topology lambda)
    
    for(meta.topologies.run.id = 1; meta.topologies.run.id <= meta.topologies.max_runs; meta.topologies.run.id++) {
        
        // solver ea 5 runs ... each time it runs n gens which depends on (for a max_eval limit ...
        // num solver gens = (max_solver_evals - solver_mu) / solver_lambda
        // look at eval vs. fitness solver ea graphs
        // k gens -> fitness
        
        meta.topologies.begin(meta.topologies.run, meta);
        meta.topologies.populate(meta, topologies_populate);
        
        for(int i=0; i<meta.topologies.mu; i++) {
            
            //for(meta.topologies.run.eval.id = 1; meta.topologies.run.eval.id <= meta.topologies.max_fit_evals; meta.topologies.run.eval.id++) {
                
                meta.topologies.begin(meta.topologies.run.eval, meta);
                
                solver_begin(meta, solver, meta.topologies.population[i]);
                
                meta.topologies.end(meta.topologies.run.eval, meta);
                
            //}

        }
        
        LOG(6, 0, 0, "ISLAND %d FINISHED objective<topology> POPULATION AT SOLVER[%d,%d] META[%d,%d]\r\n", solver.variant.isle.id, solver.solutions.run.id, solver.solutions.run.eval.id, meta.topologies.run.id, meta.topologies.run.eval.id);

        // generations or cycles, what is the termination?  max_evals / num_offspring
        
        int generation = 0;
        int meta_evolution_termination = 100;
        
        for(generation = 1; generation <= meta_evolution_termination; generation++) {
            
            //meta.topologies.begin(meta.topologies.run.eval, meta);
            
            meta.topologies.evolve(solver, meta, topology_evolve);
            
           // meta.topologies.end(meta.topologies.run.eval, meta);

        }
        
        meta.topologies.end(meta.topologies.run, meta);
        
    }
    
}


#endif /* meta_h */

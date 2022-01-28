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

std::vector<topology> topology_crossover(ea_meta &meta) {
   
    std::vector<topology> children;
    
    // return an empty vector to any non-root rank as a
    // placeholder for future use
    
    if(mpi.id != 0) {
        children.resize(meta.topologies.lambda);
        return children;
    }
        
    //meta.topologies.cpd();
    
    int total_channels = 0;
    std::vector<int> channel_counts;
    
    for(int n = 0; n < meta.topologies.lambda; n++) { // loop lambda

        int child_id = meta.topologies.mu + (meta.topologies.lambda * meta.topologies.run.cycle.id + n + 1);
        
        LOG(6, mpi.id, 0, "%d * %d + %d + 1 = %d\r\n", meta.topologies.mu, meta.topologies.run.cycle.id, n, child_id);
       
        // create child skeleton ...
        
        topology child;
        
        child.init(child_id);
        child.init_multi();
        
        LOG(5, mpi.id, 0, "creating topo kids\r\n");

        std::pair<topology,topology> t1 = meta.topologies.tournament(binary_tournament<topology>);
        std::pair<topology,topology> t2 = meta.topologies.tournament(binary_tournament<topology>);
        
        // topology t1 = meta.topologies.select(parent<topology>, 2);
        // topology t2 = meta.topologies.select(parent<topology>);

        // calculate an adjacency matrix for each parent's associated topology for use in
        // generating child topology ...

        LOG(5, mpi.id, 0, "parents<topology> t1=%2.10f,t2=%2.10f ...\r\n", t1.first.fitness, t2.first.fitness);

        std::vector<std::vector<int>> m1 = topology::create::matrix(t1.first);
        std::vector<std::vector<int>> m2 = topology::create::matrix(t2.first);

        // recombine the parent adjacency matrices, initialize ...

        LOG(5, mpi.id, 0, "recombining topology %d <-> %d ...\r\n", t1.first.id, t2.first.id);

        std::vector<std::vector<int>> child_matrix;
        child_matrix.resize(mpi.size);

        int comm_count = 0;
        // int send_max = mpi.size * ((config::send_cap * 1.0) / 100.0);
        // int recv_max = mpi.size * ((config::recv_cap * 1.0) / 100.0);
        int send_max = config::send_cap;
        int recv_max = config::recv_cap;
        int rec_count[mpi.size];
        int snd_count[mpi.size];

        for(int i=0; i<mpi.size; i++) {
            rec_count[i] = 0;
            snd_count[i] = 0;
        }

        LOG(5, mpi.id, 0, "child<topology> %d initialized \r\n", child.id);

        // iterate row->column for each x,y element in the child matrix, and for each
        // gene and randomly choose a parent from which to assign the value ...

        LOG(5, mpi.id, 0, "performing child<topology> %d matrix crossover ...\r\n", child.id);

        while(comm_count == 0) { // failsafe to prevent empty matrix

            for(int i=0; i<m1.size(); i++) {  // child matrix row

                child_matrix[i].resize(mpi.size);

                for(int j=0; j<m2.size(); j++) { // child matrix column

                    //if(rec_count[j] >= config::migration_cap) {
                    if(rec_count[j] >= recv_max) {
                        LOG(6, 0, 0, "migration cap limit reached for process %d\r\n", i);
                        continue;
                    }

                    //if(snd_count[i] >= config::send_cap) {
                    if(snd_count[i] >= send_max) {
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

                //if(rec_count[j] >= config::migration_cap) {
                if(rec_count[j] >= recv_max) {
                    LOG(6, 0, 0, "migration cap limit reached for process %d\r\n", i);
                    continue;
                }

                //if(snd_count[i] >= config::send_cap) {
                if(snd_count[i] >= send_max) {
                    LOG(6, 0, 0, "send cap limit reached for process %d\r\n", i);
                    continue;
                }

                if(rand()/(RAND_MAX+1.0) < meta.topologies.mutation_rate) {

                    LOG(2, mpi.id, 0, "mutating child<topology> %d ...\r\n", child.id);

                    if(child_matrix[i][j] == 0 && rand()/(RAND_MAX+1.0) < config::sparsity) {

                        LOG(2, mpi.id, 0, "mutation adding channel at %d,%d ...\r\n", i, j);
                        
                        child_matrix[i][j] = 1;
                        rec_count[j]++;
                        snd_count[i]++;
                        comm_count++;

                    }

                    if(child_matrix[i][j] == 1 && rand()/(RAND_MAX+1.0) > config::sparsity && comm_count > 1) {
                                                
                        if(rand()%2 == 1) {
                            LOG(2, mpi.id, 0, "mutation removing channel at %d,%d ...\r\n", i, j);
                            child_matrix[i][j] = 0;
                            comm_count--;
                        } else {
                            int randj, randi;
                            do {
                                randi = rand()%mpi.size;
                                randj = rand()%mpi.size;
                            }
                            while(child_matrix[randi][randj] == 1);
                            LOG(2, mpi.id, 0, "mutation swapping channel at %d,%d with %d,%d ...\r\n", i, j, randi, randj);
                            child_matrix[i][j] = 0;
                            child_matrix[randi][randj] = 1;
                        }

                    }

                }

            } // matrix col

        } // matrix row
        
        // convert the adjaceny matrix to a sender and receiver arrays for use with MPI send\recv ...
    
        topology::create::channels(child, child_matrix);

        LOG(5, mpi.id, 0, "child<topology> %d born!\r\n", child.id);

        children.push_back(child);

        channel_counts.push_back(child.stats.total_channels);
        
        LOG(5, 0, 0, "child<topology> %d created by rank %d from matrix with %lu senders and %lu receivers\r\n", child.id, mpi.id, child.channels[n].senders.size(), child.channels[n].receivers.size());

    } // loop lambda
    
    if(mpi.id == 0) {
    
        double mean = total_channels / children.size();
        double variance = 0;
        
        for(int n = 0; n<meta.topologies.lambda; n++) {
          variance += (channel_counts[n] - mean) * (channel_counts[n] - mean);
        }
        
        variance /= meta.topologies.lambda;
        
        double std_dev = sqrt(variance);
        
        LOG(5, 0, 0, "topology generation %d size=%lu | channels: mean=%f variance=%f standard_deviation=%f\r\n", meta.topologies.run.cycle.id, children.size(), mean, variance, std_dev);
            
    }
    
    return children;
    
}

#pragma mark FUNCTION: topology_evolve()

void topology_evolve(ea_solver &solver, ea_meta &meta) {
       
    std::vector<topology> children;
    
    if(mpi.id == 0) {
        
        // the root island holds the only valid, authoritative topology population,
        // we only need to perform parent selection and recombination on that process ...
        
        children = topology_crossover(meta);
        
        LOG(6, 0, 0, "island %d topology survival, population before: %lu, fitness before: %f, best fit: %f\r\n", mpi.id, meta.topologies.population.size(), meta.topologies.aggregate.value.fitness, meta.topologies.population[0].fitness);
            
    } else {
        
        children.resize(meta.topologies.lambda);
        
    }

    LOG(6, 0, 0, "island %d <topology> survival, population after: %lu, fitness after: %f, best fit: %f\r\n", mpi.id, meta.topologies.population.size(), meta.topologies.aggregate.value.fitness, meta.topologies.population[0].fitness);

    // in order to determine a topoology fitness, perform evolutionary cycles of the solver population,
    // with all islands participating using the send\recv channels as defined in the topology
    // and then assign the aggregate time spent in migration as topology fitness ...
    
    // 𝛭𝜆 iterations ...
    
    for (std::vector<topology>::iterator it = children.begin(); it != children.end(); ++it) {

        // solver will execute topologies_lambda times, e.g. 32
        
        LOG(4, mpi.id, 0, "ISLAND %d (root) applying child<topology>[%d] for eval\r\n", mpi.id, it->id);

        // the apply function distributes a topology's send\recv channels to all islands ...

        it->apply(meta.model.isle, *it);

        // once the send\recv channels have been established, we perform a user defined number of evolutionary cycles
        // of the solution population ...
        
        // 𝑆𝑟𝑚𝑎𝑥 * 𝑆𝑒𝑚𝑎𝑥 iterations in solver_begin ...
        
        meta.ea::begin(meta.topologies, meta.topologies.run.cycle.eval, &(*it), solver.solutions.run);
      
        solver_begin(meta, solver, *it);
                
        meta.ea::end(meta.topologies, meta.topologies.run.cycle.eval, &(*it), solver.solutions.run);
        
        meta.topologies.population.push_back(*it);

    }
    
    // children evaluated for fitness, perform survival selection ...
    
    if(mpi.id == 0) {
    
        // similar to crossover, we only need to perform survival selection on the root island,
        // so only truncate the the worst individuals in the root island population ...
        
        LOG(2, 0, 0, "defining fronts, deriving from front[0..n].dom_genomes ...\r\n");
        
        std::vector<std::vector<topology*>> fronts = meta.topologies.define_fronts();
        
        LOG(2, 0, 0, "defined %lu fronts ...\r\n", fronts.size());
        
        LOG(2, 0, 0, "assigning crowding distances for fronts ...\r\n");
        
        meta.topologies.crowding_distance(fronts);
        
        LOG(2, 0, 0, "assigned crowding distances for fronts ...\r\n");
        
        LOG(2, mpi.id, 0, "FRONT %d: distance=%f, x=%f, y=%f\r\n", fronts[0][0]->id, fronts[0][0]->distance, fronts[0][0]->fitness, fronts[0][0]->fitness);
        
        meta.log_fronts(fronts);
        
        LOG(2, 0, 0, "sorting by multiple objective fitness i[0] = %f, i[mu] = %f...\r\n", meta.topologies.population[0].fitness, meta.topologies.population[meta.topologies.mu-1].fitness);

        std::sort(meta.topologies.population.begin(), meta.topologies.population.end(), compare_multi<topology>);

        LOG(2, 0, 0, "sorted by multiple objective fitness i[0] = %f, i[mu] = %f...\r\n", meta.topologies.population[0].fitness, meta.topologies.population[meta.topologies.mu-1].fitness);
        
        LOG(2, 0, 0, "truncating population size %lu...\r\n", meta.topologies.population.size());
        
        meta.topologies.population.erase(meta.topologies.population.begin()+meta.topologies.mu, meta.topologies.population.end());

        LOG(2, 0, 0, "truncated population size %lu...\r\n", meta.topologies.population.size());
        
        LOG(2, 0, 0, "calculating aggregate fitness current = %f...\r\n", meta.topologies.aggregate.value.fitness);
        
        meta.topologies.fitness();
        
        LOG(2, 0, 0, "calculated aggregate fitness new = %f...\r\n", meta.topologies.aggregate.value.fitness);
        
    }
   
}

#pragma mark FUNCTION: topology_populate()

// returns a collection of randomly generated adjaceny matrices, representing                           |
// an island (communication) topology.  accepts a reference to a list of island                         |
// (process) identifiers from @ea{@meta{@island_ids[param:0}} to use as indices.                        |

void topologies_populate(ea_meta &meta) {
    
//    int send_max = mpi.size * ((config::send_cap * 1.0) / 100.0);
//    int recv_max = mpi.size * ((config::recv_cap * 1.0) / 100.0);
    int send_max = config::send_cap;
    int recv_max = config::recv_cap;
    
    int matrix_size = mpi.size * mpi.size;
    
    LOG(2, 0, mpi.id, "sparsity=%f send_max=%d recv_max=%d prob=%f\r\n", config::sparsity, send_max, recv_max, config::sparsity * matrix_size);
    
    meta.topologies.population.clear();
    
    LOG(6, mpi.id, 0, "ISLAND %d OBJECTIVE %d (root) initializing topology population ... \r\n", mpi.id, meta.topologies.id);

    std::vector<int> channel_counts;
    
    for(int i=1; i <= meta.topologies.mu; i++) {
        
        topology t;
        t.id = i;
        
        if(mpi.id != 0) {
            LOG(8, 0, 0, "island %d (leaf) adding topology stub %d, returning ... \r\n", mpi.id, i);
            meta.topologies.population.push_back(t);
        } else {
            LOG(8, 0, 0, "island %d (root) topology %d\r\n", mpi.id, i);
            topology::create::dynamic(t);
            LOG(2, 0, 0, "dynamic topology %d created with %d channels\r\n", i, t.stats.total_channels);
            meta.topologies.population.push_back(t);
            meta.topologies.run.stats.total_channels += t.stats.total_channels;
            channel_counts.push_back(t.stats.total_channels);
        }
        
    }
    
    if(mpi.id == 0) {
    
        double mean = meta.topologies.run.stats.total_channels / meta.topologies.mu;
        double variance = 0;
        
        for(int n = 0; n<meta.topologies.mu; n++) {
          variance += (channel_counts[n] - mean) * (channel_counts[n] - mean);
        }
        
        variance /= meta.topologies.mu;
        
        double std_dev = sqrt(variance);
        
        LOG(2, 0, 0, "created %lu topologies | channels: mean=%f var=%f stdev=%f\r\n", meta.topologies.population.size(), mean, variance, std_dev);
        
    }
    
    // we have our initial topology population, *NOT* evaluated for fitness
    
    LOG(5, mpi.id, 0, "initialized objective (topology) population size %lu \r\n", meta.topologies.population.size());
    
}

void benchmark_topology(ea_meta &meta) {
    
    meta.topologies.population.clear();
    
    topology t;
    
    t.stats = {};
    t.fitness = 0.0;
    t.selection_distribution = 0.0;
    t.channels = {};
    t.channels.resize(mpi.size);
    
    for(int i=0; i<mpi.size; i++) {
        
        t.channels[i].senders = {};
        t.channels[i].receivers = {};
        
        int next = i+1 < mpi.size ? i+1 : 0;
        int prev = i-1 < 0 ? (int)mpi.size-1 : i-1;

        t.channels[i].senders.push_back(prev);
        t.stats.send_channels++;
        t.stats.total_channels++;
        t.channels[i].receivers.push_back(next);
        t.stats.recv_channels++;
        t.stats.total_channels++;
        
    }

    meta.topologies.population.push_back(t);

}

template<> template<typename i> void objective<topology>::end(i &interval, topology *current) {
    
    interval.end(current);
    
    this->log_end(interval);
    
    if(interval.log_population_interval != 0 && interval.id%interval.log_population_interval == 0) {
        this->log_population(interval);
    }
    
}

template<typename e> void ea_meta::begin(e &target) {
    
    // 𝛭𝑟𝑚𝑎𝑥 iterations ...
    
    for(this->topologies.run.id = 1; this->topologies.run.id <= this->topologies.run.max; this->topologies.run.id++) {
    
        this->ea::begin(this->topologies, &this->topologies.population[0]);
        
        // evaluate initial population by applying each randomly generated
        // channels accordingly and return the elapsed time to perform
        // all island migrations as the topology fitness
        
        if(config::ea_mode == 1) {
       
            this->topologies.populate(*this, topologies_populate);
            
            // 𝛭𝜇 iterations ...
            
            for(int i=0; i<this->topologies.mu; i++) {
                
                // 𝑆𝑟𝑚𝑎𝑥 * 𝑆𝑒𝑚𝑎𝑥 iterations in solver_begin  ...
                
                this->ea::begin(this->topologies, this->topologies.run.cycle.eval, &this->topologies.population[i], target.solutions.run);
                
                solver_begin(*this, target, this->topologies.population[i]);
                
                this->ea::end(this->topologies, this->topologies.run.cycle.eval, &this->topologies.population[i], target.solutions.run);
                
            } // 𝛭𝜇 * 𝑆𝑟𝑚𝑎𝑥 * 𝑆𝑒𝑚𝑎𝑥 iterations

            // evaluate initial population for multi-objective fitness metrics
            
            std::vector<std::vector<topology*>> fronts = this->topologies.define_fronts();
            
            this->topologies.crowding_distance(fronts);
            
            // 𝛭𝑒𝑚𝑎𝑥 iterations in solver_begin  ...
            
            for(this->topologies.run.cycle.id = 1; this->topologies.run.cycle.id <= this->topologies.run.cycle.max; this->topologies.run.cycle.id++) {
               
                // 𝛭𝜆 * 𝑆𝑟𝑚𝑎𝑥 * 𝑆𝑒𝑚𝑎𝑥 iterations in topology_evolve
                
                this->ea::begin(this->topologies, this->topologies.run.cycle, this->topologies.run.local);
                
                this->topologies.evolve(target, *this, topology_evolve);
                
                this->ea::end(this->topologies, this->topologies.run.cycle);
                
            }  // 𝛭𝑒𝑚𝑎𝑥 * 𝛭𝜆 * 𝑆𝑟𝑚𝑎𝑥 * 𝑆𝑒𝑚𝑎𝑥 iterations
        
            // (𝛭𝜇 * 𝑆𝑟𝑚𝑎𝑥 * 𝑆𝑒𝑚𝑎𝑥) + (𝛭𝑒𝑚𝑎𝑥 * 𝛭𝜆 * 𝑆𝑟𝑚𝑎𝑥 * 𝑆𝑒𝑚𝑎𝑥) iterations
            
            this->ea::end(this->topologies, this->topologies.run, this->topologies.run.local);
            
        } else {

            //benchmark_topology(*this);
            
            // 𝛭𝑒𝑚𝑎𝑥 iterations ...
            
            for(this->topologies.run.cycle.id = 1; this->topologies.run.cycle.id <= this->topologies.run.cycle.max; this->topologies.run.cycle.id++) {
             
                this->ea::begin(this->topologies, this->topologies.run.cycle, this->topologies.run.cycle.local);
                
                // we don't need to evolve since the benchmark uses a static topology (id=1)
                // but, simulate the same iteration logic so that logging is more consistent
                // between the benchmark and experimental output
                
                for(int i=0; i<this->topologies.lambda; i++) {
                  
                    // 𝑆𝑟𝑚𝑎𝑥 * 𝑆𝑒𝑚𝑎𝑥 iterations in solver_begin  ...
                    benchmark_topology(*this);
                    
                    this->ea::begin(this->topologies, this->topologies.run.cycle.eval, &this->topologies.population[0], target.solutions.run);
                   
                    solver_begin(*this, target, this->topologies.population[0]);
                    
                    this->ea::end(this->topologies, this->topologies.run.cycle.eval, this->topologies.run.cycle.eval.local, target.solutions.run);
                    
                }
                
                if(mpi.id == 0) {
                
                    std::fprintf(config::ea_2_multi_out, "%d," "%d," "%d," "%s," "%s," "%s," "%s," "%f," "%f\r\n", this->topologies.run.id, this->topologies.run.cycle.id, this->topologies.run.cycle.eval.id, "NA", "NA", "NA", "NA", topologies.population[0].fitness_multi.first, topologies.population[0].fitness_multi.second);
                    
                    fflush(config::ea_2_multi_out);
                    
                }
                
                this->ea::end(this->topologies, this->topologies.run.cycle, this->topologies.run.cycle.local);
                
            }
            
        }
        
        this->ea::end(this->topologies);
        
    }
    
    // 𝛭𝑟𝑚𝑎𝑥 * ((𝛭𝜇 * 𝑆𝑟𝑚𝑎𝑥 * 𝑆𝑒𝑚𝑎𝑥) + (𝛭𝑒𝑚𝑎𝑥 * 𝛭𝜆 * 𝑆𝑟𝑚𝑎𝑥 * 𝑆𝑒𝑚𝑎𝑥)) iterations
    
}


#endif /* meta_h */

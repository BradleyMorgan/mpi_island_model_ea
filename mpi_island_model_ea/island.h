//
//  island.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 2/26/20.
//  Copyright © 2020 Bradley Morgan. All rights reserved.
//

#ifndef island_h
#define island_h

#include "config.h"
#include "utility.h"

struct individual {
    
    std::array<double, DIM> input;
    
    double fitness;
    double selection_distribution;
    
};

struct island {
    
    int id;
    
    double total_fitness;
    
    std::vector<individual> population;
    std::vector<double> cpd;
    std::vector<int> senders;
    std::vector<int> receivers;
    
    // calculate the island's total fitness for distribution ...
    
    void calc_total_fitness() {

        this->total_fitness = 0;
        
        std::vector<individual>::iterator it;
        
        for(it = this->population.begin(); it != this->population.end(); ++it) {
            this->total_fitness += it->fitness;
        }
        
    }
    
    // this function calculates the cumulative probability distribution to be used by
    // the fitness proportional (roulette wheel) selection ...
    
    void calc_cpd() {
        
        double cumulative_probability = 0.0;
        
        this->calc_total_fitness();
        this->cpd.clear();
        
        for(int i=0; i<this->population.size(); i++) {

            this->population[i].selection_distribution = (double)this->population[i].fitness / this->total_fitness;

            cumulative_probability += this->population[i].selection_distribution;
            this->cpd.push_back(cumulative_probability);
            
        }
        
    }
    
    // create a receive channel for every island in this island's senders list ...
    
    void receive_migrant(MPI_Comm &comm) {
        
        std::array<double, DIM> x;
        
        MPI_Status migrant_status;
        
        for(int i=0; i<this->senders.size(); i++) {
            LOG(10, 0, 0, "island %d waiting for migrant from island %d ... \r\n", this->id, this->senders[i]);
            MPI_Recv(&x, DIM, MPI_DOUBLE, this->senders[i], 0, comm, &migrant_status);
            this->population[rand()%population.size()].input = x;
            LOG(10, 0, 0, "island %d received migrant from island %d: [%f,%f] with status %d\r\n", this->id, migrant_status.MPI_SOURCE, this->population[0].input[0], this->population[0].input[0], migrant_status.MPI_ERROR);
        }
        
    }
    
    // create a send channel for every island in this island's receivers list ...
    
    void send_migrant(MPI_Comm &comm) {
        
        for(int i=0; i<this->receivers.size(); i++) {
            LOG(10, 0, 0, "island %d sending migrant to island %d ... \r\n", this->id, this->receivers[i]);
            MPI_Send(&this->population[i].input, DIM, MPI_DOUBLE, this->receivers[i], 0, comm);
            LOG(10, 0, 0, "island %d sent migrant to island %d: [%f,%f]\r\n", this->id, this->receivers[i], this->population[i].input[0], this->population[i].input[1]);
        }
                
    }
    
    double average_fitness() {
        
        this->calc_total_fitness();
        
        return this->total_fitness / this->population.size();
        
    }
    
};

// struct to hold a node's (or island's) senders and receivers

struct group {
    
    int node;
    
    double fitness = 0.0;
    
    std::vector<int> senders;
    std::vector<int> receivers;
    
};

// struct to hold the full island topology using a vector of the individual group objects

struct topology {
  
    int rounds = 0;
    
    double fitness = 0.0;
    double round_fitness = 0.0;
    double selection_distribution;
    
    std::vector<group> comm;
    
};

// from a provided adjacency matrix, create sender and receiver arrays for
// each island for easier use with MPI send and receive ...

std::vector<group> create_group(std::vector<std::vector<int>> &matrix) {
    
    LOG(10, 0, 0, "creating group...\r\n");
    
    std::vector<group> comm;
    
    comm.resize(matrix.size());
    
    // iterate through the adjacency matrix and assign any marked channels
    // to corresponding island sender or receiver array ...
    
    for(int i=0; i<matrix.size(); i++) {  // matrix row
        
        for(int j=0; j<matrix[i].size(); j++) { // matrix column
        
            LOG(10, 0, 0, "%d ", matrix[i][j]);
            
            if(matrix[i][j] == 1) { // i sends to j ...
                
                LOG(10, 0, 0, "adding receiver %d to %d and sender %d to %d\r\n", j, i, i, j);
                
                comm[i].receivers.push_back(j); // add island j to island i receivers list
                comm[j].senders.push_back(i); // add island i to island j senders list
                
            }
            
        }
        
        LOG(10, 0, 0, "\r\n");
        
    }
    
    return comm;
    
}

// with a provided vector of topology groups (which holds an island's senders and receivers),
// create an adjacency matrix of size world_size*world_size, populating the coordinate
// index with 1 whenever a receiver is found ...
//
//    1  2  3  4
// 1 [0][1][0][0]  ->  1 sends to 2
// 2 [0][0][1][0]  ->  2 sends to 3
// 3 [0][0][0][0]  ->  3 sends to nobody
// 4 [1][0][0][0]  ->  4 sends to 1

std::vector<std::vector<int>> create_adjaceny_matrix(std::vector<group> &comm, int world_size) {
    
    // the passed vector comm should contain world_size group objects,
    // tracking the senders and receivers for each island
    
    std::vector<std::vector<int>> matrix;
    matrix.resize(world_size);

    for(int i=0; i<world_size; i++) { // matrix row, sized with world_size indices

        matrix[i].resize(world_size);
        
        for(int j=0; j<world_size; j++) { // matrix column [0..world size]
            
            // here we only mark the receiving nodes for each island row with a 1
            // a boolean value of true within a a row's index will tell us what we need to know
            // to track the communication channel on both sides
            
            if(std::find(comm[i].receivers.begin(), comm[i].receivers.end(), j) != comm[i].receivers.end()) {
                // if island i's receivers array contains island j, mark that as a sender\receiver channel ...
                matrix[i][j] = 1;
            } else {
                // otherwise, no sender\receiver channel exists, so 0 ...
                matrix[i][j] = 0;
            }
            
        }
        
    }
    
    return matrix;
    
}

bool prob_true(double p){
    return rand()/(RAND_MAX+1.0) < p;
}

// create a randomly populated adjacency matrix of size world_size*world_size,
// communication neighbors with probability determined by the sparsity parameter ...

std::vector<std::vector<int>> create_dyn_adjaceny_matrix(int world_size) {
    
    std::vector<std::vector<int>> matrix;
    matrix.resize(world_size);

    int comm_count = 0; // failsafe for avoiding an empty matrix if the sparsity probability is low
    int rec_count[world_size];
    
    while(comm_count == 0) {
    
        for(int i=0; i<world_size; i++) {

            matrix[i].resize(world_size);
            
            for(int j=0; j<world_size; j++) {
                
                if(rec_count[j] > config::migration_cap) {
                    LOG(4, 0, 0, "migration cap limit reached for process %d\r\n", i);
                    continue;
                }
                
                if(prob_true(config::sparsity) && i != j) {
                    matrix[i][j] = 1;
                    comm_count++;
                    rec_count[j]++;
                } else {
                    matrix[i][j] = 0;
                }

                
            }
            
        }
        
    }
    
    
    return matrix;
    
}

// calculate the cumulative probability distribution for the topology population ...

std::vector<double> topo_cpd(std::vector<topology> &topologies) {
    
    LOG(10, 0, 0, "generating cpd\r\n");
           
    double total_fitness = 0.0;
    double cumulative_probability = 0.0;
    
    std::vector<double> cpd;
    
    for(int i=0; i<topologies.size(); i++) {
    
        total_fitness += topologies[i].fitness;
        
    }
    
    std::reverse(topologies.begin(), topologies.end());
    
    for(int i=0; i<topologies.size(); i++) {

        topologies[i].selection_distribution = (double)topologies[i].fitness / total_fitness;

        cumulative_probability += topologies[i].selection_distribution;
        cpd.push_back(cumulative_probability);
        
    }
    
    return cpd;
    
}

// parent selection from the topology population, based on the cumulative probability distribution ...

topology select_topo_parent(std::vector<double> &cpd, std::vector<topology> &topologies) {
    
    LOG(7, 0, 0, "selecting topo parent: ");
    
    topology t;
    
    int i = 1;
    
    double r = ((double)rand()/(double)RAND_MAX);
    
    while (cpd[i] < r ) { i++; }
    
    t = topologies[i];
    
    LOG(7, 0, 0, "%d\r\n", i);
    
    return t;
    
}

// create a new generation of topologies using parent selection ...

std::vector<topology> topo_gen(std::vector<topology> &topologies, int world_size) {
    
    std::vector<double> cpd = topo_cpd(topologies);
    
    std::reverse(topologies.begin(), topologies.end());
    
    std::vector<topology> children;

    for(int n = 0; n < config::topo_lambda; n++) {

        LOG(10, 0, 0, "creating topo kids\r\n");

        topology t1 = select_topo_parent(cpd, topologies);
        topology t2 = select_topo_parent(cpd, topologies);
        
        LOG(8, 0, 0, "parents t1=%2.10f,t2=%2.10f ...\r\n", t1.fitness, t2.fitness);
        
        topology child;
        
        child.fitness = (t1.fitness + t2.fitness) / 2; // an estimated fitness, average of the parents
        
        // calculate an adjacency matrix for each parent's associated topology for use in
        // generating child topology ...
        
        std::vector<std::vector<int>> m1 = create_adjaceny_matrix(t1.comm, world_size);
        std::vector<std::vector<int>> m2 = create_adjaceny_matrix(t2.comm, world_size);
        
        std::vector<std::vector<int>> child_matrix;
        child_matrix.resize(world_size);
        
        int comm_count = 0;
        
        while(comm_count == 0) {
    
            for(int i=0; i<m1.size(); i++) {  // child matrix row
                
                child_matrix[i].resize(world_size);
                
                for(int j=0; j<m2.size(); j++) { // child matrix column
                    
                    if(i != j) {  // we don't want an island sending migrants to itself
                        
                        if(rand()%2 == 1) { // coin flip, heads take the row index value from parent 1
                            LOG(10, 0, 0, "assigning child m1[%d][%d] -> %d\r\n", i, j, m1[i][j]);
                            child_matrix[i][j] = m1[i][j];
                        } else { // tails, take it from parent 2
                            LOG(10, 0, 0, "assigning child m2[%d][%d] -> %d\r\n", i, j, m2[i][j]);
                            child_matrix[i][j] = m2[i][j];
                        }
                        
                        if(child_matrix[i][j] == 1) { comm_count++; }
                        
                    } else {
                        
                        child_matrix[i][j] = 0;
                        
                    }
                    
                }
                
                LOG(10, 0, 0, "------\r\n");
                
            }
            
        }
        
        // mutation, interate through the matrix and flip the bit with probability m,
        // also considering any sparsity constraints ...
        
        for(int i=0; i<child_matrix.size(); i++) {

            for(int j=0; j<child_matrix[i].size(); j++) {

                if(rand()/(RAND_MAX+1.0) < config::mutation_rate) {
                
                    if(child_matrix[i][j] == 0 && rand()/(RAND_MAX+1.0) < config::sparsity) {
                        
                        child_matrix[i][j] = 1;
                        comm_count++;
                        
                    }
                    
                    if(child_matrix[i][j] == 1 && rand()/(RAND_MAX+1.0) > config::sparsity && comm_count > 1) {
                        
                        child_matrix[i][j] = 0;
                        comm_count--;
                        
                    }
                    
                }

            }

        }
        
        // convert the adjaceny matrix to a sender and receiver arrays for use with MPI send\recv ...
        
        child.comm = create_group(child_matrix);

        LOG(8, 0, 0, "child %d created from matrix with %lu senders and %lu receivers\r\n", n, child.comm[n].senders.size(), child.comm[n].receivers.size());

        children.push_back(child);

    }
    
    return children;
    
}

// used to create the initial topology population ...

std::vector<group> create_dyn_topology(std::vector<int> ids) {
    
    std::vector<group> topology;
    
    std::vector<std::vector<int>> matrix = create_dyn_adjaceny_matrix((int)ids.size());
    topology = create_group(matrix);
    
    return topology;
    
}


#endif /* island_h */

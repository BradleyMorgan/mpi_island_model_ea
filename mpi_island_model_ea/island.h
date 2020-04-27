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

void create_topology(island &isle, int world_size) {

    int next = isle.id+1 < world_size ? isle.id+1 : 0;
    int prev = isle.id-1 < 0 ? (int)world_size-1 : isle.id-1;
    
    isle.receivers.push_back(next);
    isle.senders.push_back(prev);
    
    LOG(8, 0, 0, "%d -> %d -> %d\r\n", prev, isle.id, next);
    
}

struct group {
    
    int node;
    
    double fitness = 0.0;
    
    std::vector<int> senders;
    std::vector<int> receivers;
    
};

struct topology {
  
    int rounds = 0;
    
    double fitness = 0.0;
    double round_fitness = 0.0;
    double selection_distribution;
    
    std::vector<group> comm;
    
};

std::vector<group> create_group(std::vector<std::vector<int>> &matrix) {
    
    LOG(10, 0, 0, "creating group...\r\n");
    
    std::vector<group> comm;
    
    comm.resize(matrix.size());
    
    for(int i=0; i<matrix.size(); i++) {
        
        //std::vector<int> *receivers = new std::vector<int>();
        
        for(int j=0; j<matrix[i].size(); j++) {
        
            LOG(10, 0, 0, "%d ", matrix[i][j]);
            
            if(matrix[i][j] == 1) {
                
                LOG(10, 0, 0, "adding receiver to %d and sender to %d\r\n", i, j);
                //receivers->push_back(j);
                
                comm[i].receivers.push_back(j);
                comm[j].senders.push_back(i);
                
            }
            
        }
        
        LOG(10, 0, 0, "\r\n");
        
    }
    
    return comm;
    
}

std::vector<std::vector<int>> create_adjaceny_matrix(std::vector<group> &comm, int world_size) {
    
    std::vector<std::vector<int>> matrix;
    matrix.resize(world_size);

    for(int i=0; i<world_size; i++) {

        matrix[i].resize(world_size);
        
        for(int j=0; j<world_size; j++) {
            
            if(std::find(comm[i].receivers.begin(), comm[i].receivers.end(), j) != comm[i].receivers.end()) {
                matrix[i][j] = 1;
            } else {
                matrix[i][j] = 0;
            }
            
        }
        
    }
    
    return matrix;
    
}

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
        
        child.fitness = (t1.fitness + t2.fitness) / 2;
        std::vector<std::vector<int>> m1 = create_adjaceny_matrix(t1.comm, world_size);
        std::vector<std::vector<int>> m2 = create_adjaceny_matrix(t2.comm, world_size);
        
        std::vector<std::vector<int>> child_matrix;
        child_matrix.resize(world_size);
        
        for(int i=0; i<m1.size(); i++) {
            
            child_matrix[i].resize(world_size);
            
            for(int j=0; j<m2.size(); j++) {
                
                if(i != j) {
                    
                    if(rand()%2 == 1) {
                        LOG(10, 0, 0, "assigning child m1[%d][%d] -> %d\r\n", i, j, m1[i][j]);
                        child_matrix[i][j] = m1[i][j];
                    } else {
                        LOG(10, 0, 0, "assigning child m2[%d][%d] -> %d\r\n", i, j, m2[i][j]);
                        child_matrix[i][j] = m1[i][j];
                    }
                    
                } else {
                    
                    child_matrix[i][j] = 0;
                    
                }
                
            }
            
            LOG(10, 0, 0, "------\r\n");
            
        }
        
        if(rand()/(RAND_MAX+1.0) < config::mutation_rate) {

            for(int i=0; i<child_matrix.size(); i++) {

                for(int j=0; j<child_matrix[i].size(); j++) {

                    if(rand()/(RAND_MAX+1.0) < config::mutation_rate) {

                        child_matrix[i][j] == 0 ? child_matrix[i][j] = 1 : child_matrix[i][j] = 0;

                    }

                }

            }

        }
        
        child.comm = create_group(child_matrix);

        LOG(8, 0, 0, "child %d created from matrix with %lu senders and %lu receivers\r\n", n, child.comm[n].senders.size(), child.comm[n].receivers.size());
        children.push_back(child);

    }
    
    return children;
    
}

void add_neighbors(int node, int size, int world_size, std::vector<group> &topology) {
    
    LOG(6, 0, 0, "entered add_neighbors with size %d for node %d...\r\n", size, node);
    
    if(size > 0) {
        
        std::vector<group>::iterator it;
        
        topology[node].node = node;
        
        bool rfound = false;
        
        for(it=topology.begin(); it!=topology.end(); ++it) {
            // search for this node in the set of senders for the currently iterated node ...
            std::vector<int>::iterator s = std::find(it->senders.begin(), it->senders.end(), node);
            // if this node is found as a sender, and it does not already exist, add it to the receivers list ...
            if(s != it->senders.end() && std::find(topology[node].receivers.begin(), topology[node].receivers.end(), it->node) == topology[node].receivers.end()) {
                LOG(6, 0, 0, "assigning found receiver %d -> %d\r\n", node, it->node);
                topology[node].receivers.push_back(it->node);
                rfound = true;
            }
        }
        
        bool sfound = false;
        
        for(it=topology.begin(); it!=topology.end(); ++it) {
            std::vector<int>::iterator s = std::find(it->receivers.begin(), it->receivers.end(), node);
            if(s != it->receivers.end() && std::find(topology[node].senders.begin(), topology[node].senders.end(), it->node) == topology[node].senders.end()) {
                LOG(6, 0, 0, "assigning found sender %d -> %d\r\n", it->node, node);
                topology[node].senders.push_back(it->node);
                sfound = true;
            }
        }
        
        
        if(!rfound) {
    
            if(topology[node].senders.size() == world_size-1) { return; }
            
            int rnd_source = node;
            
            while(rnd_source == node || std::find(topology[node].senders.begin(), topology[node].senders.end(), rnd_source) != topology[node].senders.end()) {
                rnd_source = rand()%topology.size();
            }
            
            LOG(6, 0, 0, "assigning random sender %d -> %d\r\n", rnd_source, node);
            topology[node].senders.push_back(rnd_source);
            
            add_neighbors(rnd_source, size-1, world_size, topology);
                
        }
        
        if(!sfound) {
            
            if(topology[node].receivers.size() == world_size-1) { return; }
            
            int rnd_target = node;
            
            while(rnd_target == node || std::find(topology[node].receivers.begin(), topology[node].receivers.end(), rnd_target) != topology[node].receivers.end()) {
                rnd_target = rand()%topology.size();
            }
            
            LOG(6, 0, 0, "assigning random receiver %d -> %d\r\n", node, rnd_target);
            topology[node].receivers.push_back(rnd_target);
            
            add_neighbors(rnd_target, size-1, world_size, topology);
        
        }
        
    } else {
        
        std::vector<group>::iterator it;
        
        for(it=topology.begin(); it!=topology.end(); ++it) {
            std::vector<int>::iterator s = std::find(it->senders.begin(), it->senders.end(), node);
            if(s != it->senders.end() && std::find(topology[node].receivers.begin(), topology[node].receivers.end(), it->node) == topology[node].receivers.end()) {
                LOG(6, 0, 0, "assigning found receiver %d -> %d\r\n", node, it->node);
                topology[node].receivers.push_back(it->node);
            }
        }
        
        for(it=topology.begin(); it!=topology.end(); ++it) {
            std::vector<int>::iterator s = std::find(it->receivers.begin(), it->receivers.end(), node);
            if(s != it->receivers.end() && std::find(topology[node].senders.begin(), topology[node].senders.end(), it->node) == topology[node].senders.end()) {
                LOG(6, 0, 0, "assigning found sender %d -> %d\r\n", it->node, node);
                topology[node].senders.push_back(it->node);
            }
        }
        
        LOG(6, 0, 0, "returning from add_neighbors ...\r\n");
        
        return;
        
    }
    
}

std::vector<group> create_dyn_topology(std::vector<int> ids) {
    
    std::vector<group> topology;
    topology.resize(ids.size());
    
    add_neighbors(0, (int)ids.size(), (int)ids.size(), topology);
    
    return topology;
    
}


#endif /* island_h */

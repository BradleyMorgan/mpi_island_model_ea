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
            LOG(6, "island %d waiting for migrant from island %d ... \r\n", this->id, this->senders[i]);
            MPI_Recv(&x, DIM, MPI_DOUBLE, this->senders[i], 0, MPI_COMM_WORLD, &migrant_status);
            this->population[rand()%population.size()].input = x;
            LOG(6, "island %d received migrant from island %d: [%f,%f] with status %d\r\n", this->id, migrant_status.MPI_SOURCE, this->population[0].input[0], this->population[0].input[0], migrant_status.MPI_ERROR);
        }
        
    }
    
    void send_migrant(MPI_Comm &comm) {
        
        for(int i=0; i<this->receivers.size(); i++) {
            LOG(6, "island %d sending migrant to island %d ... \r\n", this->id, this->receivers[i]);
            MPI_Send(&this->population[i].input, DIM, MPI_DOUBLE, this->receivers[i], 0, comm);
            LOG(6, "island %d sent migrant to island %d: [%f,%f]\r\n", this->id, this->receivers[i], this->population[i].input[0], this->population[i].input[1]);
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
    
    LOG(8, "%d -> %d -> %d\r\n", prev, isle.id, next);
    
}

struct group {
    
    int node;
    
    double fitness;
    
    std::vector<int> senders;
    std::vector<int> receivers;
    
};

void add_neighbors(int node, int size, std::vector<group> &topology) {
    
    if(size > 0) {
        
        std::vector<group>::iterator it;
        
        topology[node].node = node;
        
        bool rfound = false;
        
        for(it=topology.begin(); it!=topology.end(); ++it) {
            std::vector<int>::iterator s = std::find(it->senders.begin(), it->senders.end(), node);
            if(s != it->senders.end() && std::find(topology[node].receivers.begin(), topology[node].receivers.end(), it->node) == topology[node].receivers.end()) {
                LOG(6, "assigning found receiver %d -> %d\r\n", node, it->node);
                topology[node].receivers.push_back(it->node);
                rfound = true;
            }
        }
        
        bool sfound = false;
        
        for(it=topology.begin(); it!=topology.end(); ++it) {
            std::vector<int>::iterator s = std::find(it->receivers.begin(), it->receivers.end(), node);
            if(s != it->receivers.end() && std::find(topology[node].senders.begin(), topology[node].senders.end(), it->node) == topology[node].senders.end()) {
                LOG(6, "assigning found sender %d -> %d\r\n", it->node, node);
                topology[node].senders.push_back(it->node);
                sfound = true;
            }
        }
        
        
        if(!rfound) {
    
            int rnd_source = node;
            
            while(rnd_source == node || std::find(topology[node].senders.begin(), topology[node].senders.end(), rnd_source) != topology[node].senders.end()) {
                rnd_source = rand()%topology.size();
            }
            
            LOG(6, "assigning random sender %d -> %d\r\n", rnd_source, node);
            topology[node].senders.push_back(rnd_source);
            
            add_neighbors(rnd_source, size-1, topology);
                
        }
        
        if(!sfound) {
            
            int rnd_target = node;
            
            while(rnd_target == node || std::find(topology[node].receivers.begin(), topology[node].receivers.end(), rnd_target) != topology[node].receivers.end()) {
                rnd_target = rand()%topology.size();
            }
            
            LOG(6, "assigning random receiver %d -> %d\r\n", node, rnd_target);
            topology[node].receivers.push_back(rnd_target);
            
            add_neighbors(rnd_target, size-1, topology);
        
        }
        
    } else {
        
        std::vector<group>::iterator it;
        
        for(it=topology.begin(); it!=topology.end(); ++it) {
            std::vector<int>::iterator s = std::find(it->senders.begin(), it->senders.end(), node);
            if(s != it->senders.end() && std::find(topology[node].receivers.begin(), topology[node].receivers.end(), it->node) == topology[node].receivers.end()) {
                LOG(6, "assigning found receiver %d -> %d\r\n", node, it->node);
                topology[node].receivers.push_back(it->node);
            }
        }
        
        for(it=topology.begin(); it!=topology.end(); ++it) {
            std::vector<int>::iterator s = std::find(it->receivers.begin(), it->receivers.end(), node);
            if(s != it->receivers.end() && std::find(topology[node].senders.begin(), topology[node].senders.end(), it->node) == topology[node].senders.end()) {
                LOG(6, "assigning found sender %d -> %d\r\n", it->node, node);
                topology[node].senders.push_back(it->node);
            }
        }
        
        return;
        
    }
    
}

std::vector<group> create_dyn_topology(std::vector<int> *ids) {
    
    std::vector<group> topology;
    topology.resize(ids->size());
    
    add_neighbors(0, (int)ids->size(), topology);
    
    return topology;
    
}

std::vector<group> create_dynamic_topology(std::vector<int> *ids) {
    
    std::vector<group> topology;
    
    for(int i=0; i<ids->size(); i++) {
        
        group g;
        
        g.node = i;
        
        LOG(6, "assigning group %d\r\n", g.node);
        
        std::vector<group>::iterator it;
        
        bool sfound = false;
        bool rfound = false;
        
        for(it=topology.begin(); it!=topology.end(); ++it) {
            std::vector<int>::iterator s = std::find(it->senders.begin(), it->senders.end(), g.node);
            if(s != it->senders.end()) {
                LOG(6, "assigning found receiver %d\r\n", it->node);
                g.receivers.push_back(it->node);
                rfound = true;
            }
            
            std::vector<int>::iterator r = std::find(it->receivers.begin(), it->receivers.end(), g.node);
            if(r != it->receivers.end()) {
                LOG(6, "assigning found sender %d\r\n", it->node);
                g.senders.push_back(it->node);
                sfound = true;
            }
        }
        
        if(!sfound) {
        
            int rnd_source = i;
            
            while(rnd_source == i) {
                rnd_source = rand()%ids->size();
            }
            
            LOG(6, "assigning random sender %d\r\n", rnd_source);
            g.senders.push_back(rnd_source);
            
        }
        
        if(!rfound) {
        
            int rnd_target = i;
            
            while(rnd_target == i) {
                rnd_target = rand()%ids->size();
            }
            
            LOG(6, "assigning random receiver %d\r\n", rnd_target);
            g.receivers.push_back(rnd_target);
            
        }
        
        LOG(10, "adding node %d to topology\r\n", g.node);
        
        topology.push_back(g);
    
        
    }
    
    for(int i=0; i<ids->size(); i++) {
    
        for(int k=0; k<topology[i].senders.size(); k++) {
            
            printf("%d -> ", topology[i].senders[k]);
            
        }
        
        printf("[%d] -> ", topology[i].node);
        
        for(int k=0; k<topology[i].receivers.size(); k++) {
            
            printf("%d -> ", topology[i].receivers[k]);
            
        }
        
        
        printf("\r\n");
        
    }
    
    return topology;
    
}

#endif /* island_h */

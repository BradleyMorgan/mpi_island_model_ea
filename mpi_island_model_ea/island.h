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
    
    void receive_migrant() {
        
        std::array<double, DIM> x;
        
        MPI_Status migrant_status;
        
        MPI_Recv(&x, DIM, MPI_DOUBLE, this->senders[0], 0, MPI_COMM_WORLD, &migrant_status);
                
        this->population[rand()%population.size()].input = x;
                
        printf("island %d received migrant from island %d: [%f,%f] with status %d\r\n", this->id, migrant_status.MPI_SOURCE, this->population[0].input[0], this->population[0].input[0], migrant_status.MPI_ERROR);
        
    }
    
    void send_migrant() {
                
        MPI_Send(&this->population[0].input, DIM, MPI_DOUBLE, this->receivers[0], 0, MPI_COMM_WORLD);
                
        printf("island %d sending migrant to island %d: [%f,%f]\r\n", this->id, this->receivers[0], this->population[0].input[0], this->population[0].input[1]);
                
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
    
    //printf("%d -> %d -> %d\r\n", prev, isle.id, next);
    
}

struct group {
    
    int node;
    
    std::vector<int> senders;
    std::vector<int> receivers;
    
};

std::vector<group> create_dynamic_topology(std::vector<int> *ids) {
    
    std::vector<group> topology;
    
    for(int i=0; i<ids->size(); i++) {
        
        group g;
        
        g.node = i;
        
        printf("assigning group %d\r\n", g.node);
        
        std::vector<group>::iterator it;
        
        bool sfound = false;
        bool rfound = false;
        
        for(it=topology.begin(); it!=topology.end(); ++it) {
            std::vector<int>::iterator s = std::find(it->senders.begin(), it->senders.end(), g.node);
            if(s != it->senders.end()) {
                printf("assigning found receiver %d\r\n", it->node);
                g.receivers.push_back(it->node);
                rfound = true;
            }
            
            std::vector<int>::iterator r = std::find(it->receivers.begin(), it->receivers.end(), g.node);
            if(r != it->receivers.end()) {
                printf("assigning found sender %d\r\n", it->node);
                g.senders.push_back(it->node);
                sfound = true;
            }
        }
        
        if(!sfound) {
        
            int rnd_source = i;
            
            while(rnd_source == i) {
                rnd_source = rand()%ids->size();
            }
            
            printf("assigning random sender %d\r\n", rnd_source);
            g.senders.push_back(rnd_source);
            
        }
        
        if(!rfound) {
        
            int rnd_target = i;
            
            while(rnd_target == i) {
                rnd_target = rand()%ids->size();
            }
            
            printf("assigning random receiver %d\r\n", rnd_target);
            g.receivers.push_back(rnd_target);
            
        }
        
        printf( "adding node %d to topology\r\n", g.node);
        
        topology.push_back(g);
        
//        printf("erasing %d\r\n", g.node);
//
//        ids->erase(std::find(ids->begin(), ids->end(), g.node));
            
        printf("%d -> %d -> %d\r\n", g.senders[0], g.node, g.receivers[0]);
        
    }

    return topology;
    
}

#endif /* island_h */

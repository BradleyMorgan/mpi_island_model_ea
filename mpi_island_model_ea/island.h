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
    
    int send[1];
    int receive[1];
    
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
                
        //printf("island %d received migrant from island %d: [%f,%f] with status %d\r\n", this->id, migrant_status.MPI_SOURCE, this->population[0].input[0], this->population[0].input[0], migrant_status.MPI_ERROR);
        
    }
    
    void send_migrant() {
                
        MPI_Send(&this->population[0].input, DIM, MPI_DOUBLE, this->receivers[0], 0, MPI_COMM_WORLD);
                
        //printf("island %d sending migrant to island %d: [%f,%f]\r\n", this->id, this->receivers[0], this->population[0].input[0], this->population[0].input[1]);
                
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
    
    isle.receive[0] = next;
    isle.send[0] = prev;
    
    //printf("%d -> %d -> %d\r\n", prev, isle.id, next);
    
}

#endif /* island_h */

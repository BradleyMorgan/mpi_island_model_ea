//
//  island.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 2/26/20.
//  Copyright © 2020 Bradley Morgan. All rights reserved.
//

#ifndef island_h
#define island_h

#pragma mark DATATYPE: island{}

// representation of an island model EA entity that performs evaluations on
// a subpopulation of the primary (rastrigin) objective function solutions.
// islands perform migrations of individuals using designated communication
// channels as implemented in @island{receive_migrant}, @island{send_migrant}

struct island {
    
    int id;
    
    double total_fitness;
    double average_fitness;
    
    std::vector<solution> population;
    std::vector<double> cpd;
    std::vector<int> senders;
    std::vector<int> receivers;
    
    MPI_Comm tcomm;
    
    struct calculate {
    
        static void total_fitness(island &p);
        static void average_fitness(island &p);
        static void cpd(island &p);
        
    };
    
    struct migration {
        
        static void send(island &p, MPI_Comm &comm);
        static void receive(island &p, MPI_Comm &comm);
        
    };
    
    void init() {
        
        LOG(6, 0, 0, "initializing island %d\r\n", this->id);
            
        this->total_fitness = 0.0;
        this->average_fitness = 0.0;
    
        this->tcomm = MPI_COMM_WORLD;
        
        this->senders.clear();
        this->receivers.clear();
        this->population.clear();
        this->cpd.clear();

        LOG(6, 0, 0, "island %d initalized\r\n", this->id);
        
    }

};

#pragma mark FUNCTION: island::calculate::total_fitness()

// calculate the island's total fitness for distribution ...

void island::calculate::total_fitness(island &p) {
    
    p.total_fitness = 0.0;
    
    std::vector<solution>::iterator it;
    
    for(it = p.population.begin(); it != p.population.end(); ++it) {
        p.total_fitness += it->fitness;
    }
    
    LOG(6, 0, 0, "island::calculate::total_fitness island %d = %f\r\n", p.id, p.total_fitness);
    
}

#pragma mark FUNCTION: island::calculate::average_fitness()

void island::calculate::average_fitness(island &p) {
    
    island::calculate::total_fitness(p);
    
    p.average_fitness = p.total_fitness / p.population.size();
    
    LOG(6, 0, 0, "island::calculate::average_fitness island %d = %f\r\n", p.id, p.average_fitness);
    
}

#pragma mark FUNCTION: island::calculate::cpd()

// this function calculates the cumulative probability distribution to be used by
// the fitness proportional (roulette wheel) selection ...

void island::calculate::cpd(island &p) {
    
    LOG(6, 0, 0, "island %d, population size %lu calculating cpd ...\r\n", p.id, p.population.size());
    
    double cumulative_probability = 0.0;
    
    LOG(6, 0, 0, "calculating total fitness ...\r\n");
    
    island::calculate::total_fitness(p);
    
    // std::sort(p.population.begin(), p.population.end(), compare_fitness);
    // std::reverse(p.population.begin(), p.population.end());
    
    LOG(6, 0, 0, "island %d total fitness = %f\r\n", p.id, p.total_fitness);
    
    p.cpd.clear();
    
    LOG(7, 0, 0, "calculating island %d (population size = %lu) selection distribution\r\n", p.id, p.population.size());
    
    for(int i=0; i<p.population.size(); i++) {

        p.population[i].selection_distribution = (double)p.population[i].fitness / p.total_fitness;

        LOG(8, 0, 0, "calculating island %d solution %d fitness %f selection distribution = %f\r\n", p.id, i, p.population[i].fitness, p.population[i].selection_distribution);
        
        cumulative_probability += p.population[i].selection_distribution;
        
        LOG(8, 0, 0, "island %d solution %d cumulative prob = %f\r\n", p.id, i, cumulative_probability);
        
        p.cpd.push_back(cumulative_probability);
        
    }
    
}


#pragma mark FUNCTION: island::migration::receive()

// initiate a receive operation for every island in this island's senders list ...

void island::migration::receive(island &p, MPI_Comm &comm) {
    
    LOG(5, 0, 0, "island::migration::receive(%d) operations for %lu senders\r\n", p.id, p.senders.size());
    
    std::array<double, DIM> x;
    
    MPI_Status migrant_status;
    
    for(int i=0; i<p.senders.size(); i++) {
        LOG(5, 0, 0, "island %d waiting for migrant from island %d ... \r\n", p.id, p.senders[i]);
        MPI_Recv(&x, DIM, MPI_DOUBLE, p.senders[i], 0, comm, &migrant_status);
        p.population[rand()%p.population.size()].input = x;
        LOG(5, 0, 0, "island %d received migrant from island %d: [%f,%f] with status %d\r\n", p.id, migrant_status.MPI_SOURCE, p.population[0].fitness, p.population[0].fitness, migrant_status.MPI_ERROR);
    }
    
}

#pragma mark FUNCTION: island::migration::send()

// initiate a send operation for every island in this island's receivers list ...

void island::migration::send(island &p, MPI_Comm &comm) {
    
    LOG(5, 0, 0, "island::migration::send(%d) initiating operations for %lu receivers\r\n", p.id, p.receivers.size());
    
    for(int i=0; i<p.receivers.size(); i++) {
        LOG(6, 0, 0, "island %d sending migrant to island %d ... \r\n", p.id, p.receivers[i]);
        MPI_Send(&p.population[i].input, DIM, MPI_DOUBLE, p.receivers[i], 0, comm);
        LOG(6, 0, 0, "island %d sent migrant to island %d: [%f,%f]\r\n", p.id, p.receivers[i], p.population[i].input[0], p.population[i].input[1]);
    }
            
}


#endif /* island_h */

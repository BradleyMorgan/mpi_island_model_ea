//
//  dtype_island.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 7/22/21.
//  Copyright © 2021 Bradley Morgan. All rights reserved.
//

#ifndef dtype_island_h
#define dtype_island_h

struct visa {

    int eval;
    int source;
    int destination;
    
    char genome_id[64];
    
    visa() : eval(0), source(0), destination(0) {}
    visa(int e, int s, int d, char g[64]) : eval(e), source(s), destination(d) { strcpy(genome_id, g); }
    
};

#pragma mark DATATYPE: island{}

// representation of an island model EA entity that performs evaluations on
// a subpopulation of the primary (rastrigin) objective function solutions.
// islands perform migrations of individuals using designated communication
// channels as implemented in @island{receive_migrant}, @island{send_migrant}

struct island {
    
    int id = 0;
    
    double total_fitness = 0.0;
    double average_fitness = 0.0;
    
    std::vector<genome> population;
    
    std::vector<double> cpd = {};
    std::vector<int> senders = {};
    std::vector<int> receivers = {};
    std::vector<visa> visas;
    
    MPI_Comm tcomm;
    
    struct calculate {
    
        static void total_fitness(island &p);
        static void average_fitness(island &p);
        static void cpd(island &p);
        
    };
    
    struct migration {
        
        static void send(island &p, MPI_Datatype &d, int &eval);
        static void receive(island &p, MPI_Datatype &d, int &eval);
        
    };
    
    void init() {
        
        LOG(6, 0, 0, "initializing island %d\r\n", this->id);
        
        this->population.clear();
        this->population.resize(config::island_mu);
        this->total_fitness = 0.0;
        this->average_fitness = 0.0;
        this->tcomm = MPI_COMM_WORLD;
        this->senders.clear();
        this->receivers.clear();
        this->cpd.clear();

        LOG(6, 0, 0, "island %d initalized\r\n", this->id);
        
    }

};


#endif /* dtype_island_h */

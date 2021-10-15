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
// a subpopulation of the primary (rastrigin) objective functionsolver.solutions.population.
// islands perform migrations of individuals using designated communication
// channels as implemented in @island{receive_migrant}, @island{send_migrant}

struct island {
    
    int id = 0;
    
    double total_fitness = 0.0;
    double average_fitness = 0.0;
    
    std::vector<solution> population;
    
    std::vector<double> cpd = {};
    std::vector<int> senders = {};
    std::vector<int> receivers = {};
    std::vector<visa> visas;
    
    MPI_Comm tcomm;
    
    struct calculate {
    
        void total_fitness(island &p);
        void average_fitness(island &p);
        void cpd(island &p);
        
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
    
    calculate calculator;
    
};

// comparator for parent fitness values ...

template<typename genome> bool compare_fitness(const genome &p1, const genome &p2) {
    return p1.fitness < p2.fitness;
}


void island::calculate::cpd(island &p) {
    
    LOG(10, 0, 0, "generating cpd\r\n");

    double cumulative_probability = 0.0;

    LOG(6, 0, 0, "calculating total fitness ...\r\n");

    LOG(6, 0, 0, "sorting population descending fitness ...\r\n");

    std::sort(p.population.begin(), p.population.end(), compare_fitness<solution>);
    std::reverse(p.population.begin(), p.population.end());

    LOG(6, 0, 0, "objective %d cpd calculation total fitness = %f", p.id, p.total_fitness);

    for(int i=0; i < p.population.size(); i++) {

        p.population[i].selection_distribution = p.population[i].fitness / p.total_fitness;

        LOG(8, 0, 0, "calculated island %d solution %d fitness %f selection distribution = %f\r\n", p.id, i, p.population[i].fitness, p.population[i].selection_distribution);

        cumulative_probability += p.population[i].selection_distribution;

        LOG(8, 0, 0, "solution %d cumulative prob = %f\r\n", i, cumulative_probability);

        p.cpd.push_back(cumulative_probability);

    }
    
}


#endif /* dtype_island_h */

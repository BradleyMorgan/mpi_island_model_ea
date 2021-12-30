//
//  dtype_island.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 7/22/21.
//  Copyright © 2021 Bradley Morgan. All rights reserved.
//

#ifndef dtype_island_h
#define dtype_island_h

struct island_stats {
    
    int log_interval = 1;
    
    double local_t = 0.0;
    double min_run_t = 0.0;
    double max_run_t = 0.0;
    double avg_run_t = 0.0;
    
    int departures = 0;
    int arrivals = 0;
    int migrations = 0;
    
    char island_out[64];
    
    FILE *island_log;
    
};

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
    
    int id = mpi.id;
    
    std::vector<solution> population;
    
    std::vector<int> senders = {};
    std::vector<int> receivers = {};
    
    std::vector<visa> visas;
    
    island_stats stats;
    
    MPI_Comm tcomm;
    
    struct metric {
        
        struct calculate {

            void total_fitness();
            void average_fitness();
            void cpd();

        };
    
        struct values {
            
            std::vector<double> cpd = {};
            
            double total_fitness = 0.0;
            double average_fitness = 0.0;
            double fitness = 0.0;
            
            values() : cpd(0.0), fitness(0.0) {}
            
        };
        
        metric(void) {}
        
        values value;
        
    };
    
    struct migration {
        
        static void send(island &p, MPI_Datatype &d, int &eval);
        static void receive(island &p, MPI_Datatype &d, int &eval);
        
    };
    
    void init() {
        
        LOG(6, 0, 0, "ISLAND %d entered initialization\r\n", this->id);
        
        this->population.clear();
        this->population.resize(config::mu_sub);
        this->tcomm = MPI_COMM_WORLD;
        this->senders.clear();
        this->receivers.clear();

        sprintf(this->stats.island_out, "%s/island_%d.csv", config::stats_subpath_ea_1_o1, this->id);
        
        LOG(2, 0, 0, "ISLAND %d initalized with empty population [0,n] => [%f,%f] sized %lu log %s\r\n", this->id, this->population[0].fitness, this->population[config::mu_sub-1].fitness, this->population.size(), this->stats.island_out);
        
    }
    
    island() : id(mpi.id) {};
    
    metric metrics;
    
    void log();
    void cpd();
    void total_fitness();
    void average_fitness();
    
};

// comparator for parent fitness values ...

template<typename genome> bool compare_fitness(const genome &p1, const genome &p2) {
    return p1.fitness < p2.fitness;
}


void island::cpd() {
    
    LOG(10, 0, 0, "generating cpd\r\n");

    double cumulative_probability = 0.0;

    LOG(6, 0, 0, "calculating total fitness ...\r\n");

    LOG(6, 0, 0, "sorting population descending fitness ...\r\n");

    std::sort(this->population.begin(), this->population.end(), compare_fitness<solution>);
    std::reverse(this->population.begin(), this->population.end());

    LOG(6, 0, 0, "objective %d cpd calculation total fitness = %f", this->id, this->metrics.value.total_fitness);

    for(int i=0; i < this->population.size(); i++) {

        this->population[i].selection_distribution = this->population[i].fitness / this->metrics.value.total_fitness;

        LOG(8, 0, 0, "calculated island %d solution %d fitness %f selection distribution = %f\r\n", this->id, i, this->population[i].fitness, this->population[i].selection_distribution);

        cumulative_probability += this->population[i].selection_distribution;

        LOG(8, 0, 0, "solution %d cumulative prob = %f\r\n", i, cumulative_probability);

        this->metrics.value.cpd.push_back(cumulative_probability);

    }
    
}

void island::log() {

    this->stats.island_log = fopen(this->stats.island_out, "w");
    
    LOG(2,0,0, "\r\n\r\nIIIISSSSLAND %d,%f,%f,%lu,%d,%d,%d,%f\r\n\r\n", this->id, this->metrics.value.average_fitness, this->metrics.value.total_fitness, this->population.size(), this->stats.arrivals, this->stats.departures, this->stats.migrations, this->stats.local_t);
    
    fprintf(this->stats.island_log, "%d,%f,%f,%lu,%d,%d,%d,%f", this->id, this->metrics.value.average_fitness, this->metrics.value.total_fitness, this->population.size(), this->stats.arrivals, this->stats.departures, this->stats.migrations, this->stats.local_t);
    
    fclose(this->stats.island_log);
    
}

#endif /* dtype_island_h */

//
//  objective.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 8/18/21.
//  Copyright © 2021 Bradley Morgan. All rights reserved.
//

#ifndef objective_h
#define objective_h

#pragma mark DATATYPE: objective

template<typename genome> struct objective {
    
    int id = 0;
    int max_runs = 0;
    int max_evo_evals = 0;
    int max_fit_evals = 0;
    int run = 0;
    int eval = 0;
    int mu = 0;
    int lambda = 0;
    int islands = 0;
    int workers = 0;
    
    double mutation_rate = 0.0;
    
    std::vector<genome> population = {};
    
    template<typename f, typename e> void populate(e &multi, f function) { function(multi); }
    template<typename f, typename e> void distribute(e &multi, f function) { function(multi); }
    template<typename f, typename e> void evaluate(e &multi, genome &individual, f function) { function(multi, individual); }
    template<typename f> void calculate(f function) { function(*this); }
    
    template<typename f> genome select(f function) { return function(*this); }
    template<typename f> genome crossover(f function) { return function(*this); }
    
    template<typename f, typename e> void evolve(e &multi, f function) { function(multi); };
    template<typename f, typename e, typename d> void coevolve(e &multi, f function, d &dependent) {
        function(multi, dependent);
    };
    
    #pragma mark DATATYPE: metrics
    
    struct metrics {
        
        void cpd(objective<genome> &o);
        void fitness(objective<genome> &o);
      
        struct values {
            
            std::vector<double> cpd = {};
            double fitness = 0.0;
            
            values() : cpd(0.0), fitness(0.0) {}
            
        };
        
        metrics(void) {}
        
        values value;
        
    };
    
    metrics aggregate;
    
};

// comparator for parent fitness values ...

//template<typename genome> bool compare_fitness(const genome &p1, const genome &p2) {
//    return p1.fitness < p2.fitness;
//}


//template<typename genome> void fitness(objective<genome> &o) {
template<typename genome> void objective<genome>::metrics::fitness(objective<genome> &o) {

    o.aggregate.value.fitness = 0.0;
    
    for(typename std::vector<genome>::iterator it = o.population.begin(); it != o.population.end(); ++it) {
        LOG(10, 0, 0, "genome %d fitness = %f\r\n", o.id, it->fitness);
        o.aggregate.value.fitness += it->fitness;
    }

}

//template<typename genome> void cpd(objective<genome> &o) {
template<typename genome> void objective<genome>::metrics::cpd(objective<genome> &o) {

    LOG(10, 0, 0, "generating cpd\r\n");

    double cumulative_probability = 0.0;

    o.aggregate.value.cpd.clear();

    LOG(6, 0, 0, "calculating total fitness ...\r\n");

    //o.aggregate.cpd(o);

    LOG(6, 0, 0, "sorting population descending fitness ...\r\n");

    std::sort(o.population.begin(), o.population.end(), compare_fitness<genome>);
    std::reverse(o.population.begin(), o.population.end());

    LOG(6, 0, 0, "objective %d cpd calculation total fitness = %f", o.id, o.aggregate.value.fitness);

    for(int i=0; i < o.population.size(); i++) {

        o.population[i].selection_distribution = o.population[i].fitness / o.aggregate.value.fitness;

        LOG(8, 0, 0, "calculated island %d solution %d fitness %f selection distribution = %f\r\n", o.id, i, o.population[i].fitness, o.population[i].selection_distribution);

        cumulative_probability += o.population[i].selection_distribution;

        LOG(8, 0, 0, "solution %d cumulative prob = %f\r\n", i, cumulative_probability);

        o.aggregate.value.cpd.push_back(cumulative_probability);

    }

}


#endif /* objective_h */

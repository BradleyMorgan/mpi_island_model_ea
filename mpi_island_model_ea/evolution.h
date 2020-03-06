//
//  evolution.h
//  hello_mpi_topology
//
//  Created by Bradley Morgan on 2/24/20.
//  Copyright © 2020 Eight Bit Software Co. All rights reserved.
//

#ifndef evolution_h
#define evolution_h

#include <algorithm>
#include "island.h"
#include "config.h"
#include "utility.h"

// pick a parent from the population using the provided method ...

individual select_parent(const island &isle) {
    
    individual p;

    // implementation uses the single armed roulette wheel approach to select
    // μ parents from the population ...
    
    int i = 1;
    
    double r = ((double)rand()/(double)RAND_MAX);
    
    while (isle.cpd[i] < r ) { i++; }
    
    p = isle.population[i];
    
    return p;
        
}

// comparator for parent fitness values ...

bool compare_fitness(const individual &p1, const individual &p2) {
    
    return p1.fitness < p2.fitness;
    
}

void select_survivors(island &isle, std::vector<individual> &children, int island_mu) {

    // truncation: add new children to the population, and then kill the weakest
    
    isle.population.insert(isle.population.end(), children.begin(), children.end());
    
    std::sort(isle.population.begin(), isle.population.end(), compare_fitness);
    std::reverse(isle.population.begin(), isle.population.end());
    
    isle.population.erase(isle.population.begin()+island_mu, isle.population.end());
    
}

void mutate(individual &mutant) {

    mutant.input[0] = drand(-5.12, 5.12);
    
}

std::vector<individual> crossover(const island &isle, std::array<double, DIM> &offsets) {

    std::vector<individual> children;
    
    for(int i = 0; i < config::lambda; i++) {
        
        individual p1 = select_parent(isle);
        individual p2 = select_parent(isle);
        
        individual child;
        
        for(int j=0; j<DIM; j++) {
            if(rand()%2 == 1) {
                child.input[j] = p1.input[rand()%DIM];
            } else {
                child.input[j] = p2.input[rand()%DIM];
            }
        }
        
        if(rand()/(RAND_MAX+1.0) < config::mutation_rate) { mutate(child); }
        
        child.fitness = offset_rastrigin(child.input, offsets);
        
        children.push_back(child);
        
    }
    
    return children;
    
}

// select μ parents from the population using a probability based on the
// individual fitness relative to the entire population ...

std::vector<individual> fitness_proportional_selection(const island &isle) {
    
    // implementation uses the single armed roulette wheel approach to select
    // μ parents from the population ...
    
    std::vector<individual> parents;
   
    int curr_member = 0;
    
    while(curr_member < config::mu) {
        
        int i = 1;
        
        double r = ((double)rand()/(double)RAND_MAX);
        
        while(isle.cpd[i] < r ) { i++; }
        
        parents.push_back(isle.population[i]);
        
        curr_member++;
        
    }
    
    return parents;
    
}

#endif /* evolution_h */

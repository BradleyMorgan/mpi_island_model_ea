//
//  solution.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 4/12/21.
//  Copyright © 2021 Bradley Morgan. All rights reserved.
//

#ifndef solution_h
#define solution_h

#include "utility.h"

#pragma mark DATATYPE: solution{}

// representation of rastrigin solution

struct solution {
    
    std::array<double, DIM> input;
    
    double fitness;
    double selection_distribution;
    
};

struct offset {
    
    std::array<double, DIM> input;
    
    void generate(double min, double max, double step);
    
};

void offset::generate(double min, double max, double step) {
    
    std::array<double, DIM> offsets;
    std::vector<double> increments;
    
    for (double increment = min; increment <= max; increment += step) {
        increments.push_back(increment);
    }
    
    for(int i=0; i<DIM; i++) {
        offsets[i] = increments[rand()%increments.size()];
    }
    
    this->input = offsets;
    
}

std::array<double, DIM> generate_offsets(double min, double max, double step) {
    
    std::array<double, DIM> offsets;
    std::vector<double> increments;
    
    for (double increment = min; increment <= max; increment += step) {
        increments.push_back(increment);
    }
    
    for(int i=0; i<DIM; i++) {
        offsets[i] = increments[rand()%increments.size()];
    }
    
    return offsets;
    
}


// comparator for parent fitness values ...

bool compare_fitness(const solution &p1, const solution &p2) {
    
    return p1.fitness < p2.fitness;
    
}


#endif /* solution_h */

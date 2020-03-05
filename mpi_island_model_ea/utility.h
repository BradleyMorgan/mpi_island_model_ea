//
//  utility.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 2/26/20.
//  Copyright © 2020 Bradley Morgan. All rights reserved.
//

#ifndef utility_h
#define utility_h

#include <random>
#include "config.h"

// return a random double between min and max ...

double drand(double min, double max) {
    
    double f = (double)rand() / RAND_MAX;
    
    return min + f * (max - min);
    
}

// the rastrigin function is a non-linear, multimodal function
// with a large search space and a large number of local minima ...

double rastrigin(std::array<double, DIM> x) {
    
    double sum = 10 * DIM;
    
    for (unsigned int i = 0; i < DIM; i++) {
        sum += (std::pow(x[i],2) - (10 * std::cos(2 * M_PI * x[i])));
    }
    
    return sum;
}

double offset_rastrigin(std::array<double, DIM> x) {
    
    double sum = 10 * DIM;
    
    for (unsigned int i = 0; i < DIM; i++) {
        sum += (std::pow(x[i],2) - (10 * std::cos(2 * M_PI * x[i])));
    }
    
    return sum;
    
}

double random_uniform() {
    
    // create a normal distribution from which to select random integers ...
    
    unsigned seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count();
    
    std::default_random_engine generator(seed);
    
    std::uniform_real_distribution<double> d(-2.5, 2.5);
    
    double nd = d(generator);
    
    return nd;
    
}

std::array<double, DIM> generate_offsets(double min, double max, double step) {
    
    std::array<double, DIM> offsets;
    
    for (unsigned int i = 0; i < DIM; i++) {
        double offset = (min + ((double)rand() / RAND_MAX) * (max - min)) * step + min;
        offsets[i] = offset;
        random_uniform();
    }
    
    return offsets;
    
}

#endif /* utility_h */

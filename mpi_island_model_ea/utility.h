//
//  utility.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 2/26/20.
//  Copyright © 2020 Bradley Morgan. All rights reserved.
//

#ifndef utility_h
#define utility_h

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

#endif /* utility_h */

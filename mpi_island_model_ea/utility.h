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
#include <cpuid.h>
#include <array>

#pragma mark FUNCTION prob_true() --------

// utility function to find arbitrary event probability

bool prob_true(double p){
    return rand()/(RAND_MAX+1.0) < p;
}

#pragma mark FUNCTION: drand() --------

// return a random double between min and max ...

double drand(double min, double max) {
    
    double f = (double)rand() / RAND_MAX;
    
    return min + f * (max - min);
    
}

#pragma mark FUNCTION: prob_true() --------

// the rastrigin function is a non-linear, multimodal function
// with a large search space and a large number of local minima ...

double rastrigin(std::array<double, DIM> x) {
    
    double sum = 10 * DIM;
    
    for (unsigned int i = 0; i < DIM; i++) {
        sum += (std::pow(x[i],2) - (10 * std::cos(2 * M_PI * x[i])));
    }
    
    return sum * -1;
}

# pragma mark EA::OBJECTIVE::BENCHMARK ------------
# pragma mark FUNCTION: offset_rastrigin() --------

// accepts an array of floating point numbers @x[param:1] to be used as the dimensional
// input values for the offset rastrigin fitness calculation, and returns
// the calculated fitness using a corresponding array of offset values @offsets[param:2].
// the function applies the rastrigin calculation, iterating @dim[config.txt:12] for
// each dimension.

double offset_rastrigin(std::array<double, DIM> x, std::array<double, DIM> &offsets) {

    // rastrigin function: f(x) = A*n + [Σ(1,n,𝜆𝑖↦(x[i]^2-Acos(2πx[i])]
    
    // rastrigin function constants:
    // A=10
    // x[i] ∈ [-5.12,5.12]
    
    double sum = 10 * DIM;
    
    for (unsigned int i = 0; i < DIM; i++) {
    
        // add the corresponding offset to the gene value, which should make this
        // a little more difficult to solve ...
        
        double offset_gene = x[i] + offsets[i];
        
        // apply the rastrigin function to the gene and add the value to the
        // fitness ...
        
        sum += (std::pow(offset_gene,2) - (10 * std::cos(2 * M_PI * offset_gene)));
    
    }
    
    // rastrigin is a minimization problem, so fitness values closer to zero are
    // considered better. we negate the fitness so that we can perform more
    // logical sorting and comparison using better_fitness > not_as_good_fitness
    
    return sum * -1;
    
}

#pragma mark FUNCTION: uniqid()

// unique genome tags containing with a specifed fixed length prefix

char* uniqid(unsigned long long int instance) {
    
    unsigned long long int count = instance <= 0 ? 1 : instance;
    
    int prefix = config::id_field_prefix1.second <= 0 ? 3 : config::id_field_prefix1.second;
    int rank=config::id_field_prefix1.first+1;
    int fixed_len=prefix;
    int instance_len=ceil(log10(count));
    
    static char buffer[64];
    
    sprintf(buffer, "%d%0*d%0*llu", rank, fixed_len, 0, instance_len, count);
    
    return buffer;
}

#endif /* utility_h */

#pragma mark MACRO: LOG()

#define LOG(level, rank, target, format, ...) if(level <= stoi(config::items["loglevel"]) && rank==target){ fprintf(stderr, format, ##__VA_ARGS__);}


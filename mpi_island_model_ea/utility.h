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

#pragma mark FUNCTION prob_true()

// utility function to find arbitrary event probability

bool prob_true(double p){
    return rand()/(RAND_MAX+1.0) < p;
}


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
    
    return sum * -1;
}

# pragma mark FUNCTION: offset_rastrigin

// accepts an array of floating point numbers @x[param:1] to be used as the dimensional     |
// input values for the offset rastrigin fitness calculation, and returns                   |
// the calculated fitness using a corresponding array of offset values @offsets[param:2].   |
// the function applies the rastrigin calculation, iterating @dim[config.txt:12] for        |
// each dimension.                                                                          |

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

#define CPUID(INFO, LEAF, SUBLEAF) __cpuid_count(LEAF, SUBLEAF, INFO[0], INFO[1], INFO[2], INFO[3])

#define GETCPU(CPU) {                                   \
        uint32_t CPUInfo[4];                            \
        CPUID(CPUInfo, 1, 0);                           \
        /* CPUInfo[1] is EBX, bits 24-31 are APIC ID */ \
        if ( (CPUInfo[3] & (1 << 9)) == 0) {            \
          CPU = -1;  /* no APIC on chip */              \
        }                                               \
        else {                                          \
          CPU = (unsigned)CPUInfo[1] >> 24;             \
        }                                               \
        if (CPU < 0) CPU = 0;                           \
      }

#endif /* utility_h */

//#define LOG(level, format, ...) if(level <= stoi(config::items["loglevel"])){ fprintf(stderr, "%s -- %d -- ", __FUNCTION__, __LINE__); fprintf(stderr, format, ##__VA_ARGS__);}
#define LOG(level, rank, target, format, ...) if(level <= stoi(config::items["loglevel"]) && rank==target){ fprintf(stderr, format, ##__VA_ARGS__);}

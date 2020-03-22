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

double offset_rastrigin(std::array<double, DIM> x, std::array<double, DIM> &offsets) {
    
    double sum = 10 * DIM;
    
    for (unsigned int i = 0; i < DIM; i++) {
        double offset_gene = x[i] + offsets[i];
        sum += (std::pow(offset_gene,2) - (10 * std::cos(2 * M_PI * offset_gene)));
    }
    
    return sum * -1;
    
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

//void log(const char *message, int level) {
//    
//    if(level <= stoi(config::items["loglevel"])) {
//        
//        printf("%s", message);
//        
//    }
//    
//}

#define CPUID(INFO, LEAF, SUBLEAF) __cpuid_count(LEAF, SUBLEAF, INFO[0], INFO[1], INFO[2], INFO[3])

#define GETCPU(CPU) {                              \
        uint32_t CPUInfo[4];                           \
        CPUID(CPUInfo, 1, 0);                          \
        /* CPUInfo[1] is EBX, bits 24-31 are APIC ID */ \
        if ( (CPUInfo[3] & (1 << 9)) == 0) {           \
          CPU = -1;  /* no APIC on chip */             \
        }                                              \
        else {                                         \
          CPU = (unsigned)CPUInfo[1] >> 24;                    \
        }                                              \
        if (CPU < 0) CPU = 0;                          \
      }

#endif /* utility_h */

//#define LOG(level, format, ...) if(level <= stoi(config::items["loglevel"])){ fprintf(stderr, "%s -- %d -- ", __FUNCTION__, __LINE__); fprintf(stderr, format, ##__VA_ARGS__);}
#define LOG(level, format, ...) if(level <= stoi(config::items["loglevel"])){ fprintf(stderr, format, ##__VA_ARGS__);}

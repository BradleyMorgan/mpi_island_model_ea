//
//  ea.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 3/30/20.
//  Copyright © 2020 Bradley Morgan. All rights reserved.
//

#ifndef ea_h
#define ea_h

#include <sys/time.h>
#include "dtype_heirarchy.h"
#include "dtype_ea.h"
#include "dtype_solver.h"
#include "dtype_meta.h"
#include "island.h"
#include "topology.h"
#include "stats.h"

template<typename genome> genome parent(objective<genome> &o) {
    
    genome p;

    o.cpd();
    
    // implementation uses the single armed roulette wheel approach to select
    // an individual from the population
    
    int i = 1;
    
    // random double uniformly distributed between 0 and 1
    
    double r = ((double)rand()/(double)RAND_MAX);
    
    // spin the wheel
    
    while (o.aggregate.value.cpd[i] < r ) { i++; }
    
    p = o.population[i];
    
    return p;
        
}

template<typename genome> std::pair<genome, genome> binary_tournament(objective<genome> &o) {

    genome p1 = o.population[rand()%o.population.size()];
    genome p2 = o.population[rand()%o.population.size()];

    //std::sort(o.population.begin(), o.population.end(), compare_multi<topology>);
    
    std::pair<genome, genome> result = std::minmax(p1, p2, compare_multi<topology>);
    
    return result;
        
}

//template<typename variant> void ea_end(ea<variant> &solver, objective<solution> &obj) {
//    
//    LOG(5, 0, 0, "\r\n--- (%d) END EA ---\r\n", mpi.id);
//    
//}
//
//template<typename variant> void ea_end(ea<variant> &meta, objective<topology> &obj) {
//
//    LOG(5, 0, 0, "\r\n--- (%d) END EA ---\r\n", mpi.id);
//
//    if(mpi.id == 0) {
//
//        char canary[128];
//
//        sprintf(canary, "%s/end.txt", config::logs_subpath);
//        FILE *eaend = fopen(canary, "w");
//
//        fprintf(eaend, "ended at %lu", time(0));
//
//        fclose(eaend);
//
//    }
//
//    MPI_Finalize();
//   
//}


#endif /* ea_h */

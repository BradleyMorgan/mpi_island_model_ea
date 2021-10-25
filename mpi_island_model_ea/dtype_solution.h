//
//  dtype_solution.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 7/22/21.
//  Copyright © 2021 Bradley Morgan. All rights reserved.
//

#ifndef dtype_solution_h
#define dtype_solution_h

#pragma mark DATATYPE: @solution{}

struct solution {
    
    char id[64];
    
    std::array<double, DIM> input = {};
    
    double fitness = 0.0;
    double selection_distribution = 0.0;
    double group = 0.0;

    int source = 0;
    int locale = 0;
    int migrations = 0;
    int selected = 0;
    int survival = 0;
    
    std::array<char[64], 2> parents;
    
    solution() { strcpy(id, uniqid(sinstances++)); }

};

#endif /* dtype_solution_h */

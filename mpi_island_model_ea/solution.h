//
//  solution.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 4/12/21.
//  Copyright © 2021 Bradley Morgan. All rights reserved.
//

#ifndef solution_h
#define solution_h

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
    
    LOG(6, 0, 0, "generating offsets...\r\n");
    
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

/*----------------------------------------------------------------------------*/

#pragma mark EA::FUNCTION::TEMPLATES:INITIALIZATION: <solution>

// specialization of generic @objective{} methods for performing pre and post
// interval initialization and cleanup

#endif /* solution_h */

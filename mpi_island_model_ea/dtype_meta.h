//
//  dtype_meta.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 10/17/21.
//  Copyright © 2021 Bradley Morgan. All rights reserved.
//

#ifndef dtype_meta_h
#define dtype_meta_h

struct meta : ea {
    
    objective<topology> topologies;
    
    meta(solver &solver) {
        
        this->topologies.mu = config::ea_2_mu;
        this->topologies.lambda = config::ea_2_lambda;
        this->topologies.max_runs = config::ea_2_runs;
        this->topologies.mutation_rate = config::ea_2_mutation_rate;
        this->topologies.max_evo_cycles = config::ea_2_max_evo_cycles;
        this->topologies.max_fit_evals = config::ea_2_max_fit_evals;
        this->variant = solver.variant;
        this->init_duration += (MPI_Wtime() - this->start);
        
    }
    
    //void end(objective<topology> &obj) { obj.end(*this); }
    
};

#endif /* dtype_meta_h */

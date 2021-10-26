//
//  dtype_solver.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 10/17/21.
//  Copyright © 2021 Bradley Morgan. All rights reserved.
//

#ifndef dtype_solver_h
#define dtype_solver_h

struct solver: ea {
    
    std::array<double, DIM> offsets;
    
    objective<solution> solutions;
    
    solver() {
    
        this->variant.init(*this);
    
        LOG(2, this->variant.isle.id, 0, "SOLVER EA INIT\r\n");
        
        this->solutions.mu = config::ea_1_mu;
        this->solutions.lambda = config::ea_1_lambda;
        this->solutions.max_runs = config::ea_1_runs;
        this->solutions.mutation_rate = config::ea_1_mutation_rate;
        this->solutions.max_evo_cycles = config::ea_1_max_evo_cycles;
        this->solutions.max_fit_evals = config::ea_1_max_fit_evals;
        
        if(this->variant.isle.id == 0) {
            this->offsets = generate_offsets(-2.5, 2.5, .5);
        }
        
        MPI_Bcast(&this->offsets, DIM, MPI_DOUBLE, 0, this->variant.tcomm);
        
        double init_end = MPI_Wtime();
        
        this->init_duration = (init_end - this->start);
        
    }
    
    //template<typename o> void end(objective<o> &obj) { obj.end(*this); }
    
};

#endif /* dtype_solver_h */

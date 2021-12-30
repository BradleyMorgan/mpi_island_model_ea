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

    int source = mpi.id;
    int locale;
    int migrations = 0;
    int selected = 0;
    int survival = 0;
    
    std::array<char[64], 2> parents;

    int numeric_id();
    
    template<typename i> void log(i &interval);
    template<typename i> void measure(i &interval);
    
    solution() : source(mpi.id), locale(mpi.id) { strcpy(id, uniqid(sinstances++)); }

};

template<typename i> void solution::measure(i &interval) { }
 
int solution::numeric_id() {
    return (int)std::strtol(this->id, nullptr, 10);
}

template<typename i> void solution::log(i &interval) {
 
    
    if(mpi.id==0) {
        
        fprintf(config::ea_1_genome_out, "%s,%d,%s,", interval.name, interval.id, this->id);
        
        for(int n=0; n<DIM; n++) {
            
            fprintf(config::ea_1_genome_out, "%f;", this->input[n]);
            
        }
        
        fprintf(config::ea_1_genome_out, "\r\n");
        
    }

    fflush(config::ea_1_genome_out);
    
}

#endif /* dtype_solution_h */

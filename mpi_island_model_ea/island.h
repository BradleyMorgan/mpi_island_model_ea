//
//  island.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 2/26/20.
//  Copyright © 2020 Bradley Morgan. All rights reserved.
//

#ifndef island_h
#define island_h

#pragma mark FUNCTION: island::calculate::total_fitness()

// calculate the island's total fitness for distribution ...

struct ea;

void island::calculate::total_fitness(island &p) {
    
    p.total_fitness = 0.0;
    
    int i = 0;
    
    for(std::vector<genome>::iterator it = p.population.begin(); it != p.population.end(); ++it) {
        LOG(10, 0, 0, "isle %d solution %d fitness = %f\r\n", p.id, i, it->fitness);
        p.total_fitness += it->fitness;
        //printf("CALC id %llu p1 %llu p2 %llu\r\n", it->id, it->parents[0], it->parents[1]);
        i++;
    }
    
    LOG(6, 0, 0, "island::calculate::total_fitness island %d = %f\r\n", p.id, p.total_fitness);
    
}

#pragma mark FUNCTION: island::calculate::average_fitness()

void island::calculate::average_fitness(island &p) {
    
    island::calculate::total_fitness(p);
    
    p.average_fitness = p.total_fitness / p.population.size();
    
    LOG(6, 0, 0, "island::calculate::average_fitness island %d = %f\r\n", p.id, p.average_fitness);
    
}

#pragma mark FUNCTION: island::migration::receive()

// initiate a receive operation for every island in this island's senders list ...

//void island::migration::receive(island &p, MPI_Comm &comm, int &eval) {
void island::migration::receive(island &p, MPI_Datatype &d, int &eval) {
    
    LOG(5, 0, 0, "island::migration::receive(%d) operations for %lu senders\r\n", p.id, p.senders.size());
    
    //std::array<double, DIM> x;
    
    MPI_Status migrant_status;
    
    for(int i=0; i<p.senders.size(); i++) {
        
        genome x;
        
        int tag = ((p.senders[i]+1)*10000)+eval;
        
        LOG(5, 0, 0, "island %d waiting for migrant recv tag %d from island %d ... \r\n", p.id, tag, p.senders[i]);
        
        MPI_Recv(&x, 1, d, p.senders[i], tag, p.tcomm, &migrant_status);
        
        int top = p.population.size() * 0.20;
        int idx = rand()%top;
        
        x.migrations++;
        x.locale = p.id;
        p.population[idx] = x;
        
        visa v(eval, p.id, p.senders[i], p.population[idx].id);
        p.visas.push_back(v);
            
        LOG(5, 0, 0, "island %d received migrant %s tag=%d from island %d: [%f,%f] with status %d\r\n", p.id, x.id, tag, migrant_status.MPI_SOURCE, p.population[0].fitness, p.population[0].fitness, migrant_status.MPI_ERROR);
    }
    
}
#pragma mark FUNCTION: island::migration::send()

// initiate a send operation for every island in this island's receivers list ...

void island::migration::send(island &p, MPI_Datatype &d, int &eval) {
//void island::migration::send(island &p, MPI_Datatype &d, int &eval) {
    
    LOG(5, 0, 0, "island::migration::send(%d) initiating operations for %lu receivers\r\n", p.id, p.receivers.size());

    //TODO: test different migration selection methods ...
    
    for(int i=0; i<p.receivers.size(); i++) {
        int tag = ((p.id+1)*10000+eval);
        int top = p.population.size() * 0.20;
        int idx = rand()%top;
        //LOG(6, 0, 0, "island %d sending migrant %s to island %d ... \r\n", p.id, p.population[i].id, p.receivers[i]);
        LOG(6, 0, 0, "island %d sending tag=%d sol_id=%s p1=%s p2=%s to island %d\r\n", p.id, tag, p.population[idx].id, p.population[idx].parents[0], p.population[idx].parents[1], p.receivers[i]);
        MPI_Send(&p.population[idx], 1, d, p.receivers[i], tag, p.tcomm);
        LOG(6, 0, 0, "island %d sent tag %d migrant to island %d: [%f,%f]\r\n", p.id, tag, p.receivers[i], p.population[i].input[0], p.population[i].input[1]);
    }

}


#endif /* island_h */

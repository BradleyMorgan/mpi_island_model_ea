//
//  stats.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 3/30/20.
//  Copyright © 2020 Bradley Morgan. All rights reserved.
//

#ifndef stats_h
#define stats_h

#include <mpi.h>

#define RED   "\x1B[31m"
#define GRN   "\x1B[32m"
#define YEL   "\x1B[33m"
#define BLU   "\x1B[34m"
#define MAG   "\x1B[35m"
#define CYN   "\x1B[36m"
#define WHT   "\x1B[37m"
#define RESET "\x1B[0m"

//void log_pop_stats(ea_solver &solver, island &isle, MPI_Datatype &visa_type) {
//    
//    std::vector<visa> visas;
//    unsigned long int vrecsize = 0;
//    unsigned long int vsndsize = isle.visas.size();
//
//    if(mpi.id==0) {
//        for(int i=1; i<config::world_size; i++) {
//            MPI_Recv(&vrecsize, 1, MPI_INT, i, i, isle.tcomm, MPI_STATUS_IGNORE);
//            std::vector<visa> v;
//            v.resize(vrecsize);
//            MPI_Recv(&v[0], int(vrecsize), visa_type, i, i*10, isle.tcomm, MPI_STATUS_IGNORE);
//            visas.insert(visas.end(), v.begin(), v.end());
//        }
//        visas.insert(visas.end(), isle.visas.begin(), isle.visas.end());
//    } else {
//        MPI_Send(&vsndsize, 1, MPI_INT, 0, mpi.id, isle.tcomm);
//        MPI_Send(&isle.visas[0], int(vsndsize), visa_type, 0, mpi.id*10, isle.tcomm);        
//    }
//    
//    if(mpi.id==0) {
//            
//        std::map<long long int, int> groups_e5;
//        std::map<long long int, int> groups_e4;
//        std::map<long int, int> groups_e3;
//        
//        // top 20 individuals to cut log size
//        
//        for (auto it = solver.solutions.population.begin(); it !=solver.solutions.population.begin() + 20; ++it) {
//            
//            long long int e5 = (it->group * 1000000000);
//            ++groups_e5[e5];
//            
//            long long int e4 = (it->group * 1000000);
//            ++groups_e4[e4];
//            
//            long int e3 = (it->group * 1000);
//            ++groups_e3[e3];
//            
//            //std::sort(sol_stats.begin(), sol_stats.end());
//            
//            char sol[DIM*sizeof(double)];
//            
//            int offs = 0;
//            
//            for(int i=0; i<DIM; i++) {
//                char delim = (i == DIM-1) ? '\0' : ';';
//                offs += snprintf(sol+offs, sizeof(sol)>offs?sizeof(sol)-offs:0, "%f%c", it->input[i], delim);
//            }
//    
//            char sol_id[sizeof(it->id)*sizeof(char)];
//            strcpy(sol_id, it->id);
//            
//            std::vector<visa> sol_visas;
//            
//            std::copy_if(visas.begin(), visas.end(), std::back_inserter(sol_visas), [&](visa &v) { return strcmp(sol_id, v.genome_id) == 0; });
//            
//            char vis[sizeof(int)*sizeof(sol_visas)+sizeof(char)*sizeof(sol_visas)];
//            int offv = 0;
//            
//            for(int v=0; v<sol_visas.size(); v++) {
//                char delim = (v == sol_visas.size()-1) ? '\0' : ';';
//                offv += snprintf(vis+offv, sizeof(vis)>offv?sizeof(vis)-offv:0, "%d%c", sol_visas[v].destination, delim);
//            }
//            
//            std::fprintf(config::ea_1_population_out, "%d," "%d," "%s," "%d," "%d," "%s," "%s," "%d," "%d,", solver.solutions.run.id, solver.solutions.run.cycle.eval.id, it->id, it->source, it->locale, it->parents[0], it->parents[1], it->selected, it->survival);
//                         
//            std::fprintf(config::ea_1_population_out, "%lld," "%d,"  "%lld," "%d," "%ld," "%d,", e5, groups_e5[e5], e4, groups_e4[e4], e3, groups_e3[e3]);
//                         
//            std::fprintf(config::ea_1_population_out, "%f," "%f," "%d," "%s," "%s\r\n", it->fitness, it->selection_distribution, it->migrations, vis, sol);
//            
//        }
//        
//    }
//    
//}

#endif /* stats_h */

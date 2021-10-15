//
//  dtype_ea.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 10/13/21.
//  Copyright © 2021 Bradley Morgan. All rights reserved.
//

#ifndef dtype_ea_h
#define dtype_ea_h

#pragma mark DATATYPE: @ea_meta{}

// describes the properties needed for parallel (island model) EA

struct ea_variant {

    unsigned long start;

    int islands;
    int island_size;
    int root = 0;

    // isle holds context specific population and migration data
    // from the perspective of *this* island (mpi rank)
    // see @island{} datatype
    
    island isle;
    
    // TODO: maintaining a list of island_ids doesn't seem necessary
    
    std::vector<int> island_ids;
    
    MPI_Datatype solution_type;
    MPI_Datatype visa_type;
    MPI_Datatype topology_type;
    MPI_Datatype stats_type;
    MPI_Datatype comm_type;
    
    MPI_Comm tcomm;
    
};

//template<typename o> struct ea_objectives {
//    //objective<genome> objective;
//    std::vector<o> genomes;
//};

#pragma mark DATATYPE: @ea{}

// wrapper type to serve as the trunk for ea properties

struct ea {
    
    ea_run run;
    ea_variant variant;
    
    unsigned long int start = 0;
    unsigned long int duration = 0;
    
    unsigned long int init_start = 0;
    unsigned long int init_duration = 0;
    
    template<typename f> void init(f function) { function(*this); }
    template<typename f, typename v> void init(f function, v &variant) { function(*this, variant); }
    
    //template<typename o, typename v, typename genome> void begin(ea_run &run, v &variant, objective<genome> &obj) { run.begin(*this, obj); }
    template<typename o> void begin(ea_run &run, objective<o> &obj) { run.begin(*this, obj); log_begin(obj.run, obj, *this); }
    template<typename o> void end(ea_run &run, objective<o> &obj) { run.end(*this, obj); log_end(obj.run, obj, *this); }
    
    template<typename o> void begin(ea_eval &eval, objective<o> &obj) { eval.begin(); log_begin(obj.run.eval, obj, *this);  }
    template<typename o> void end(ea_eval &eval, objective<o> &obj) { eval.end(); log_end(obj.run.eval, obj, *this); }
    
    template<typename o, typename g> void log(ea_run &run, objective<o> &obj, ea &meta, g &genome) { log(run, obj, *this, meta, genome); }
    template<typename o, typename g> void log(ea_eval &eval, objective<o> &obj, ea &meta, g &genome) { log(eval, obj, *this, meta, genome); }
    
    template<typename o, typename g> void log(objective_run &run, objective<o> &obj, ea &meta, g &genome) { obj.log(run, obj, *this, meta, genome); }
    template<typename o, typename g> void log(objective_eval &eval, objective<o> &obj, ea &meta, g &genome) { obj.log(eval, obj, *this, meta, genome); }
    
//    template<typename o> void begin(ea_eval &eval, objective<o> &obj) { this->run.eval.begin(*this, obj); }
//    template<typename o> void end(ea_eval &eval, objective<o> &obj) { this->run.eval.end(*this, obj); }
    
//    template<typename o> void begin(ea_eval &eval, objective<o> &obj) { eval.begin(*this, obj); }
//    template<typename o> void end(ea_eval &eval, objective<o> &obj) { eval.end(*this, obj); }
//
//    template<typename f, typename o> void begin(ea_run &run, objective<o> &obj) { run.begin(*this, obj); }
//    template<typename f, typename o> void end(ea_run &run, objective<o> &obj) { run.end(*this, obj); }
        
//    template<typename f> void run_init(f function) { function(*this); }
//    template<typename f, typename v> void run_init(f function, v &variant) { function(*this, variant); }
//    
//    template<typename f> void run_end(f function) { function(*this); }
//    template<typename f, typename v> void run_end(f function, v &variant) { function(*this, variant); }
//    
//    template<typename f> void eval_init(f function) { function(*this); }
//    template<typename f, typename v> void eval_init(f function, v &variant) { function(*this, variant); }
//    
//    template<typename f> void eval_end(f function) { function(*this); }
//    template<typename f, typename v> void eval_end(f function, v &variant) { function(*this, variant); }
    
    std::array<double, DIM> offsets;
    
    // TODO: generalize collection type to store multiple objectives
    
    objective<topology> topologies;
    objective<solution> solutions;
    
//    void ea_end() {
//        
//        LOG(8, 0, 0, "island %d end EA\r\n", this->variant.isle.id);
//        
//        if(this->variant.isle.id == 0) {
//
//            char canary[128];
//
//            sprintf(canary, "%s/end.txt", config::logs_subpath);
//            FILE *eaend = fopen(canary, "w");
//
//            fprintf(eaend, "ended at %lu", time(0));
//
//            fclose(eaend);
//
//        }
//        
//        MPI_Finalize();
//        
//    }
    
    template<typename f, typename o> void populate(objective<o> &obj) { obj.populate(*this); }
    template<typename f, typename o> void populate(objective<o> &obj, f function) { obj.populate(*this, function); }
    template<typename f, typename o, typename v> void populate(objective<o> &obj, v &variant, f function) { obj.populate(*this, variant, function); }
    template<typename f, typename o, typename v, typename m> void populate(objective<o> &obj, v &variant, m &meta, f function) { obj.populate(*this, variant, meta, function); }
    
    template<typename f, typename o> void distribute(objective<o> &obj, f function) { obj.distribute(*this, function); }
    template<typename f, typename o, typename v> void distribute(objective<o> &obj, v &variant, f function) { obj.distribute(*this, variant, function); }
    
    template<typename f, typename o, typename g> void evaluate(objective<o> &obj, g &individual, f function) { obj.evaluate(*this, individual, function); }
    
    template<typename o> void evolve(objective<o> &obj) { obj.evolve(*this); }
    template<typename f, typename o> void evolve(objective<o> &obj, f function) { obj.evolve(*this, function); }
    template<typename f, typename o, typename m> void evolve(objective<o> &obj, m &meta, f function) { obj.evolve(*this, meta, function); }
    template<typename f, typename o, typename m, typename g> void evolve(objective<o> &obj, m &meta, g &individual, f function) { obj.evolve(*this, meta, individual, function); }
    template<typename f, typename o, typename v, typename m, typename g> void evolve(objective<o> &obj, v &variant, m &meta, g &individual, f function) { obj.evolve(*this, variant, meta, individual, function); }
    
    template<typename o> void end(objective<o> &obj) { obj.end(*this, obj); }
    template<typename o> void end(ea &ea, objective<o> &obj) { ea_end(*this, obj); }
    
};

void log_eval() {
    
    //void meta_run_end(ea &meta) {
    //
    //    if(meta.variant.isle.id != 0) { return; }
    //
    //    LOG(2, meta.variant.isle.id, 0, "\r\n--- END META RUN %d ---\r\n", meta.run.id);
    //
    //    double run_end = MPI_Wtime();
    //
    //    meta.run.stats.run_duration = run_end - meta.run.start;
    //
    //    fprintf(config::topo_run_stats_out, "average_topo_fitness, global_best_topo_id, global_best_topo_rounds, global_best_topo_channels, global_best_topo_round_fitness, global_best_topo_fitness1, local_best_topo_fitness, global_best_topo_fitness2, average_local_best_topo_fitness, average_global_best_topo_fitness, t_id, t_rounds, t_channels, t_fitness\r\n");
    //
    //    std::fprintf(config::topo_run_stats_out, "%d,%f,%f,%f,%f,%f,%d", meta.run.id, meta.run.stats.run_duration, meta.run.eval.stats.average_local_best_topo_fitness, meta.run.eval.stats.average_global_best_topo_fitness, meta.run.eval.stats.global_best_topo_fitness, meta.run.eval.stats.total_migrate_time, meta.run.stats.total_channels);
    //
    //    fflush(config::topo_run_stats_out);
    //
    //}

    
}

#endif /* dtype_ea_h */

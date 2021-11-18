//
//  config.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 2/26/20.
//  Copyright © 2020 Bradley Morgan. All rights reserved.
//

#ifndef config_h
#define config_h

#include <fstream>
#include <map>
#include <sys/stat.h>
#include <sys/time.h>

static const int DIM = 10;

static unsigned long long int sinstances = 0;
static unsigned long long int tinstances = 0;

# pragma mark DATATYPE @config{}

// global configuration items

namespace config {

    std::map<std::string, std::string> items;

    FILE *log_out;
    FILE *sol_stats_out;
    FILE *topo_stats_out;
    FILE *topo_run_stats_out;
    FILE *run_stats_out;
    FILE *solution_out;
    FILE *topo_out;
    FILE *solpop_out;

    int dim;
    int world_size;
    int seed = 0;
    int migration_cap = 0;
    int send_cap = 0;
    int ea_mode = 0;
    int migration_interval = 1;
    int island_mu = 0;
    int island_lambda = 0;
    int mu_mode = 0;

    int ea_1_runs = 0;
    int ea_1_lambda = 0;
    int ea_1_mu = 0;
    int ea_1_max_evo_cycles = 0;
    int ea_1_max_fit_evals = 0;
    int ea_1_log_interval = 0;
    int ea_1_population_log_interval = 0;

    double ea_1_mutation_rate = 0.0;

    int ea_2_runs = 0;
    int ea_2_mu = 0;
    int ea_2_lambda = 0;
    int ea_2_max_evo_cycles = 0;
    int ea_2_max_fit_runs = 0;
    int ea_2_max_fit_evals = 0;
    int ea_2_log_interval = 0;
    int ea_2_population_log_interval = 0;

    double ea_2_mutation_rate = 0.0;

    double island_lambda_pct = 0.0;
    double sparsity = 0.0;

    char log_fname[100];
    char stats_fname[100];
    char solution_fname[100];
    char topo_fname[100];
    char topo_run_fname[100];
    char run_stats_fname[100];
    char solpop_fname[100];
    char logs_subpath[100];
    char stats_subpath[100];
    char topos_subpath[100];
    char solpop_subpath[100];

    char mu_msg[64];
    char lambda_msg[64];
    char subpop_msg[64];

    std::pair<int,int> id_field_prefix1;

    void load(const char *input, const int world_size, const int world_rank);

};

void config::load(const char *input, const int world_size, const int world_rank) {
    
    if(input == NULL) { input = (char *)"config.txt"; }
    
    std::ifstream config_file(input);

    std::string line;
    std::string key;
    std::string value;
    
    while (std::getline(config_file, line)) {
        
        line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
        
        if(line[0] == '#' || line.empty() || line.length()==0 ) { continue; }
        
        unsigned long int delimiterPos = line.find(':');
        key = line.substr(0, delimiterPos);
        value = line.substr(delimiterPos + 1);
        
        if(key != "") {
           //printf("%s => %s\r\n", key.c_str(), value.c_str());
           config::items[key] = value;
       }
    
    }
    
    config::world_size = world_size;
    config::dim = stoi(config::items["dim"]);
    config::sparsity = stod(config::items["sparsity"]);
    config::migration_cap = stoi(config::items["migration_cap"]);
    config::send_cap = stoi(config::items["migration_cap"]);
    config::ea_mode = stoi(config::items["ea_mode"]);
    config::migration_interval = stoi(config::items["migration_interval"]);
    config::mu_mode = stoi(config::items["mu_mode"]);
    
    config::ea_1_runs = stoi(config::items["ea_1_runs"]);
    config::ea_1_mu = stoi(config::items["ea_1_mu"]);
    config::ea_1_lambda = stoi(config::items["ea_1_lambda"]);
    config::ea_1_mutation_rate = stod(config::items["ea_1_mutation_rate"]);
    config::ea_1_max_evo_cycles = stoi(config::items["ea_1_max_evo_cycles"]);
    config::ea_1_max_fit_evals = stoi(config::items["ea_1_max_fit_evals"]);
    config::ea_1_log_interval = stoi(config::items["ea_1_log_interval"]);
    config::ea_1_population_log_interval = stoi(config::items["ea_1_population_log_interval"]);
    
    config::ea_2_runs = stoi(config::items["ea_2_runs"]);
    config::ea_2_mu = stoi(config::items["ea_2_mu"]);
    config::ea_2_lambda = stoi(config::items["ea_2_lambda"]);
    config::ea_2_mutation_rate = stod(config::items["ea_2_mutation_rate"]);
    config::ea_2_max_evo_cycles = stoi(config::items["ea_2_max_evo_cycles"]);
    config::ea_2_max_fit_evals = stoi(config::items["ea_2_max_fit_evals"]);
    config::ea_2_max_fit_runs = stoi(config::items["ea_2_max_fit_runs"]);
    config::ea_2_log_interval = stoi(config::items["ea_2_log_interval"]);
    config::ea_2_population_log_interval = stoi(config::items["ea_2_population_log_interval"]);
    
    sprintf(config::mu_msg, ": ");
    sprintf(config::lambda_msg,": ");
    sprintf(config::subpop_msg, ": ");
    
    // mu_mode
    // 0: fixed global mu, comm size based island mu, fixed island lambda
    // 1: comm size based global mu, fixed island mu, fixed island lambda
    // 2: comm size based global mu, fixed island mu, percentage island lambda
    //
    // global lambda is fairly insignificant as the current implementation
    // performs crossover at the island level, but including it for placeholder
    
    if(config::mu_mode == 1) {
        config::island_mu = stoi(config::items["island_mu"]);
        config::island_lambda = stod(config::items["island_lambda"]);
        config::ea_1_mu = world_size * config::island_mu;
    } else if(config::mu_mode == 2) {
        config::island_mu = stoi(config::items["island_mu"]);
        config::island_lambda = config::island_mu * stod(config::items["island_lambda_pct"]);
        config::ea_1_mu = world_size * config::island_mu;
    } else {
        config::ea_1_mu = stoi(config::items["mu"]);
        config::ea_1_lambda = stoi(config::items["lambda"]);
        config::island_mu = config::ea_1_mu / world_size;
        config::island_lambda = stod(config::items["island_lambda"]);   
    }
    
    char mode[5];
    
    if(config::ea_mode == 0) {
        sprintf(mode, "%s", "bench");
    } else {
        sprintf(mode, "%s", "evo");
    }
    
    char mumode[8];
    
    if(config::mu_mode < 1) {
        sprintf(mumode, "%s", "fixed" );
    } else {
        sprintf(mumode, "%s", "dynamic" );
    }
    
    // create unique pathing for the supplied configuration so we can keep better track of results ...
    
    if(world_rank == 0) {
        
        // collect current time data
        
        time_t rawtime;
        struct tm * timeinfo;

        time (&rawtime);
        timeinfo = localtime (&rawtime);
        
        char dstr[6];
        
        strftime(dstr, 80, "%m%d%y", timeinfo);
        
        // create log path heirarchy to organize output by run date, communicator size, and ea mode
        
        //char sub1[12];
        char sub2[12];
        char sub3[58];
        //char sub4[68];
        
        //sprintf(sub1, "logs/%s", dstr);
        //mkdir(sub1,0740);
        //sprintf(sub2, "%s/%d", sub1, world_size);
        sprintf(sub2, "logs/%d", world_size);
        mkdir(sub2,0740);
        sprintf(sub3, "%s/%s", sub2, mode);
        mkdir(sub3,0740);
        //sprintf(sub4, "%s/%s", sub3, mumode);
        //mkdir(sub4,0740);
        
        int fn = 1;
        bool found = true;
        while(found) {
            sprintf(config::logs_subpath, "%s/%s_%03d", sub3, config::items["config_name"].c_str(), fn);
            std::ifstream ifile;
            ifile.open(config::logs_subpath);
            if(ifile) {
                fn++;
            } else {
                found = false;
            }
        }
        
        mkdir(config::logs_subpath, 0740);
        
        sprintf(config::stats_subpath, "%s/stats", logs_subpath);
        mkdir(config::stats_subpath, 0740);
    
        // open a log file with the parsed configuration values ...
        
        sprintf(config::log_fname, "%s/%s_%d_%s.txt", config::logs_subpath, config::items["log_file"].c_str(), world_size, mode);
        config::log_out = fopen(config::log_fname, "w");
        
    }
    
    // seed the rng ...
    
    if(stoi(config::items["seed"]) == 0) {
        timeval time;
        gettimeofday(&time, NULL);
        config::seed = (unsigned int)(time.tv_sec * 1000) + time.tv_usec;
    } else {
        config::seed = stoi(config::items["seed"]);
    }
                            
    srand(config::seed);
    
    if(world_rank == 0) {
    
        fprintf(config::log_out, "seed: %u\r\n", config::seed);
    
        std::map <std::string, std::string>::iterator items_iterator;
        
        for(items_iterator = std::next(items.begin()); items_iterator != items.end(); ++items_iterator) {
            fprintf(config::log_out, "%s: %s\r\n", items_iterator->first.c_str(), items_iterator->second.c_str());
        }
        
        if(config::mu_mode == 0) {
            sprintf(subpop_msg, "(calculated: %d/%d): ", config::ea_1_mu, world_size);
        }
        
        if(config::mu_mode > 0) {
           sprintf(mu_msg, "(calculated: %d*%d): ", world_size, config::island_mu);
        }
        
        if(config::mu_mode == 2) {
            sprintf(lambda_msg, "(calculated: %d*%f): ", config::island_mu, stod(config::items["island_lambda_pct"]));
        }
        
        fprintf(config::log_out, "log file: %s\r\n", config::log_fname);
        fprintf(config::log_out, "world size: %d\r\n", world_size);
        fprintf(config::log_out, "subpopulation size %s %d\r\n", config::subpop_msg, config::island_mu);
        fprintf(config::log_out, "global mu %s %d\r\n", config::mu_msg, config::ea_1_mu);
        fprintf(config::log_out, "island lambda %s %d\r\n", config::lambda_msg, stoi(config::items["island_lambda"]));
        
    }
    
    if(world_rank == 0) {
    
        sprintf(config::stats_fname, "%s/%s_sol_%d_%s.csv", config::stats_subpath, config::items["stats_file"].c_str(), world_size, mode);
        config::sol_stats_out = fopen(config::stats_fname, "w");
        
        sprintf(config::topo_fname, "%s/%s_topo_%d_%s.csv", config::stats_subpath, config::items["topo_file"].c_str(), world_size, mode);
        config::topo_stats_out = fopen(config::topo_fname, "w");
        
        sprintf(config::topo_run_fname, "%s/%s_topo_run_%d_%s.csv", config::stats_subpath, config::items["topo_file"].c_str(), world_size, mode);
        config::topo_run_stats_out = fopen(config::topo_run_fname, "w");
        
        sprintf(config::run_stats_fname, "%s/%s_run_%d_%s.csv", config::stats_subpath, config::items["stats_file"].c_str(), world_size, mode);
        config::run_stats_out = fopen(config::run_stats_fname, "w");
        
        sprintf(config::solpop_fname, "%s/%s_%d_%s.csv", config::stats_subpath, "solpop", world_size, mode);
        config::solpop_out = fopen(config::solpop_fname, "w");
        
        //sprintf(config::solution_fname, "%s/%s_solution_%d_%ld.txt", config::stats_subpath, config::items["stats_file"].c_str(), world_size, time(0));
        //config::solution_out = fopen(config::solution_fname, "w");
        
        fprintf(config::sol_stats_out, "run,cycle,eval,average_fitness,local_best_fitness,global_best_fitness,average_local_best_fitness,average_global_best_fitness,topo_id,topo_fitness,topo_channels,average_scatter_time,average_gather_time,average_migrate_time,init_duration,ea_elapsed_t,run_elapsed_t,cycle_elapsed_t,avg_cycle_t,eval_elapsed_t\r\n");
        
        fprintf(config::topo_stats_out, "run, cycle, eval, average_topo_fitness, global_best_topo_id, global_best_topo_rounds, global_best_topo_channels, global_best_topo_round_fitness, global_best_topo_fitness1, local_best_topo_fitness, global_best_topo_fitness2, average_local_best_topo_fitness, average_global_best_topo_fitness, t_id, t_rounds, t_channels, t_fitness, ea_elapsed_t, run_elapsed_t, cycle_elapsed_t, eval_elapsed_t, migration_t, avg_migr_time, total_cycle_t, avg_cycle_t\r\n");
        
        fprintf(config::topo_run_stats_out, "run, cycle, eval, average_topo_fitness, global_best_topo_id, global_best_topo_fitness, global_best_topo_rounds, global_best_topo_channels, global_best_topo_round_fitness, average_local_best_topo_fitness, average_global_best_topo_fitness, total_channels, ea_elapsed_t, run_elapsed_t, cycle_elapsed_t, eval_elapsed_t, migration_t, avg_migration_t\r\n");
        
        fprintf(config::run_stats_out, "run,cycle,eval,global_best_fitness,average_local_best_fitness,average_global_best_fitness,total_scatter_time,total_gather_time,total_migration_time,run_duration,init_duration,world_size,subpopulation_size,ea_elapsed_t,run_elapsed_t,cycle_elapsed_t,avg_cycle_t,eval_elapsed_t\r\n\r\n");
        
        fprintf(config::solpop_out, "run,cycle,eval,id,origin,locale,parent1,parent2,pselected,survival,10e5_fit_group,10e5_fit_count,10e4_fit_group,10e4_fit_count,10e3_fit_group,10e3_fit_count,fitness,selection_dist,migration_count,visas,genes\r\n");
        
        fprintf(config::log_out, "stats file: %s\r\n", config::stats_fname);
        fprintf(config::log_out, "run stats file: %s\r\n", config::run_stats_fname);
        fprintf(config::log_out, "topology stats file: %s\r\n", config::topo_fname);
        fprintf(config::log_out, "topology run stats file: %s\r\n", config::topo_run_fname);
        
        sprintf(topos_subpath, "%s/%s", config::stats_subpath, "topos");
        
        mkdir(topos_subpath, 0775);
    
        fflush(config::log_out);
        fflush(config::topo_stats_out);
        fflush(config::topo_run_stats_out);
        
    }
    
    // associate each mpi rank (island) with
    
    config::id_field_prefix1 = std::make_pair(world_rank, 3);
    
}

#endif /* config_h */

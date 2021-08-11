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
  
namespace config {

    std::map<std::string, std::string> items;

    FILE *log_out;
    FILE *sol_stats_out;
    FILE *topo_stats_out;
    FILE *run_stats_out;
    FILE *solution_out;
    FILE *topo_out;
    FILE *solpop_out;

    int dim;
    int world_size;
    int evals = 0;
    int runs = 0;
    int lambda = 0;
    int mu = 0;
    int topo_lambda = 0;
    int topo_mu = 0;
    int seed = 0;
    int migration_cap = 0;
    int send_cap = 0;
    int topo_evals = 0;
    int ea_mode = 0;
    int migration_interval = 1;
    int island_mu = 0;
    int island_lambda = 0;
    int mu_mode = 0;

    double island_lambda_pct = 0.0;
    double mutation_rate = 0.0;
    double topo_mutation_rate = 0.0;
    double sparsity = 0.0;

    char log_fname[100];
    char stats_fname[100];
    char solution_fname[100];
    char topo_fname[100];
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

    std::string key;
    std::string value;
    
    while (config_file.good()) {
       
       getline(config_file, key, ':');
       getline(config_file, value, '\n');
       
       if(key != "") {
           //printf("%s => %s\r\n", key.c_str(), value.c_str());
           config::items[key] = value;
       }
       
    }
    
    config::world_size = world_size;
    config::dim = stoi(config::items["dim"]);
    config::runs = stoi(config::items["runs"]);
    config::evals = stoi(config::items["evals"]);
    config::mutation_rate = stod(config::items["mutation_rate"]);
    config::sparsity = stod(config::items["sparsity"]);
    config::migration_cap = stoi(config::items["migration_cap"]);
    config::send_cap = stoi(config::items["migration_cap"]);
    config::topo_evals = stoi(config::items["topo_evals"]);
    config::ea_mode = stoi(config::items["ea_mode"]);
    config::migration_interval = stoi(config::items["migration_interval"]);
    config::topo_mutation_rate = stod(config::items["topo_mutation_rate"]);
    config::mu_mode = stoi(config::items["mu_mode"]);
    config::topo_mu = stoi(config::items["topo_mu"]);
    config::topo_lambda = stoi(config::items["topo_lambda"]);
        
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
        config::mu = world_size * config::island_mu;
    } else if(config::mu_mode == 2) {
        config::island_mu = stoi(config::items["island_mu"]);
        config::island_lambda = config::island_mu * stod(config::items["island_lambda_pct"]);
        config::mu = world_size * config::island_mu;
    } else {
        config::mu = stoi(config::items["mu"]);
        config::lambda = stoi(config::items["lambda"]);
        config::island_mu = config::mu / world_size;
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
        sprintf(mumode, "%s", "fixedmu" );
    } else {
        sprintf(mumode, "%s", "dynmu" );
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
        
        char sub1[12];
        char sub2[4];
        char sub3[58];
        char sub4[68];
        
        sprintf(sub1, "logs/%s", dstr);
        mkdir(sub1,0740);
        sprintf(sub2, "%s/%d", sub1, world_size);
        mkdir(sub2,0740);
        sprintf(sub3, "%s/%s", sub2, mode);
        mkdir(sub3,0740);
        sprintf(sub4, "%s/%s", sub3, mumode);
        mkdir(sub4,0740);
        
        sprintf(config::logs_subpath, "%s/%s_%s_%ld", sub4, config::items["config_name"].c_str(), mode, time(0));
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
            sprintf(subpop_msg, "(calculated: %d/%d): ", config::mu, world_size);
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
        fprintf(config::log_out, "global mu %s %d\r\n", config::mu_msg, config::mu);
        fprintf(config::log_out, "island lambda %s %d\r\n", config::lambda_msg, stoi(config::items["island_lambda"]));
        
    }
    
    if(world_rank == 0) {
    
        sprintf(config::stats_fname, "%s/%s_sol_%d_%s.csv", config::stats_subpath, config::items["stats_file"].c_str(), world_size, mode);
        config::sol_stats_out = fopen(config::stats_fname, "w");
        
        sprintf(config::topo_fname, "%s/%s_topo_%d_%s.csv", config::stats_subpath, config::items["topo_file"].c_str(), world_size, mode);
        config::topo_stats_out = fopen(config::topo_fname, "w");
        
        sprintf(config::run_stats_fname, "%s/%s_run_%d_%s.csv", config::stats_subpath, config::items["stats_file"].c_str(), world_size, mode);
        config::run_stats_out = fopen(config::run_stats_fname, "w");
        
        sprintf(config::solpop_fname, "%s/%s_%d_%s.csv", config::stats_subpath, "solpop", world_size, mode);
        config::solpop_out = fopen(config::solpop_fname, "w");
        
        //sprintf(config::solution_fname, "%s/%s_solution_%d_%ld.txt", config::stats_subpath, config::items["stats_file"].c_str(), world_size, time(0));
        //config::solution_out = fopen(config::solution_fname, "w");
        
        fprintf(config::sol_stats_out, "run,eval,average_fitness,local_best_fitness,global_best_fitness,average_local_best_fitness,average_global_best_fitness,average_scatter_time,average_gather_time,average_migrate_time,init_duration,eval_duration\r\n");
        
        fprintf(config::topo_stats_out, "average_topo_fitness, global_best_topo_id, global_best_topo_rounds, global_best_topo_channels, global_best_topo_round_fitness, global_best_topo_fitness1, local_best_topo_fitness, global_best_topo_fitness2, average_local_best_topo_fitness, average_global_best_topo_fitness, t_id, t_rounds, t_channels, t_fitness\r\n");
        
        fprintf(config::run_stats_out, "run,global_best_fitness,average_local_best_fitness,average_global_best_fitness,total_scatter_time,total_gather_time,total_migration_time,run_duration,init_duration,world_size,subpopulation_size, global_best_topo_fitness, average_local_best_topo_fitness, average_global_best_topo_fitness\r\n");
        
        fprintf(config::solpop_out, "run, eval, id, origin, locale, parent1, parent2, pselected, survival, 10e5_fit_group, 10e5_fit_count, 10e4_fit_group, 10e4_fit_count, 10e3_fit_group, 10e3_fit_count, fitness, selection_dist, migration_count, visas, genes\r\n");
        
        fprintf(config::log_out, "stats file: %s\r\n", config::stats_fname);
        fprintf(config::log_out, "run stats file: %s\r\n", config::run_stats_fname);
        fprintf(config::log_out, "topology files: %s\r\n", config::topo_fname);
        
        sprintf(topos_subpath, "%s/%s", config::stats_subpath, "topos");
        
        mkdir(topos_subpath, 0775);
    
        fflush(config::log_out);
        
    }
    
    config::id_field_prefix1 = std::make_pair(world_rank, 3);
    
}

#endif /* config_h */

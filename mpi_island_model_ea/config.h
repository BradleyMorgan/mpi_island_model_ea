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

# pragma mark PROGRAM DATATYPE @config{}

// global configuration items

namespace config {

    std::map<std::string, std::string> items;

    // log handles

    FILE *ea_out;
    FILE *system_out;
    FILE *program_out;

    FILE *ea_1_run_out;
    FILE *ea_1_cycle_out;
    FILE *ea_1_eval_out;
    FILE *ea_1_population_out;
    FILE *ea_1_genome_out;

    FILE *ea_2_run_out;
    FILE *ea_2_cycle_out;
    FILE *ea_2_eval_out;
    FILE *ea_2_population_out;
    FILE *ea_2_genome_out;

    // log filenames

    char ea_log_out[64];
    char system_log_out[64];
    char program_log_out[64];
    
    char ea_1_stats_run[64];
    char ea_1_stats_cycle[64];
    char ea_1_stats_eval[64];
    char ea_1_stats_population[64];
    char ea_1_stats_genome[64];
    
    char ea_2_stats_run[64];
    char ea_2_stats_cycle[64];
    char ea_2_stats_eval[64];
    char ea_2_stats_population[64];
    char ea_2_stats_genome[64];
    
    // log subpaths

    char logs_subpath[64];
    char stats_subpath_ea_1_o1[64];
    char stats_subpath_ea_2_o1[64];
    char population_subpath[64];
    char genome_subpath[64];

    // program control

    int ea_mode = 0;

    // algorithm 1

    char ea_1_name[128];

    int ea_1_max_runs = 0;
    int ea_1_max_cycles = 0;
    int ea_1_max_evals = 0;

    int ea_1_log_run_interval = 0;
    int ea_1_log_cycle_interval = 0;
    int ea_1_log_eval_interval = 0;

    int ea_1_log_population_interval = 0;
    int ea_1_log_island_interval = 0;
    int ea_1_log_genome_interval = 0;

    // algorithm 1, objective 1

    char ea_1_o1_name[128];

    int ea_1_o1_max_runs = 0;
    int ea_1_o1_max_cycles = 0;
    int ea_1_o1_max_evals = 0;
    int ea_1_o1_max_fitness_evals = 0;

    int ea_1_o1_log_run_interval = 0;
    int ea_1_o1_log_cycle_interval = 0;
    int ea_1_o1_log_eval_interval = 0;
    int ea_1_o1_log_population_interval = 0;
    int ea_1_o1_log_island_interval = 0;
    int ea_1_o1_log_genome_interval = 0;

    int ea_1_mu = 0;
    int ea_1_lambda = 0;

    double ea_1_mutation_rate = 0.0;

    // algorithm 2

    char ea_2_name[128];

    int ea_2_max_runs = 0;
    int ea_2_max_cycles = 0;
    int ea_2_max_evals = 0;

    int ea_2_log_run_interval = 0;
    int ea_2_log_cycle_interval = 0;
    int ea_2_log_eval_interval = 0;
    int ea_2_log_population_interval = 0;
    int ea_2_log_island_interval = 0;
    int ea_2_log_genome_interval = 0;

    // algorithm 2, objective 1

    char ea_2_o1_name[128];

    int ea_2_o1_max_runs = 0;
    int ea_2_o1_max_cycles = 0;
    int ea_2_o1_max_evals = 0;
    int ea_2_o1_max_fitness_evals = 0;
    
    int ea_2_o1_log_run_interval = 0;
    int ea_2_o1_log_cycle_interval = 0;
    int ea_2_o1_log_eval_interval = 0;
    int ea_2_o1_log_population_interval = 0;
    int ea_2_o1_log_island_interval = 0;
    int ea_2_o1_log_genome_interval = 0;

    int ea_2_mu = 0;
    int ea_2_lambda = 0;
    
    double ea_2_mutation_rate = 0.0;

    // parallel

    int world_size = 0;
    int recv_cap = 0;
    int send_cap = 0;
    int mu_mode = 0;
    int mu_sub = 0;
    int lambda_sub = 0;
    int migration_interval = 1;

    // variant

    int seed = 0;
    int dim = 0;

    double sparsity = 0.0;

    std::pair<int,int> id_field_prefix1;

    // config output strings

    char mu_log[160];
    char sub_log[160];

    void load(const char *input, const int world_size, const int world_rank);

};

# pragma mark PROGRAM DATATYPE @mpi_data{}

// for parallel EAs using message passing interface

// initialize the MPI runtime environment, communicator properties
// each island (rank) given a unique id

struct mpi_data {

    int id, size, name;
    char host[MPI_MAX_PROCESSOR_NAME];
    
    mpi_data() {
          MPI_Init(NULL, NULL);
          MPI_Comm_size(MPI_COMM_WORLD, &size);
          MPI_Comm_rank(MPI_COMM_WORLD, &id);
          MPI_Get_processor_name(host, &name);
    }
    
    ~mpi_data() {
          MPI_Finalize();
    }
    
};

# pragma mark PROGRAM GLOBAL @mpi_data{}

// global, single instance for mpi runtime properties

extern mpi_data mpi;

mpi_data mpi;


#pragma mark PROGRAM CONFIG: FUNCTION: config::load()

void config::load(const char *input, const int world_size, const int world_rank) {
    
    mpi.id = world_rank;
    mpi.size = world_size;
    
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
    
    // program benchmarking -> 0: experiment, 1: control
    
    config::ea_mode = stoi(config::items["ea_mode"]);

    // algorithm 1
    
    sprintf(config::ea_1_name, "%s", config::items["ea_1_name"].c_str());
    
    config::ea_1_max_runs = stoi(config::items["ea_1_max_runs"]);
    config::ea_1_max_cycles = stoi(config::items["ea_1_max_cycles"]);
    config::ea_1_max_evals = stoi(config::items["ea_1_max_evals"]);
    
    config::ea_1_mu = stoi(config::items["ea_1_mu"]);
    config::ea_1_lambda = stoi(config::items["ea_1_lambda"]);
    config::ea_1_mutation_rate = stod(config::items["ea_1_mutation_rate"]);
    
    config::ea_1_log_run_interval = stoi(config::items["ea_1_log_run_interval"]);
    config::ea_1_log_cycle_interval = stoi(config::items["ea_1_log_cycle_interval"]);
    config::ea_1_log_eval_interval = stoi(config::items["ea_1_log_eval_interval"]);
    config::ea_1_log_population_interval = stoi(config::items["ea_1_log_population_interval"]);
    config::ea_1_log_island_interval = stoi(config::items["ea_1_log_island_interval"]);
    config::ea_1_log_genome_interval = stoi(config::items["ea_1_log_genome_interval"]);
    
    sprintf(config::ea_1_o1_name, "%s", config::items["ea_1_o1_name"].c_str());
    
    config::ea_1_o1_max_runs = stoi(config::items["ea_1_o1_max_runs"]);
    config::ea_1_o1_max_cycles = stoi(config::items["ea_1_o1_max_cycles"]);
    config::ea_1_o1_max_evals = stoi(config::items["ea_1_o1_max_evals"]);
    config::ea_1_o1_max_fitness_evals = stoi(config::items["ea_1_o1_max_fitness_evals"]);
    
    config::ea_1_o1_log_run_interval = stoi(config::items["ea_1_o1_log_run_interval"]);
    config::ea_1_o1_log_cycle_interval = stoi(config::items["ea_1_o1_log_cycle_interval"]);
    config::ea_1_o1_log_eval_interval = stoi(config::items["ea_1_o1_log_eval_interval"]);
    config::ea_1_o1_log_population_interval = stoi(config::items["ea_1_o1_log_population_interval"]);
    config::ea_1_o1_log_island_interval = stoi(config::items["ea_1_o1_log_island_interval"]);
    config::ea_1_o1_log_genome_interval = stoi(config::items["ea_1_o1_log_genome_interval"]);
    
    // algorithm 2
    
    sprintf(config::ea_2_name, "%s", config::items["ea_2_name"].c_str());
    
    config::ea_2_max_runs = stoi(config::items["ea_2_max_runs"]);
    config::ea_2_max_cycles = stoi(config::items["ea_2_max_cycles"]);
    config::ea_2_max_evals = stoi(config::items["ea_2_max_evals"]);

    config::ea_2_mu = stoi(config::items["ea_2_mu"]);
    config::ea_2_lambda = stoi(config::items["ea_2_lambda"]);
    config::ea_2_mutation_rate = stod(config::items["ea_2_mutation_rate"]);
    
    // algorithm 2, objective 1
    
    sprintf(config::ea_2_o1_name, "%s", config::items["ea_2_o1_name"].c_str());
    
    config::ea_2_o1_max_runs = stoi(config::items["ea_2_o1_max_runs"]);
    config::ea_2_o1_max_cycles = stoi(config::items["ea_2_o1_max_cycles"]);
    config::ea_2_o1_max_evals = stoi(config::items["ea_2_o1_max_evals"]);
    config::ea_2_o1_max_fitness_evals = stoi(config::items["ea_2_o1_max_fitness_evals"]);
    
    config::ea_2_log_run_interval = stoi(config::items["ea_2_log_run_interval"]);
    config::ea_2_log_cycle_interval = stoi(config::items["ea_2_log_cycle_interval"]);
    config::ea_2_log_eval_interval = stoi(config::items["ea_2_log_eval_interval"]);
    
    config::ea_2_log_population_interval = stoi(config::items["ea_2_log_population_interval"]);
    config::ea_2_log_island_interval = stoi(config::items["ea_2_log_island_interval"]);
    config::ea_2_log_genome_interval = stoi(config::items["ea_2_log_genome_interval"]);
    
    config::ea_2_o1_log_run_interval = stoi(config::items["ea_2_o1_log_run_interval"]);
    config::ea_2_o1_log_cycle_interval = stoi(config::items["ea_2_o1_log_cycle_interval"]);
    config::ea_2_o1_log_eval_interval = stoi(config::items["ea_2_o1_log_eval_interval"]);
    config::ea_2_o1_log_population_interval = stoi(config::items["ea_2_o1_log_population_interval"]);
    config::ea_2_o1_log_island_interval = stoi(config::items["ea_2_o1_log_island_interval"]);
    config::ea_2_o1_log_genome_interval = stoi(config::items["ea_2_o1_log_genome_interval"]);
    
    // parallel
    
    config::world_size = world_size;
    config::recv_cap = stoi(config::items["recv_cap"]);
    config::send_cap = stoi(config::items["send_cap"]);
    config::mu_mode = stoi(config::items["mu_mode"]);
    config::migration_interval = stoi(config::items["migration_interval"]);
    config::mu_sub = stoi(config::items["mu_sub"]);
    config::lambda_sub = stoi(config::items["lambda_sub"]);
    
//    double sub_t = stod(config::items["lambda_sub"]);
//    sub_t *= (config::mu_sub * 1.0);
//    config::lambda_sub = (int)sub_t;
    
    // parallel, mu_mode
    //
    // tuning for island population sizes
    //
    // 0: fixed global mu, relative based island mu, fixed island lambda
    // 1: proportional global mu, fixed island mu, fixed island lambda
    // 2: proportional global mu, fixed island mu, proportional island lambda
    //
    // global lambda is fairly insignificant as the current implementation
    // performs crossover at the island level, but including it for placeholder
    //
    
    sprintf(config::mu_log, "global mu: %d\r\nglobal lambda: NA", config::ea_1_mu);
    sprintf(config::sub_log,"local mu: %d\r\nlocal lambda: %d", config::mu_sub, config::lambda_sub);
    
    // TODO: generalize evolution parameters per objective
    
    if(config::mu_mode == 0) {
        // static global mu, static local lamba
        config::mu_sub = config::ea_1_mu / world_size; // relative island mu
        sprintf(config::sub_log, "local_mu: %d\r\nlocal mu calcuation: %d/%d", config::mu_sub, config::ea_1_mu, world_size);
    } else if(config::mu_mode == 1) {
        // static local mu, static local lambda
        config::ea_1_mu = world_size * config::mu_sub; // relative global mu
        sprintf(config::mu_log, "global mu: %d\r\nglobal mu calculation: %d*%d", ea_1_mu, world_size, config::mu_sub);
    } else if(config::mu_mode == 2) {
        //config::lambda_sub = config::mu_sub * stod(config::items["lambda_sub"]); // relative island lambda
        config::ea_1_mu = world_size * config::mu_sub; // relative global mu
        sprintf(config::mu_log, "global mu: %d\r\nglobal mu calculation: %d*%d", ea_1_mu, world_size, config::mu_sub);
        sprintf(config::sub_log, "local mu: %d\r\nlocal lambda: %d\r\nlocal lambda calculation: %d*%d", config::mu_sub, config::lambda_sub, config::mu_sub, config::lambda_sub);
    }
    
    char mumode[16];
    config::mu_mode < 1 ? sprintf(mumode, "%s", "fixed") : sprintf(mumode, "%s", "relative" );
    
    // program control
    
    char mode[16];
    config::ea_mode == 0 ? sprintf(mode, "%s", "benchmark") : sprintf(mode, "%s", "experiment");
    
    // variant
    
    config::dim = stoi(config::items["dim"]);
    config::sparsity = stod(config::items["sparsity"]);
    
    
    // create unique pathing for the supplied configuration so we can keep better track of results ...
    
    // create log path heirarchy to organize output by run date, communicator size, and ea mode
    
    // output heirarchy from present working directory
    // ./logs/<config[n]>/
    // ./logs/<config[n]>/<objective[n]>/
    
    char sub1[32];
    char sub2[32];

    sprintf(sub1, "logs/%d", world_size);
    sprintf(sub2, "%s/%s", sub1, mode);
    
    int fn = 1;
    bool found = true;
    while(found) {
        sprintf(config::logs_subpath, "%s/%s_%03d", sub2, config::items["config_name"].c_str(), fn);
        std::ifstream ifile;
        ifile.open(config::logs_subpath);
        if(ifile) {
            fn++;
        } else {
            found = false;
        }
    }
    
    // subdirectory for algorithm 1, objective 1 file output
    
    sprintf(config::stats_subpath_ea_1_o1, "%s/%s", config::logs_subpath, config::ea_1_o1_name);
  
    
    // subdirectory for algorithm 2, objective 1 file output
    
    sprintf(config::stats_subpath_ea_2_o1, "%s/%s", config::logs_subpath, config::ea_2_o1_name);
    

    if(world_rank == 0) {
        
        // collect current time data
        
        time_t rawtime;
        struct tm * timeinfo;

        time (&rawtime);
        timeinfo = localtime (&rawtime);
        
        char dstr[6];
        
        strftime(dstr, 80, "%m%d%y", timeinfo);
        
        // create log path heirarchy to organize output by run date, communicator size, and ea mode
        
        // output heirarchy from present working directory
        // ./logs/<config[n]>/
        // ./logs/<config[n]>/<objective[n]>/

        mkdir(sub1,0740);
        mkdir(sub2,0740);
        mkdir(config::logs_subpath, 0740);
        
        // subdirectory for algorithm 1, objective 1 file output
        
        mkdir(config::stats_subpath_ea_1_o1, 0740);
        
        // subdirectory for algorithm 2, objective 1 file output
        
        mkdir(config::stats_subpath_ea_2_o1, 0740);
        
        // program log filename
        
        sprintf(config::program_log_out, "%s/%s_%d_%s.txt", config::logs_subpath, config::items["log_file"].c_str(), world_size, mode);
        config::program_out = fopen(config::program_log_out, "w");
        
        // system log filename
        
        sprintf(config::system_log_out, "%s/%s_%d_%s.txt", config::logs_subpath, config::items["log_file"].c_str(), world_size, mode);
        config::program_out = fopen(config::program_log_out, "w");
        
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
    
        // experiment
        
        time_t stamp = time(NULL);
        fprintf(config::program_out, "configuration: %s\r\n", config::items["config_name"].c_str());
        fprintf(config::program_out, "stamp: %s\r\n", ctime(&stamp));
        fprintf(config::program_out, "seed: %u\r\n", config::seed);
        
        // titles
        
        fprintf(config::program_out, "ea[1]: %s\r\n", config::ea_1_name);
        fprintf(config::program_out, "ea[1] objective[1]: %s\r\n", config::ea_1_o1_name);
        fprintf(config::program_out, "ea[2]: %s\r\n", config::ea_2_name);
        fprintf(config::program_out, "ea[2] objective[1]: %s\r\n", config::ea_2_o1_name);
        
        // file locations
        
        fprintf(config::program_out, "log file: %s\r\n", config::program_log_out);
        fprintf(config::program_out, "world size: %d\r\n", world_size);
        
        // parallell
        
        fprintf(config::program_out, "mu_mode_name: %s\r\n", mode);
        fprintf(config::program_out, "mu_mode: %d ", config::mu_mode);
        fprintf(config::program_out, "%s\r\n", config::mu_log);
        fprintf(config::program_out, "%s\r\n", config::sub_log);
        
        // log all parameter values as read from file into items collection
        
        std::map <std::string, std::string>::iterator items_iterator;
        
        for(items_iterator = std::next(items.begin()); items_iterator != items.end(); ++items_iterator) {
            fprintf(config::program_out, "%s: %s\r\n", items_iterator->first.c_str(), items_iterator->second.c_str());
        }
        
    }
    
    // build filenames from parameter values and prepare for output
    
    // TODO: create or use an appropriate data type to house each
    // ea and corresponding objective parameters to allow for more
    // automated specification and improve the code efficiency
    
    // as of now, the code readily accomodates two eas, each with one
    // objective (hardcoded), additional instances require source modification
    
    
    if(world_rank == 0) {
        
        if(config::ea_1_log_run_interval != 0) {
            sprintf(config::ea_1_stats_run, "%s/runs_%d.csv", config::stats_subpath_ea_1_o1, world_size);
            config::ea_1_run_out = fopen(config::ea_1_stats_run, "w");
        }
        
        if(config::ea_1_log_cycle_interval != 0) {
            sprintf(config::ea_1_stats_cycle, "%s/cycles_%d.csv", config::stats_subpath_ea_1_o1, world_size);
            config::ea_1_cycle_out = fopen(config::ea_1_stats_cycle, "w");
        }
        
        if(config::ea_1_log_eval_interval != 0) {
            sprintf(config::ea_1_stats_eval, "%s/evals_%d.csv", config::stats_subpath_ea_1_o1, world_size);
            config::ea_1_eval_out = fopen(config::ea_1_stats_eval, "w");
        }
        
        if(config::ea_1_log_population_interval != 0) {
            sprintf(config::ea_1_stats_population, "%s/evo_%d.csv", config::stats_subpath_ea_1_o1, world_size);
            config::ea_1_population_out = fopen(config::ea_1_stats_population, "w");
        }
        
        if(config::ea_1_log_genome_interval != 0) {
            sprintf(config::ea_1_stats_genome, "%s/genome_%d.csv", config::stats_subpath_ea_1_o1, world_size);
            config::ea_1_genome_out = fopen(config::ea_1_stats_genome, "w");
        }
            
        if(config::ea_2_log_run_interval != 0) {
            sprintf(config::ea_2_stats_run, "%s/runs_%d.csv", config::stats_subpath_ea_2_o1, world_size);
            config::ea_2_run_out = fopen(config::ea_2_stats_run, "w");
        }
        
        if(config::ea_2_log_cycle_interval != 0) {
            sprintf(config::ea_2_stats_cycle, "%s/cycles_%d.csv", config::stats_subpath_ea_2_o1, world_size);
            config::ea_2_cycle_out = fopen(config::ea_2_stats_cycle, "w");
        }
        
        if(config::ea_2_log_eval_interval != 0) {
            sprintf(config::ea_2_stats_eval, "%s/evals_%d.csv", config::stats_subpath_ea_2_o1, world_size);
            config::ea_2_eval_out = fopen(config::ea_2_stats_eval, "w");
        }
        
        if(config::ea_2_log_population_interval != 0) {
            sprintf(config::ea_2_stats_population, "%s/evo_%d.csv", config::stats_subpath_ea_2_o1,  world_size);
            config::ea_2_population_out = fopen(config::ea_2_stats_population, "w");
        }
        
        if(config::ea_2_log_genome_interval != 0) {
            sprintf(config::ea_2_stats_genome, "%s/genome_%d.csv", config::stats_subpath_ea_2_o1, world_size);
            config::ea_2_genome_out = fopen(config::ea_2_stats_genome, "w");
        }
        
        if (ea_1_log_run_interval != 0) {
            fprintf(config::ea_1_run_out, "run,cycle,eval,average_fitness,best_fitness,avg_best_fitness,min_scatter_t,min_scatter_i,max_scatter_t,max_scatter_i,sum_scatter_t,avg_scatter_t,min_gather_t,min_gather_i,max_gather_t,max_gather_i,sum_gather_t,avg_gather_t,min_migration_t,min_migration_i,max_migration_t,max_migration_i,sum_migration_t,avg_migration_t,min_run_t,min_run_i,max_run_t,max_run_i,sum_run_t,avg_run_t,start,duration\r\n");
        }
        
        if(ea_1_log_cycle_interval != 0) {
            fprintf(config::ea_1_cycle_out, "run,cycle,eval,average_fitness,best_fitness,avg_best_fitness,min_scatter_t,min_scatter_i,max_scatter_t,max_scatter_i,sum_scatter_t,avg_scatter_t,min_gather_t,min_gather_i,max_gather_t,max_gather_i,sum_gather_t,avg_gather_t,min_migration_t,min_migration_i,max_migration_t,max_migration_i,sum_migration_t,avg_migration_t,min_cycle_t,min_cycle_i,max_cycle_t,max_cycle_i,sum_cycle_t,avg_cycle_t,start,duration\r\n");
        }
        
        if(config::ea_1_log_eval_interval != 0) {
            fprintf(config::ea_1_eval_out, "run,cycle,eval,average_fitness,best_fitness,avg_best_fitness,min_scatter_t,min_scatter_i,max_scatter_t,max_scatter_i,sum_scatter_t,avg_scatter_t,min_gather_t,min_gather_i,max_gather_t,max_gather_i,sum_gather_t,avg_gather_t,min_migration_t,min_migration_i,max_migration_t,max_migration_i,sum_migration_t,avg_migration_t,min_eval_t,min_eval_i,max_eval_t,max_eval_i,sum_eval_t,avg_eval_t,start,duration\r\n");
        }
        
        // TODO: add stats: global_min_fitness,global_min_genome_id,
        
        if(config::ea_1_log_population_interval != 0) {
            //fprintf(config::ea_1_population_out, "run,cycle,eval,average_fitness,global_best_fitness,global_best_genome_id,average_global_best_fitness,ea_t,run_t,cycle_t,eval_t,min_scatter_t,min_scatter_i,max_scatter_t,max_scatter_i,sum_scatter_t,avg_scatter_t,min_gather_t,min_gather_i,max_gather_t,max_gather_i,sum_gather_t,avg_gather_t,min_migration_t,min_migration_i,max_migration_t,max_migration_i,sum_migration_t,avg_migration_t,min_eval_t,min_eval_i,max_eval_t,max_eval_i,sum_eval_t,avg_eval_t\r\n");
            fprintf(config::ea_1_population_out, "run,cycle,eval,genome_id,fitness,source,locale,parent1,parent2,selected,survived\r\n");
                         
        }
        
        //fprintf(config::ea_1_genome_out, "run,cycle,eval,id,rank,fitness,n_best,n_min,n_selected,n_survived,distribution,birth_cycle,parents,genes,meta_genome_id,origin,locale,n_migrations,visas,\r\n");
        
        if(config::ea_1_log_genome_interval != 0) {
            fprintf(config::ea_1_genome_out, "interval, interval_id, genome_id, genome_genes\r\n");
        }
        
        if(config::ea_2_log_run_interval != 0) {
            fprintf(config::ea_2_run_out, "run,cycle,eval,average_fitness,best_fitness,avg_best_fitness,min_scatter_t,min_scatter_i,max_scatter_t,max_scatter_i,sum_scatter_t,avg_scatter_t,min_gather_t,min_gather_i,max_gather_t,max_gather_i,sum_gather_t,avg_gather_t,min_migration_t,min_migration_i,max_migration_t,max_migration_i,sum_migration_t,avg_migration_t,min_run_t,min_run_i,max_run_t,max_run_i,sum_run_t,avg_run_t,start,duration\r\n");
        }
        
        if(config::ea_2_log_cycle_interval != 0) {
            fprintf(config::ea_2_cycle_out, "run,cycle,eval,average_fitness,best_fitness,avg_best_fitness,min_scatter_t,min_scatter_i,max_scatter_t,max_scatter_i,sum_scatter_t,avg_scatter_t,min_gather_t,min_gather_i,max_gather_t,max_gather_i,sum_gather_t,avg_gather_t,min_migration_t,min_migration_i,max_migration_t,max_migration_i,sum_migration_t,avg_migration_t,min_cycle_t,min_cycle_i,max_cycle_t,max_cycle_i,sum_cycle_t,avg_cycle_t,start,duration\r\n");
        }
        
        if(config::ea_2_log_eval_interval != 0) {
            fprintf(config::ea_2_eval_out, "run,cycle,eval,average_fitness,best_fitness,avg_best_fitness,min_scatter_t,min_scatter_i,max_scatter_t,max_scatter_i,sum_scatter_t,avg_scatter_t,min_gather_t,min_gather_i,max_gather_t,max_gather_i,sum_gather_t,avg_gather_t,min_migration_t,min_migration_i,max_migration_t,max_migration_i,sum_migration_t,avg_migration_t,min_eval_t,min_eval_i,max_eval_t,max_eval_i,sum_eval_t,avg_eval_t,start,duration\r\n");
        }
        
        if(config::ea_2_log_population_interval != 0) {
            fprintf(config::ea_2_population_out, "run,cycle,eval,average_fitness,best_fitness,best_id,avg_best_fitness,ea_t,run_t,cycle_t,eval_t,min_scatter_t,min_scatter_i,max_scatter_t,max_scatter_i,sum_scatter_t,avg_scatter_t,min_gather_t,min_gather_i,max_gather_t,max_gather_i,sum_gather_t,avg_gather_t,min_migration_t,min_migration_i,max_migration_t,max_migration_i,sum_migration_t,avg_migration_t,min_eval_t,min_eval_i,max_eval_t,max_eval_i,sum_eval_t,avg_eval_t\r\n");
        }
        
        if(config::ea_2_log_genome_interval != 0) {
            fprintf(config::ea_2_genome_out, "run,cycle,eval,id,rank,fitness,n_best,n_min,n_selected,n_survived,distribution,birth_cycle,parents,genes,meta_genome_id,origin,locale,n_migrations,visas,\r\n");
        }
        
        // TODO: add per-island logs with some interval with stats e.g., island_id, island_xxx_min|max|sum|avg_t, island_avg_fit, island_best_fit, island_min_fit, island_channels, island_total_migrations, island_hostname, island_nprocs
        
        fprintf(config::program_out, "ea1 run log: %s\r\n", config::ea_1_stats_run);
        fprintf(config::program_out, "ea1 cycle log: %s\r\n", config::ea_1_stats_cycle);
        fprintf(config::program_out, "ea1 eval log: %s\r\n", config::ea_1_stats_eval);
        fprintf(config::program_out, "ea1 population log: %s\r\n", config::ea_1_stats_population);
        fprintf(config::program_out, "ea1 genome log: %s\r\n", config::ea_1_stats_genome);

        fprintf(config::program_out, "ea2 run log: %s\r\n", config::ea_2_stats_run);
        fprintf(config::program_out, "ea2 cycle log: %s\r\n", config::ea_2_stats_cycle);
        fprintf(config::program_out, "ea2 eval log: %s\r\n", config::ea_2_stats_eval);
        fprintf(config::program_out, "ea2 population log: %s\r\n", config::ea_2_stats_population);
        fprintf(config::program_out, "ea2 genome log: %s\r\n", config::ea_2_stats_genome);

        // fflush(config::ea_out);
        // fflush(config::system_out);
        fflush(config::program_out);
        
        fflush(config::ea_1_run_out);
        fflush(config::ea_1_cycle_out);
        fflush(config::ea_1_eval_out);
        fflush(config::ea_1_population_out);
        fflush(config::ea_1_genome_out);
        
        fflush(config::ea_out);
        fflush(config::ea_2_run_out);
        fflush(config::ea_2_cycle_out);
        fflush(config::ea_2_eval_out);
        fflush(config::ea_2_population_out);
        fflush(config::ea_2_genome_out);
        
    }
    
    // associate each mpi rank (island) with
    
    config::id_field_prefix1 = std::make_pair(world_rank, 3);
    
}

#endif /* config_h */

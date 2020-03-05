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

namespace config {

    std::map<std::string, std::string> items;

    FILE *log_out;
    FILE *stats_out;

    int evals = 0;
    int runs = 0;
    int lambda = 0;
    int mu = 0;
    int seed = 0;

    double mutation_rate = 0.0;

    char log_fname[100];
    char stats_fname[100];
    char logs_subpath[100];
    char stats_subpath[100];

    void load(const char *input, int world_size);

};

void config::load(const char *input, int world_size) {
    
    if(input == NULL) { input = (char *)"config.txt"; }
    
    std::ifstream config_file(input);

    std::string key;
    std::string value;
    
    while (config_file.good()) {
       
       getline(config_file, key, ':');
       getline(config_file, value, '\n');
       
       if(key != "") {
           printf("%s => %s\r\n", key.c_str(), value.c_str());
           config::items[key] = value;
       }
       
    }
    
    config::runs = stoi(config::items["runs"]);
    config::evals = stoi(config::items["evals"]);
    config::lambda = stoi(config::items["lambda"]);
    config::mu = stoi(config::items["mu"]);
    config::mutation_rate = stod(config::items["mutation_rate"]);
    
    // create unique pathing for the supplied configuration so we can keep better track of results ...
    
    sprintf(config::logs_subpath, "logs/%s_%ld", config::items["config_name"].c_str(), time(0));
    mkdir(config::logs_subpath, 0740);
    sprintf(config::stats_subpath, "stats/%s_%ld", config::items["config_name"].c_str(), time(0));
    mkdir(config::stats_subpath, 0740);
    
    // create and initialize the log file with the parsed configuration values ...
    
    sprintf(config::log_fname, "%s/%s_%d_%ld.txt", config::logs_subpath, config::items["log_file"].c_str(), world_size, time(0));
    config::log_out = fopen(config::log_fname, "w");
    
    // seed the rng ...
    
    if(stoi(config::items["seed"]) == 0) {
        timeval time;
        gettimeofday(&time, NULL);
        config::seed = (unsigned int)(time.tv_sec * 1000) + time.tv_usec;
    } else {
        config::seed = stoi(config::items["seed"]);
    }
                            
    srand(config::seed);
    
    fprintf(config::log_out, "seed: %u\r\n", config::seed);
    
    std::map <std::string, std::string>::iterator items_iterator;
    
    for(items_iterator = std::next(items.begin()); items_iterator != items.end(); ++items_iterator) {
        fprintf(config::log_out, "%s: %s\r\n", items_iterator->first.c_str(), items_iterator->second.c_str());
    }
    
    fprintf(config::log_out, "log file: %s\r\n", config::log_fname);
    
    fflush(config::log_out);
    
    sprintf(config::stats_fname, "%s/%s_%d_%ld.txt", config::stats_subpath, config::items["stats_file"].c_str(), world_size, time(0));
    config::stats_out = fopen(config::stats_fname, "w");
    
    fprintf(config::stats_out, "run\teval\taverage_global_fitness\taverage_best_fitness\taverage_island_comm_time\r\n");
    
}

const static unsigned int DIM = 100;

#endif /* config_h */

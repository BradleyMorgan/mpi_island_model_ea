//
//  dtype_heirarchy.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 7/22/21.
//  Copyright © 2021 Bradley Morgan. All rights reserved.
//

#ifndef dtype_heirarchy_h
#define dtype_heirarchy_h


char* uniqid(unsigned long long int instance) {
    //char log[256];
    unsigned long long int count = instance<=0 ? 1 : instance;
    //printf("field prefix= [%d,%d] ", config::id_field_prefix1.first, config::id_field_prefix1.second);
    int prefix = config::id_field_prefix1.second <= 0 ? 3 : config::id_field_prefix1.second;
    //printf("prefix=%d ", prefix);
    int rank=config::id_field_prefix1.first+1;
    //printf("rank=%d ", rank);
    int fixed_len=prefix;
    //if(fixed_len != 3) { printf("hmm %d\r\n", fixed_len);  }
    //printf("fixed_len=%d ", fixed_len);
    //int rank_len=ceil(log10(rank));
    //printf("rank_len=%d ", rank_len);
    int instance_len=ceil(log10(count));
    //printf("instance_len=%d\r\n", instance_len);
    static char buffer[64];
    sprintf(buffer, "%d%0*d%0*llu", rank, fixed_len, 0, instance_len, count);
    //printf("1%0*d%0*d\r\n", fixed_len, 0, instance_len, 0);
    //unsigned long long int id_base = std::stoll(buffer);
    //unsigned long long int unique_id = (rank * id_base) + instance;
    //printf("base=%llu result=%llu\r\n", id_base, unique_id);
    //printf("%s", log);
    //if((unique_id * .00000001) > 17) { printf("hmm %llu\r\n", unique_id); }
    return buffer;
}


#endif /* dtype_heirarchy_h */

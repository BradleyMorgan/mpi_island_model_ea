//
//  dtype_topology.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 7/22/21.
//  Copyright © 2021 Bradley Morgan. All rights reserved.
//

#ifndef dtype_topology_h
#define dtype_topology_h

//
// the @comm_group{} datatype is used primarily to track the required MPI_Send() and MPI_Recv()
// operations and the required parameter values needed to perform the island (process)
// execution context, e.g. from the perspective of island @comm_group{island_id}.
//
//  i   1    2    3    4
//  1 [ 0 ][ 1 ][ 0 ][ 0 ] => island 1 sends to island 2
//  2 [ 1 ][ 0 ][ 1 ][ 0 ] => island 2 sends to island 1 and 3
//  3 [ 0 ][ 0 ][ 0 ][ 0 ] => island 3 no channels
//  4 [ 1 ][ 0 ][ 0 ][ 0 ] => island 4 sends to island 1
//
// +----+-----------------------------------------+--------------------------------------+
// | id |        assigned send channels           |     assigned receive channels        |
// +----+-----------------------------------------+--------------------------------------+
// | 1  |  @comm_group[id=1]{receivers}[0] = 2    |  @comm_group[id=2]{senders}[0] = 1   |
// +----+-----------------------------------------+--------------------------------------+
// | 2  |  @comm_group[id=2]{receivers}[0] = 1    |  @comm_group[id=1]{senders}[0] = 2   |
// |    |  @comm_group[id=2]{receivers}[1] = 3    |  @comm_group[id=3]{senders}[0] = 2   |
// +----+-----------------------------------------+--------------------------------------+
// | 3  |  @comm_group[id=3]{receivers][0] = NIL  |  @comm_group[id=?]{} NOOP            |
// +----+-----------------------------------------+--------------------------------------+
// | 4  |  @comm_group[id=4]{receivers}[0] = 1    |  @comm_group[id=1]{senders}[1] = 4   |
// +----+-----------------------------------------+--------------------------------------+
//

struct channel {
    
    int id;
    int count;
    
    // track fitness in terms of island communcation operations (time)
    
    double fitness = 0.0;
    
    std::vector<int> senders;
    std::vector<int> receivers;
    
};

#pragma mark DATATYPE: topology{}

// represents the full island topology calculated from adjacency                            |
// matrices and mapped to into the @channel{} datatype.                                     |

struct topology {
  
    int id = 0;
    int rounds = 0;
    int world_size;
    int channel_count = 0;
    
    double fitness = 0.0; // track in terms of aggregate communication time
    double round_fitness = 0.0; // number of evaluations performed
    double total_migration_time = 0.0;
    double selection_distribution = 0.0; // fitness measurement for selection method
    
    // an array of @comm_group{} describing the full context MPI_Send() and MPI_Recv()
    // operations forms the mpi-suitable topology representation
    
    std::vector<channel> channels;
    
    struct create {
     
        static void dynamic(topology &t);
        static void channels(topology &t, std::vector<std::vector<int>> &matrix);
        
        static std::vector<std::vector<int>> dynamic_matrix(const int world_size);
        static std::vector<std::vector<int>> matrix(const topology &t);
        
    };
    
    void validate(island &isle);
    void distribute(island &isle);
    void apply(island &isle, topology &t);
    //void evaluate(island &isle);
    
    topology() {}
    //printf("[b=%llu|d=%llu] DEFAULT CONSTRUCTING topology<%s> [%p]\r\n", this->id, this->id, typeid(this).name(), this);
    
};


#endif /* topology_h */

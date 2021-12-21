//
//  dtype_topology.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 7/22/21.
//  Copyright © 2021 Bradley Morgan. All rights reserved.
//

#ifndef dtype_topology_h
#define dtype_topology_h

// this datatype uses an adjacency matrix (2d array) to represent and measure a communication topology

// the design is likely general enough to use as a representation for any experimentation that
// requires a way to manipulate or dynamically generate communication patterns, however the
// intent is to find optimal communication patterns through applying the represented topology
// to a communication-heavy algorithm (e.g. an island-model ea), with optimality measured in
// terms of time incurred during data transfer.

#pragma mark DATATYPE: @channel{}

// the @channel{} datatype is used primarily to track the required MPI_Send() and MPI_Recv()
// operations and the required parameter values needed to perform the island (process)
// execution context, e.g. from the perspective of island @channel{island_id}.
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
// | 1  |  @channel[id=1]{receivers}[0] = 2       |  @channel[id=2]{senders}[0] = 1      |
// +----+-----------------------------------------+--------------------------------------+
// | 2  |  @channel[id=2]{receivers}[0] = 1       |  @channel[id=1]{senders}[0] = 2      |
// |    |  @channel[id=2]{receivers}[1] = 3       |  @channel[id=3]{senders}[0] = 2      |
// +----+-----------------------------------------+--------------------------------------+
// | 3  |  @channel[id=3]{receivers][0] = NIL     |  @channel[id=?]{} NOOP               |
// +----+-----------------------------------------+--------------------------------------+
// | 4  |  @channel[id=4]{receivers}[0] = 1       |  @channel[id=1]{senders}[1] = 4      |
// +----+-----------------------------------------+--------------------------------------+
//
// this is convenient for protocols that require explicit, symmetric endpoint communication,
// where the initialization of a matching source sender and target receiver is required
// before data transfer occurs, e.g. MPI
//

struct channel {
    
    int id;
    int count;
    
    // track fitness in terms of island communcation operations (time)
    
    double fitness = 0.0;
    
    std::vector<int> senders;
    std::vector<int> receivers;
    
};

struct topology_stats : time_stats, parallel_stats, fitness_stats {
    
    void minmax(mpi_local topology_stats::*field, mpi_local result, double const &(*func)(double const&, double const&)) {
        this->*field = (this->*field).value == 0.0 ? result : func((this->*field).value, result.value) != (this->*field).value ? result : this->*field;
    }
    
};

#pragma mark EA::META::OBJECTIVE::DATATYPE: @topology{}

// represents a communication topology, calculated from adjacency
// matrices and mapped to into the @channel{} datatype.

struct topology {
  
    int id = 1;
    
    int evaluations = 0;
    
    int world_size = 0;
    int channel_count = 0;
    
    double fitness = 0.0;
    
    topology_stats stats;
    
    double selection_distribution = 0.0;
    
    // an array of @channel{} describing the full context MPI_Send() and MPI_Recv()
    // operations forms the mpi-suitable topology representation
    
    std::vector<channel> channels;
    
    // encapsulation for initialization methods
    
    struct create {
     
        static void dynamic(topology &t);
        static void channels(topology &t, std::vector<std::vector<int>> &matrix);
        
        static std::vector<std::vector<int>> dynamic_matrix(const int world_size);
        static std::vector<std::vector<int>> matrix(const topology &t);
        
    };
    
    // primary methods, see topology.h
    
    void validate(island &isle); // experimental method for sanity-checking the representation
    void distribute(island &isle); // from a single source (root), initate the distribution of channels to each communicating node
    void apply(island &isle, topology &t); // for all communicating nodes, receive the corresponding channels as defined in the topology
    
    template<typename i> void measure(i &interval);
    template<typename o> void minmax(mpi_local topology_stats::*field, mpi_local result, o *op);
    
    int numeric_id();
    
    topology(): world_size(config::world_size) {};
    
};

int topology::numeric_id() { return this->id; }


template<typename i> void topology::measure(i &interval) {
    
    if(mpi.id != 0) { return; }
    this->fitness = (interval.stats.sum_t.value / interval.max) * -1;
    
    this->stats.minmax(&topology_stats::min_t, interval.stats.min_t, std::min<double>);
    this->stats.minmax(&topology_stats::max_t, interval.stats.max_t, std::max<double>);
    this->stats.sum_t.value += interval.stats.sum_t.value;
    this->stats.avg_t = this->stats.sum_t.value / mpi.size;
    
    this->stats.minmax(&topology_stats::min_scatter_t, interval.stats.min_scatter_t, std::min<double>);
    this->stats.minmax(&topology_stats::max_scatter_t, interval.stats.max_scatter_t, std::max<double>);
    this->stats.sum_scatter_t.value += interval.stats.sum_scatter_t.value;
    this->stats.avg_scatter_t = this->stats.sum_scatter_t.value / mpi.size;
    
    this->stats.minmax(&topology_stats::min_gather_t, interval.stats.min_gather_t, std::min<double>);
    this->stats.minmax(&topology_stats::max_gather_t, interval.stats.max_gather_t, std::max<double>);
    this->stats.sum_gather_t.value += interval.stats.sum_gather_t.value;
    this->stats.avg_gather_t = this->stats.sum_gather_t.value / mpi.size;
    
    this->stats.minmax(&topology_stats::min_migration_t, interval.stats.min_migration_t, std::min<double>);
    this->stats.minmax(&topology_stats::max_migration_t, interval.stats.max_migration_t, std::max<double>);
    this->stats.sum_migration_t.value += interval.stats.sum_migration_t.value;
    this->stats.avg_migration_t = this->stats.sum_migration_t.value / mpi.size;
    
}

#endif /* topology_h */

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
    
    bool log_head = false;
    bool log_tail = false;
    
    int head_interval = 1;
    int tail_interval = 1;
 
    double total_o1_fitness = 0.0;
    double total_o2_fitness = 0.0;
    
    double avg_o1_fitness = 0.0;
    double avg_o2_fitness = 0.0;
    
    double avg_distance;
    
    int send_channels = 0;
    int recv_channels = 0;
    int total_channels = 0;
    
    int target_runs = 0;
    int migrations = 0;
    int departures = 0;
    int arrivals = 0;
    
    void minmax(mpi_local topology_stats::*field, mpi_local result, double const &(*func)(double const&, double const&)) {
        this->*field = (this->*field).value == 0.0 ? result : func((this->*field).value, result.value) != (this->*field).value ? result : this->*field;
    }
    
    topology_stats(bool head, bool tail) : log_head(head), log_tail(tail) {};
    topology_stats() : log_head(false), log_tail(false) {};
    
};

#pragma mark EA::META::OBJECTIVE::DATATYPE: @topology{}

// represents a communication topology, calculated from adjacency
// matrices and mapped to into the @channel{} datatype.

struct topology {
  
    int id = 1;
    
    topology_stats stats;
    
    double fitness = 0.0;
    double selection_distribution = 0.0;
    
    // multi-objective fitness
    
    double distance = 0.0;
    
    std::pair<double, double> fitness_multi = { 0.0, 0.0 };
    
    // domination
    
    int dom_rank = 0;
    int dom_count = 0;
    
    std::vector<topology*> dom_genomes;

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
    
    template<typename i> void log(i &interval);
    template<typename i> void measure(i &interval);
    template<typename o> void minmax(mpi_local topology_stats::*field, mpi_local result, o *op);
    
    int numeric_id();
    
    bool dominates(topology &t);
    
    void init(int &genome_id);
    void init_multi();
    
    topology() {};
    
};

bool topology::dominates(topology &t) {
    
    // the equalities used here are for a maximization problem
    
    // to dominate, an individual must be at least as good as the individual
    // to which it is compared for all objectives ...
    // and strictly better for one (or more) of the solution objectives
    //
    // so for maximization, if any of the objective fitness values in this
    // solution (this) are less than (<) that of the individual (t) to which
    // we are comparing, we can disqualify it as dominate.
    //
    // if that criteria is met, this solution (this) is dominate if it is better (>)
    // than the provided solution (t) for any objective
    //
    
    if(this->fitness_multi.first < t.fitness_multi.first || this->fitness_multi.second < t.fitness_multi.second) { return false; }
    if(this->fitness_multi.first > t.fitness_multi.first || this->fitness_multi.second > t.fitness_multi.second) { return true; }
    
    return false;
    
}

void topology::init(int &genome_id) {
    
    this->id = genome_id;
    this->fitness = 0.0;
    this->stats = {};
    this->selection_distribution = 0.0;
    this->channels.resize(mpi.size);
    this->channels.clear();
    
}

void topology::init_multi() {
    
    if(mpi.id == 0) {
        
        this->dom_genomes.clear();
        this->dom_rank = 0;
        this->dom_count = 0;
        this->distance = 0.0;
        
    }
    
}

int topology::numeric_id() { return this->id; }

template<typename i> void topology::measure(i &interval) {
    
    if(mpi.id != 0) { return; }
    
    this->fitness += (interval.stats.sum_t.value * -1);
    this->fitness_multi.first = this->fitness;
    
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

template<typename i> void topology::log(i &interval) {
    
    //if(interval.log_fout == true) {

        if(mpi.id != 0) { return; }
    
        fprintf(config::ea_2_genome_out, "%s,%d,%d,%f,%f,%d,%lu,%d,%d,%d,%d,%d,%d,%f\r\n",

            interval.name,
            interval.id,
            this->id,
            this->fitness_multi.first,
            this->fitness_multi.second,
            this->dom_rank,
            this->dom_genomes.size(),
            this->stats.send_channels,
            this->stats.recv_channels,
            this->stats.total_channels,
            this->stats.departures,
            this->stats.arrivals,
            this->stats.migrations,
            this->selection_distribution);

    //}
    
    fflush(config::ea_2_genome_out);

}

#endif /* topology_h */

//
//  topology.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 4/12/21.
//  Copyright © 2021 Bradley Morgan. All rights reserved.
//

#ifndef topology_h
#define topology_h

#pragma mark DATATYPE: channel{}

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
    double selection_distribution; // fitness measurement for selection method
    
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
    void evaluate(island &isle);
            
};

#pragma mark FUNCTION: create::channels()

// from a provided adjacency matrix, create sender and receiver arrays for                  |
// each island for easier use with MPI send and receive ...                                 |

void topology::create::channels(topology &t, std::vector<std::vector<int>> &matrix) {
    
    LOG(10, 0, 0, "creating group...\r\n");
    
    // indices in this array of @comm_group{} map directly to a specific island (process)
    
    t.channels.resize(matrix.size());
    
    // iterate through the adjacency matrix and assign any marked channels
    // to corresponding island sender or receiver array ...
    
    // matrix row, y-axis, island send to marked columns
    
    for(int i=0; i<matrix.size(); i++) {
        
        // matrix row, y-axis, if 1 island i send to column -> i [0][1][0][0]
        
        for(int j=0; j<matrix[i].size(); j++) {
        
            // matrix column, x-axis, if 1 island receives from current row -> i
            
            LOG(10, 0, 0, "%d ", matrix[i][j]);
            
            if(matrix[i][j] == 1) { // i sends to j ...
                
                LOG(5, 0, 0, "adding receiver %d to %d and sender %d to %d\r\n", j, i, i, j);
                
                t.channels[i].receivers.push_back(j); // add island j to island i receivers list
                t.channels[j].senders.push_back(i); // add island i to island j senders list
                
                t.channels[i].count++;
                t.channels[j].count++;
                
                t.channel_count += 2;
                
            }
            
        }
        
        LOG(10, 0, 0, "\r\n");
        
    }
    
    // the result is an mpi-suitable representation of a full communication topology
    
    LOG(6, 0, 0, "returning communicator size %lu\r\n", t.channels.size());
    
}

#pragma mark FUNCTION: create::matrix()

// with a provided vector of topology groups (which holds an island's senders and receivers),
// create an adjacency matrix of size world_size*world_size, populating the coordinate
// index with 1 whenever a receiver is found ...
//
//    1  2  3  4
// 1 [0][1][0][0]  ->  1 sends to 2
// 2 [0][0][1][0]  ->  2 sends to 3
// 3 [0][0][0][0]  ->  3 sends to nobody
// 4 [1][0][0][0]  ->  4 sends to 1

std::vector<std::vector<int>> topology::create::matrix(const topology &t) {
    
    // the passed vector comm should contain world_size group objects,
    // tracking the senders and receivers for each island
    
    std::vector<std::vector<int>> matrix;
    matrix.resize(t.world_size);

    for(int i=0; i<t.world_size; i++) { // matrix row, sized with world_size indices

        matrix[i].resize(t.world_size);
        
        for(int j=0; j<t.world_size; j++) { // matrix column [0..world size]
            
            // here we only mark the receiving nodes for each island row with a 1
            // a boolean value of true within a a row's index will tell us what we need to know
            // to track the communication channel on both sides
            
            if(std::find(t.channels[i].receivers.begin(), t.channels[i].receivers.end(), j) != t.channels[i].receivers.end()) {
                // if island i's receivers array contains island j, mark that as a sender\receiver channel ...
                matrix[i][j] = 1;
            } else {
                // otherwise, no sender\receiver channel exists, so 0 ...
                matrix[i][j] = 0;
            }
            
        }
        
    }
    
    return matrix;
    
}

#pragma mark FUNCTION: create::dynamic_matrix()

// create a randomly populated adjacency matrix of size world_size*world_size,
// communication neighbors with probability determined by the sparsity parameter ...

std::vector<std::vector<int>> topology::create::dynamic_matrix(const int world_size) {
    
    std::vector<std::vector<int>> matrix;
    matrix.resize(world_size);
    
    int comm_count = 0;
   
    int rec_count[world_size];
    int snd_count[world_size];
    
    // intialize empty matrix ...
    
    for(int i=0; i<world_size; i++) {
        rec_count[i] = 0;
        snd_count[i] = 0;
    }
    
    while(comm_count == 0) {  // avoid empty matrix if sparsity probability is low
    
        for(int i=0; i<world_size; i++) { // row

            matrix[i].resize(world_size);
            
            for(int j=0; j<world_size; j++) { // column
                
                // we don't want an island sending migrants to itself
                
                if(i == j) { matrix[i][j] = 0; continue; }
                    
                // check @config.txt[dim] max send\recv ...
                
                if(rec_count[j] >= config::migration_cap) {
                    LOG(6, 0, 0, "receive cap limit reached for process %d\r\n", i);
                    continue;
                }
                
                if(snd_count[i] >= config::send_cap) {
                    LOG(6, 0, 0, "send cap limit reached for process %d\r\n", i);
                    continue;
                }
                
                // check @config.txt[sparsity] probability ...
                
                if(prob_true(config::sparsity) && i != j) {
                    matrix[i][j] = 1;
                    comm_count++;
                    rec_count[j]++;
                    snd_count[i]++;
                } else {
                    matrix[i][j] = 0;
                }
                
            }
            
        }
        
    }
    
    // constrained, randomly generated 0s and 1s in a @config.txt[dim] sized logical grid ...
    
    return matrix;
    
}

//#pragma mark FUNCTION: topo_cpd()
//
//// calculate the cumulative probability distribution for the topology population ...
//
//void topologies::calculate::cpd(topologies &t) {
//
//    LOG(10, 0, 0, "generating cpd\r\n");
//
//    t.cpd.clear();
//
//    double total_fitness = 0.0;
//    double cumulative_probability = 0.0;
//
//    for(int i=0; i<t.population.size(); i++) {
//
//        total_fitness += t.population[i].fitness;
//
//    }
//
//    std::reverse(t.population.begin(), t.population.end());
//
//    for(int i=0; i<t.population.size(); i++) {
//
//        t.population[i].selection_distribution = (double)t.population[i].fitness / t.total_fitness;
//
//        cumulative_probability += t.population[i].selection_distribution;
//        t.cpd.push_back(cumulative_probability);
//
//    }
//
//}
//
//#pragma mark FUNCTION: select_topo_parent()
//
//// parent selection from the topology population, based on the cumulative probability distribution ...
//
//topology select_topo_parent(std::vector<topology> &p) {
//
//    LOG(7, 0, 0, "selecting topo parent: ");
//
//    topology t;
//
//    int i = 1;
//
//    double r = ((double)rand()/(double)RAND_MAX);
//
//    while (p.cpd[i] < r ) { i++; }
//
//    t = p[i];
//
//    LOG(7, 0, 0, "%d\r\n", i);
//
//    return t;
//
//}

#pragma mark FUNCTION: create_dyn_topology()

// used to create the initial topology population ...

void topology::create::dynamic(topology &t) {
    
    LOG(5, 0, 0, "generating dynamic topology %d world size %d\r\n", t.id, t.world_size);
    std::vector<std::vector<int>> matrix = topology::create::dynamic_matrix(t.world_size);
    
    topology::create::channels(t, matrix);
    
}

#pragma mark FUNCTION: topology_apply()

// issue reciprocal MPI_Recv calls for each matching MPI_Send from the root process                     |
// (see @topology_distribute()) as determined fom the provided @island{} values.                        |
//                                                                                                      |
// on success, all ranks will have received a copy of their send\receive operation metadata             |
//                                                                                                      |
// should be called from all islands                                                                    |
//                                                                                                      |

void topology::apply(island &isle, topology &t) {
    
    // rank authoritative source parses and queues the distribution of
    // island-specific migration channels to be assigned
    
    LOG(4, isle.id, 0, "applying topology %d ...\r\n", t.id);
    
    if(isle.id == 0) {
        
        LOG(4, 0, 0, "distributing topology %d\r\n", t.id);
    
        t.distribute(isle);
    
    }
    
    if(t.world_size == 0) {
        MPI_Comm_size(MPI_COMM_WORLD, &t.world_size);
    }
    
    // initialize *this* (current MPI rank) island's send and receive queue sizes, assume it's empty ...
    
    int send_size = 0;
    int rec_size = 0;
    
    // initiate MPI_Recv calls.
    
    LOG(4, 0, 0, "rank %d receiving topology %d send\\receive metadata world size = %d\r\n", isle.id, t.id, t.world_size);
    
    // complete the MPI_Send of queue size for island export
    
    MPI_Recv(&send_size, 1, MPI_INT, 0, isle.id, isle.tcomm, MPI_STATUS_IGNORE);

    LOG(4, 0, 0, "rank %d got send size %d\r\n", isle.id, send_size);

    // initialize the set of process send channels

    isle.senders.clear();
    isle.senders.resize(send_size);

    // complete the MPI_Send of island ids to which this island will send ...
    
    MPI_Recv(&isle.senders[0], send_size, MPI_INT, 0, isle.id+t.world_size, isle.tcomm, MPI_STATUS_IGNORE);

    LOG(4, 0, 0, "topology %d: island %d got %lu imports from -> ", t.id, isle.id, isle.senders.size());
    
    if(isle.senders.empty()) {
        LOG(4, 0, 0, "[X]\r\n");
    } else {
        for(int i=0; i<isle.senders.size(); i++) { LOG(4, 0, 0, "[%d]", isle.senders[i]); }
        LOG(4, 0, 0, "\r\n");
    }

    // complete the MPI_Send of queue size for island import

    MPI_Recv(&rec_size, 1, MPI_INT, 0, isle.id+(t.world_size*2), isle.tcomm, MPI_STATUS_IGNORE);

    LOG(4, 0, 0, "rank %d got receive size %d\r\n", isle.id, send_size);

    // initialize the set of process receive channels

    isle.receivers.clear();
    isle.receivers.resize(rec_size);

    // complete the MPI_Send of island ids from which this island will receive ...

    MPI_Recv(&isle.receivers[0], rec_size, MPI_INT, 0, isle.id+(t.world_size*3), isle.tcomm, MPI_STATUS_IGNORE);

    LOG(4, 0, 0, "topology %d: island %d got %lu exports to ->   ", t.id, isle.id, isle.receivers.size());
    if(isle.receivers.empty()) {
        LOG(4, 0, 0, "[X]\r\n");
    } else {
        for(int i=0; i<isle.receivers.size(); i++) { LOG(4, 0, 0, "[%d]", isle.receivers[i]); }
        LOG(4, 0, 0, "\r\n");
    }
    
}

#pragma mark FUNCTION: topology_distribute()

// this function accepts a topology index, the associated topology as described in the                  |
// corresponding datatype (@topology{}), and an MPI communicator.                                       |
//                                                                                                      |
// successful execution of this function results in the insantiation of (blocking) MPI_Send channels    |
// for each island (process), over which the island metadata will be sent. a reciprocal MPI_Recv call   |
// is (must be) instantiated by each island to complete the distribution (see @topology_apply()).       |
//                                                                                                      |
// the @topology{} datatype stores metadata describing each island's specific send\receive operations   |
// in a vector of @channel{} with each element corresponding to a process (island) rank.                |
//                                                                                                      |
// performs an iteration for each island in the topology, and issues a send operation to                |
// each.  metadata describing the island's send\recv channels is queued                                 |
// as std::vector<int> of process (island) identifiers as mapped from the topology representation.      |
//                                                                                                      |
// @topology{}                                                                                          |
//   ↳ vector<@channel{}> channels[0..world_size]                                                       |
//      ↳ senders[0..send_size]                                                                         |
//      ↳ receivers[0..recv_size]                                                                       |
//                                                                                                      |
// *** this function is only intended to be called by a designated source process (rank 0) ***          |
//
//

void topology::distribute(island &isle) {

    // get number of mpi processes and the current process rank, and check to make sure
    // only rank 0 performs the topology distribution ...
    
    if(isle.id != 0) { return; }
    
    // issue an MPI_Send of the required i\o sockets for the provided corresponding process (island) ...
    
    for(int i=0; i<this->world_size; i++) {
        
        int send_size = (int)this->channels[i].senders.size();
        int rec_size = (int)this->channels[i].receivers.size();
            
        LOG(5, 0, 0, "sending topology %d, receive queue size %d to rank %d ...\r\n", this->id, send_size, i);
        
        MPI_Send(&send_size, 1, MPI_INT, i, i, isle.tcomm); // needed to size target senders vector
        MPI_Send(&this->channels[i].senders[0], send_size, MPI_INT, i, i+this->world_size, isle.tcomm);
        
        LOG(5, 0, 0, "sending topology %d, send queue size %d to rank %d ...\r\n", this->id, rec_size, i);
        
        MPI_Send(&rec_size, 1, MPI_INT, i, i+(this->world_size*2), isle.tcomm); // needed to size target receivers vector
        MPI_Send(&this->channels[i].receivers[0], rec_size, MPI_INT, i, i+(this->world_size*3), isle.tcomm);
        
    }
    
}

//#pragma mark FUNCTION: topologies::evolve()
//
//// create a new generation of topologies. accepts a reference to an array of
//// topology representations and a corresponding matrix size.
//// performs selection -> recombination -> mutation and returns an array of
//// evolved topology representations.
//
//std::vector<topology> topologies::evolve(island &isle) {
//    
//    // initialize parent selection ...
//    
//    std::vector<topology> children;
//    
//    if(isle.id != 0) { return children; }
//    
//    children.resize(config::topo_lambda);
//
//    std::reverse(this->population.begin(), this->population.end());
//    
//    topologies::calculate::cpd(*this);
//    
//    for(int n = 0; n < config::topo_lambda; n++) { // loop lambda
//
//            topology child;
//        
//            child.id = config::topo_mu+n;
//        
//            child.channels[n].senders.clear();
//            child.channels[n].receivers.clear();
//            
//            LOG(3, isle.id, 0, "creating topo kids\r\n");
//
//            topology t1 = topologies::select::parent(*this);
//            topology t2 = topologies::select::parent(*this);
//            
//            child.fitness = (t1.fitness + t2.fitness) / 2; // an estimated fitness, average of the parents
//            
//            // calculate an adjacency matrix for each parent's associated topology for use in
//            // generating child topology ...
//            
//            LOG(3, isle.id, 0, "parents t1=%2.10f,t2=%2.10f ...\r\n", t1.fitness, t2.fitness);
//            
//            std::vector<std::vector<int>> m1 = topology::create::matrix(t1);
//            std::vector<std::vector<int>> m2 = topology::create::matrix(t2);
//            
//            // recombine the parent adjacency matrices, initialize ...
//            
//            LOG(3, isle.id, 0, "recombining %d <-> %d ...\r\n", t1.id, t2.id);
//            
//            std::vector<std::vector<int>> child_matrix;
//            child_matrix.resize(this->world_size);
//            
//            int comm_count = 0;
//            int rec_count[this->world_size];
//            int snd_count[this->world_size];
//            
//            for(int i=0; i<this->world_size; i++) {
//                rec_count[i] = 0;
//                snd_count[i] = 0;
//            }
//            
//            LOG(3, isle.id, 0, "child %d initialized \r\n", child.id);
//            
//            // iterate row->column for each x,y element in the child matrix, and for each
//            // gene and randomly choose a parent from which to assign the value ...
//            
//            LOG(3, isle.id, 0, "performing child %d matrix crossover ...\r\n", child.id);
//            
//            while(comm_count == 0) { // failsafe to prevent empty matrix
//        
//                for(int i=0; i<m1.size(); i++) {  // child matrix row
//                    
//                    child_matrix[i].resize(this->world_size);
//                    
//                    for(int j=0; j<m2.size(); j++) { // child matrix column
//                        
//                        if(rec_count[j] >= config::migration_cap) {
//                            LOG(6, 0, 0, "migration cap limit reached for process %d\r\n", i);
//                            continue;
//                        }
//                        
//                        if(snd_count[i] >= config::send_cap) {
//                            LOG(6, 0, 0, "send cap limit reached for process %d\r\n", i);
//                            continue;
//                        }
//                        
//                        if(i != j) {  // we don't want an island sending migrants to itself
//                            
//                            if(rand()%2 == 1) { // coin flip, heads take the row index value from parent 1
//                                LOG(10, 0, 0, "assigning child m1[%d][%d] -> %d\r\n", i, j, m1[i][j]);
//                                child_matrix[i][j] = m1[i][j];
//                                if(m1[i][j] == 1) {
//                                    rec_count[j]++;
//                                    snd_count[i]++;
//                                }
//                            } else { // tails, take it from parent 2
//                                LOG(10, 0, 0, "assigning child m2[%d][%d] -> %d\r\n", i, j, m2[i][j]);
//                                child_matrix[i][j] = m2[i][j];
//                                if(m2[i][j] == 1) {
//                                    rec_count[j]++;
//                                    snd_count[i]++;
//                                }
//                            }
//                            
//                            if(child_matrix[i][j] == 1) { comm_count++; }
//                            
//                        } else {
//                            
//                            child_matrix[i][j] = 0;
//                            
//                        }
//                        
//                    } // child matrix column
//                    
//                    LOG(10, 0, 0, "------\r\n");
//                    
//                }  // child matrix row
//                
//            }  // end failsafe
//        
//            
//            // mutation, interate through the matrix and flip the bit with probability m,
//            // also considering any sparsity constraints ...
//            
//            for(int i=0; i<child_matrix.size(); i++) { // matrix row
//
//                for(int j=0; j<child_matrix[i].size(); j++) {  // matrix col
//
//                    if(rec_count[j] >= config::migration_cap) {
//                        LOG(6, 0, 0, "migration cap limit reached for process %d\r\n", i);
//                        continue;
//                    }
//                    
//                    if(snd_count[i] >= config::send_cap) {
//                        LOG(6, 0, 0, "send cap limit reached for process %d\r\n", i);
//                        continue;
//                    }
//                    
//                    if(rand()/(RAND_MAX+1.0) < config::mutation_rate) {
//                        
//                        LOG(3, isle.id, 0, "mutating child %d ...\r\n", child.id);
//                    
//                        if(child_matrix[i][j] == 0 && rand()/(RAND_MAX+1.0) < config::sparsity) {
//                            
//                            child_matrix[i][j] = 1;
//                            rec_count[j]++;
//                            snd_count[i]++;
//                            comm_count++;
//                            
//                        }
//                        
//                        if(child_matrix[i][j] == 1 && rand()/(RAND_MAX+1.0) > config::sparsity && comm_count > 1) {
//                            
//                            child_matrix[i][j] = 0;
//                            comm_count--;
//                            
//                        }
//                        
//                    }
//
//                } // matrix col
//
//            } // matrix row
//            
//            // convert the adjaceny matrix to a sender and receiver arrays for use with MPI send\recv ...
//        
//            topology::create::channels(child, child_matrix);
//
//            LOG(3, isle.id, 0, "child %d born!\r\n", child.id);
//            
//            children.push_back(child);
//            
//            //LOG(3, isle.id, 0, "evaluated child %d fitness = %f\r\n", isle.id, multi.topologies.children[n].fitness);
//        
//            LOG(3, 0, 0, "child %d created from matrix with %lu senders and %lu receivers\r\n", n, child.channels[n].senders.size(), child.channels[n].receivers.size());
//
//    } // loop lambda
//    
//    return children;
//    
//}

#pragma mark FUNCTION: topology_validate()

// sanity check to check topology distribution accuracy.

// TODO: incomplete implementation, needs work
// TODO: this implementation isn't ideal, because it performs a similar distribution
// TODO: find a more authoritative way to reference the topology data source

void topology::validate(island &isle) {
    
    bool result = true;
    
    // create new vectors to store migration data from the topology data source
    // for comparison ...
    
    std::vector<int> s;
    std::vector<int> r;
    
    int send_size;
    int recv_size;
    
    if(isle.id == 0) {
    
        // the root process is the authoritative source for populations, but we cannot access
        // its data context directly from the other nodes, so the current method is to send
        // the migration data again, to separate storage ...
        
        for(int i=1; i<this->world_size; i++) {

            send_size = (int)this->channels[i].senders.size();
            recv_size = (int)this->channels[i].receivers.size();
            
            s = this->channels[i].senders;
            r = this->channels[i].receivers;
            
            MPI_Send(&send_size, 1, MPI_INT, i, i, isle.tcomm);
            MPI_Send(&s[0], send_size, MPI_INT, i, i*(this->world_size), isle.tcomm);
            
            MPI_Send(&recv_size, 1, MPI_INT, i, i*(this->world_size*2), isle.tcomm);
            MPI_Send(&r[0], recv_size, MPI_INT, i, i*(this->world_size*3), isle.tcomm);
            
        }
        
        // store the process 0 (root) send and recv vectors
        
        s.resize(this->channels[0].senders.size());
        r.resize(this->channels[0].receivers.size());
        
        s = this->channels[0].senders;
        r = this->channels[0].receivers;
        
    } else {
        
        // gather the assigned senders and receivers for each non-root process
    
        MPI_Recv(&send_size, 1, MPI_INT, 0, isle.id, isle.tcomm, MPI_STATUS_IGNORE);
        s.resize(send_size);
        MPI_Recv(&s[0], send_size, MPI_INT, 0, isle.id*(this->world_size), isle.tcomm, MPI_STATUS_IGNORE);
        
        MPI_Recv(&recv_size, 1, MPI_INT, 0, isle.id*(this->world_size*2), isle.tcomm, MPI_STATUS_IGNORE);
        r.resize(recv_size);
        MPI_Recv(&r[0], recv_size, MPI_INT, 0, isle.id*(this->world_size*3), isle.tcomm, MPI_STATUS_IGNORE);
        
    }
    
    // check size of this island's receive operations against what we got from the topology source ...
    
    if(s.size() != isle.senders.size()) {
        result = false;
        LOG(4, 0, 0, "VALIDATION: topology %d failed: island %d sender count (%lu) does not match source (%lu).\r\n", this->id, isle.id, isle.senders.size(), s.size());
    }
    
    // for island imports, check to ensure there is a matching receive operation on the destination island ...
    
    for(int i=0; i<isle.senders.size(); i++) {

        int recv_from = isle.senders[i];

        if(std::find(s.begin(), s.end(), isle.id) != s.end()) {
            result = false;
            LOG(4, 0, 0, "VALIDATION: topology %d failed: island %d no matching send op found from %d\r\n", this->id, isle.id, recv_from);
        } else {
            LOG(4, 0, 0, "topology %d matching send op found [%d] -> [%d]\r\n", this->id, recv_from, isle.id)
        }
    }

    // check size of this island's send operations against what we got from the topology source ...
    
    if(r.size() != isle.receivers.size()) {
        result = false;
        LOG(4, 0, 0, "VALIDATION: topology %d failed: island %d receiver count (%lu) does not match source (%lu).\r\n", this->id, isle.id, isle.receivers.size(), r.size());
    }

    // for island exports, check to ensure there is a matching receive operation on the destination island ...
    
    for(int i=0; i<isle.receivers.size(); i++) {

        int send_to = isle.receivers[i];

        if(std::find(r.begin(), r.end(), isle.id) != r.end()) {
            result = false;
            LOG(4, 0, 0, "VALIDATION: topology %d failed: island %d no matching recv op found to %d\r\n", this->id, isle.id, send_to);
        } else {
            LOG(4, 0, 0, "topology %d matching rcv op found [%d] -> [%d]\r\n", this->id, isle.id, send_to)
        }

    }
    
    // return success if we passed all checks ...
    
    if(result == true) {
        LOG(4, 0, 0, " *** topology %d: island %d send receive operations VALIDATED ***\r\n", this->id, isle.id);
    }
    
}


//
//#pragma mark DATATYPE: @ea_objective_topology{}
//
//// population for the topology objective
//
//struct ea_objective_topology {
//    
//    std::vector<topology> population;
//    std::vector<topology> children;
//    
//};
//
//// population for the optimization benchmark objective
//
//#pragma mark DATATYPE: @ea_objective_solution{}
//
//struct ea_objective_solution {
//    
//    std::vector<solution> population;
//    std::array<double, DIM> offsets;
//    
//};

bool compare_topo_fitness(const topology &t1, const topology &t2) {
    
    return t1.fitness < t2.fitness;
    
}

bool is_zero(topology t) {
    return t.fitness == 0.0;
}

#endif /* topology_h */

//
//  topology.h
//  mpi_island_model_ea
//
//  Created by Bradley Morgan on 4/12/21.
//  Copyright © 2021 Bradley Morgan. All rights reserved.
//

#ifndef topology_h
#define topology_h

#pragma mark EA::META::OBJECTIVE::FUNCTION: <>create::channels()

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
                
                t.stats.recv_channels++;
                t.stats.send_channels++;
                
                t.stats.total_channels = t.stats.recv_channels + t.stats.send_channels;
                
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
//
// 

std::vector<std::vector<int>> topology::create::matrix(const topology &t) {
    
    // the passed vector comm should contain world_size group objects,
    // tracking the senders and receivers for each island
    
    LOG(5, 0, 0, "MATRIX: creating adjacency matrix for topology %d: ", t.id);
    
    std::vector<std::vector<int>> matrix;
    matrix.resize(mpi.size);
    
    for(int i=0; i<mpi.size; i++) { // matrix row, sized with world_size indices

        matrix[i].resize(mpi.size);
        
        LOG(9, 0, 0, "size[%d]=%lu, ", i, matrix[i].size());
        
        for(int j=0; j<mpi.size; j++) { // matrix column [0..world size]
            
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
            
            LOG(9, 0, 0, "[%d] ", matrix[i][j]);
            
        }
        
        LOG(9, 0, 0, "\r\n");
        
    }
    
    return matrix;
    
}

#pragma mark FUNCTION random_pairs()

std::vector<std::vector<int>> random_pairs() {
    
    std::vector<std::vector<int>> matrix;
    matrix.resize(mpi.size);
    int nsnd[mpi.size];
    int nrec[mpi.size];
    
    for(int i=0; i<mpi.size; i++) {
        nsnd[i] = 0;
        nrec[i] = 0;
        matrix[i].resize(mpi.size);
    }
    
    for(int i=0; i<mpi.size; i++) {
    
        int snd = 0;
        int rec = 0;
        
        while (snd == rec) {
            
            snd = rand()%(mpi.size-1);
            rec = rand()%(mpi.size-1);
            
            while(nsnd[snd] >= config::send_cap) {
                snd = rand()%(mpi.size-1);
            }
            
            while(nsnd[snd] >= config::send_cap) {
                rec = rand()%(mpi.size-1);
            }
        
        }
            
        if(prob_true(config::sparsity)) {
        
            matrix[snd][rec] = 1;
            nsnd[snd]++;
            nrec[rec]++;
            
            LOG(2, 0, 0, "%d: [%d(%d)][%d(%d)]:\r\n", i, snd, nsnd[snd], rec, nrec[rec]);
            
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
//    int send_max = world_size * ((config::send_cap * 1.0) / 100.0);
//    int recv_max = world_size * ((config::recv_cap * 1.0) / 100.0);
    int send_max = config::send_cap;
    int recv_max = config::recv_cap;
    int rec_count[world_size];
    int snd_count[world_size];
    
    // intialize empty matrix ...
    
    for(int i=0; i<world_size; i++) {
        rec_count[i] = 0;
        snd_count[i] = 0;
    }
    
    while(comm_count == 0) {  // avoid empty matrix if sparsity probability is low
    
        for(int i=0; i < world_size; i++) { // row

            matrix[i].resize(world_size);
            
            for(int j=0; j<world_size; j++) { // column
                
                // we don't want an island sending migrants to itself
                
                if(i == j) { matrix[i][j] = 0; continue; }
                    
                // check @config.txt[dim] max send\recv ...
                
                //if(rec_count[j] >= config::migration_cap) {
                if(rec_count[j] >= recv_max) {
                    LOG(6, 0, 0, "receive cap limit reached for process %d\r\n", i);
                    continue;
                }
                
                //if(snd_count[i] >= config::send_cap) {
                if(snd_count[i] >= send_max) {
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

#pragma mark FUNCTION: create_dyn_topology()

// used to create the initial topology population ...

void topology::create::dynamic(topology &t) {
    
    LOG(5, 0, 0, "generating dynamic topology %d world size %d\r\n", t.id, mpi.size);
    //std::vector<std::vector<int>> matrix = topology::create::dynamic_matrix(mpi.size);
    std::vector<std::vector<int>> matrix = random_pairs();
    
    topology::create::channels(t, matrix);
    
}

#pragma mark EA::META::OBJECTIVE::FUNCTION: topology_apply()

// issue reciprocal MPI_Recv calls for each matching MPI_Send from the root process
// (see @topology_distribute()) as determined fom the provided @island{} values.
//
// on success, all ranks will have received a copy of their send\receive operation metadata
//
// should be called from all islands
//

void topology::apply(island &isle, topology &t) {
    
    // rank authoritative source parses and queues the distribution of
    // island-specific migration channels to be assigned
    
    LOG(6, mpi.id, 0, "applying topology %d ...\r\n", t.id);
    
    if(mpi.id == 0) {
        
        LOG(6, mpi.id, 0, "distributing topology %d\r\n", t.id);
    
        t.init_multi();
        t.distribute(isle);
    
    }
    
    if(mpi.size == 0) { MPI_Comm_size(MPI_COMM_WORLD, &mpi.size); }
    
    // initialize *this* (current MPI rank) island's send and receive queue sizes, assume it's empty ...
    
    int send_size = 0;
    int rec_size = 0;
    
    // initiate MPI_Recv calls.
    
    LOG(6, 0, 0, "ISLAND %d of %d RECV INIT: inbound channel size from topology %d\r\n", mpi.id, mpi.size, t.id);
    
    // complete the MPI_Send of queue size for island export
    
    MPI_Recv(&send_size, 1, MPI_INT, 0, (mpi.id*10)+1, isle.tcomm, MPI_STATUS_IGNORE);

    LOG(4, 0, 0, "ISLAND %d of %d RECV END: %d inbound channels from topology %d\r\n", mpi.id, mpi.size, send_size, t.id);

    // initialize the set of process send channels

    isle.senders.clear();
    isle.senders.resize(send_size);

    // complete the MPI_Send of island ids to which this island will send ...
    
    LOG(6, 0, 0, "ISLAND %d of %d RECV INIT: inbound island ids from topology %d\r\n", mpi.id, mpi.size, t.id);
    
    MPI_Recv(&isle.senders[0], send_size, MPI_INT, 0, (mpi.id*10)+2, isle.tcomm, MPI_STATUS_IGNORE);

    LOG(4, 0, 0, "ISLAND %d of %d RECV END: %lu inbound island ids from topology %d\r\n", mpi.id, mpi.size, isle.senders.size(), t.id);
    
    if(isle.senders.empty()) {
        LOG(6, 0, 0, "[X]\r\n");
    } else {
        for(int i=0; i<isle.senders.size(); i++) { LOG(6, 0, 0, "[%d]", isle.senders[i]); }
        LOG(6, 0, 0, "\r\n");
    }

    // complete the MPI_Send of queue size for island import
    
    LOG(6, 0, 0, "ISLAND %d of %d RECV INIT: outbound channel size from topology %d\r\n", mpi.id, mpi.size, t.id);

    MPI_Recv(&rec_size, 1, MPI_INT, 0, (mpi.id*10)+3, isle.tcomm, MPI_STATUS_IGNORE);

    LOG(4, 0, 0, "ISLAND %d of %d RECV END: %lu outbound channels from topology %d\r\n", mpi.id, mpi.size, isle.senders.size(), t.id);

    // initialize the set of process receive channels

    isle.receivers.clear();
    isle.receivers.resize(rec_size);

    // complete the MPI_Send of island ids from which this island will receive ...

    LOG(6, 0, 0, "ISLAND %d of %d RECV INIT: outbound island ids from topology %d\r\n", mpi.id, mpi.size, t.id);
    
    MPI_Recv(&isle.receivers[0], rec_size, MPI_INT, 0, (mpi.id*10)+4, isle.tcomm, MPI_STATUS_IGNORE);
    
    LOG(4, 0, 0, "ISLAND %d of %d RECV END: %lu outbound island ids from topology %d\r\n", mpi.id, mpi.size, isle.senders.size(), t.id);

    LOG(7, 0, 0, "topology %d: island %d got %lu exports to ->   ", t.id, mpi.id, isle.receivers.size());
    if(isle.receivers.empty()) {
        LOG(6, 0, 0, "[X]\r\n");
    } else {
        for(int i=0; i<isle.receivers.size(); i++) { LOG(6, 0, 0, "[%d]", isle.receivers[i]); }
        LOG(6, 0, 0, "\r\n");
    }
    
}

#pragma mark FUNCTION: topology::distribute()

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
    
    if(mpi.id != 0) { return; }
    
    // issue an MPI_Send of the required i\o sockets for the provided corresponding process (island) ...
    
    for(int i=0; i<mpi.size; i++) {
        
        int send_size = (int)this->channels[i].senders.size();
        int rec_size = (int)this->channels[i].receivers.size();
        
        LOG(6, 0, 0, "ISLAND %d SEND INIT: sending %d send channels to %d of %d from topology %d\r\n", mpi.id, send_size, i, mpi.size, this->id);
        
        MPI_Send(&send_size, 1, MPI_INT, i, (i*10)+1, isle.tcomm); // needed to size target senders vector
        MPI_Send(&this->channels[i].senders[0], send_size, MPI_INT, i, (i*10)+2, isle.tcomm);
        
        LOG(4, 0, 0, "ISLAND %d SEND END: sent %d send channels to %d of %d from topology %d\r\n", mpi.id, send_size, i, mpi.size, this->id);
        
        LOG(6, 0, 0, "ISLAND %d SEND INIT: sending %d receive channels to %d of %d from topology %d\r\n", mpi.id, rec_size, i, mpi.size, this->id);
        
        MPI_Send(&rec_size, 1, MPI_INT, i, (i*10)+3, isle.tcomm); // needed to size target receivers vector
        MPI_Send(&this->channels[i].receivers[0], rec_size, MPI_INT, i, (i*10)+4, isle.tcomm);
        
        LOG(4, 0, 0, "ISLAND %d SEND END: sent %d receive channels to %d of %d from topology %d\r\n", mpi.id, rec_size, i, mpi.size, this->id);
        
    }
    
}

#endif /* topology_h */

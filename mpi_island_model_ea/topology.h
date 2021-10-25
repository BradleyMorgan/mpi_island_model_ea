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
//
//

std::vector<std::vector<int>> topology::create::matrix(const topology &t) {
    
    // the passed vector comm should contain world_size group objects,
    // tracking the senders and receivers for each island
    
    LOG(3, 0, 0, "MATRIX: creating adjacency matrix for topology %d: ", t.id);
    
    std::vector<std::vector<int>> matrix;
    matrix.resize(t.world_size);

    LOG(3, 0, 0, "size=%lu ", matrix.size());
    
    for(int i=0; i<t.world_size; i++) { // matrix row, sized with world_size indices

        matrix[i].resize(t.world_size);
        
        LOG(3, 0, 0, "size[%d]=%lu, ", i, matrix[i].size());
        
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
            
            LOG(3, 0, 0, "[%d] ", matrix[i][j]);
            
        }
        
        LOG(3, 0, 0, "\r\n");
        
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
    
        for(int i=0; i < world_size; i++) { // row

            matrix[i].resize(world_size);
            
            for(int j=0; j<world_size; j++) { // column
                
//                if(world_size / comm_count > config::sparsity) {
//                    return;
//                }
                
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

#pragma mark FUNCTION: create_dyn_topology()

// used to create the initial topology population ...

void topology::create::dynamic(topology &t) {
    
    LOG(5, 0, 0, "generating dynamic topology %d world size %d\r\n", t.id, t.world_size);
    std::vector<std::vector<int>> matrix = topology::create::dynamic_matrix(t.world_size);
    
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
    
    LOG(6, isle.id, 0, "applying topology %d ...\r\n", t.id);
    
    if(isle.id == 0) {
        
        LOG(6, isle.id, 0, "distributing topology %d\r\n", t.id);
    
        t.distribute(isle);
    
    }
    
    if(t.world_size == 0) { MPI_Comm_size(MPI_COMM_WORLD, &t.world_size); }
    
    // initialize *this* (current MPI rank) island's send and receive queue sizes, assume it's empty ...
    
    int send_size = 0;
    int rec_size = 0;
    
    // initiate MPI_Recv calls.
    
    LOG(6, 0, 0, "ISLAND %d of %d RECV INIT: inbound channel size from topology %d\r\n", isle.id, t.world_size, t.id);
    
    // complete the MPI_Send of queue size for island export
    
    MPI_Recv(&send_size, 1, MPI_INT, 0, (isle.id*10)+1, isle.tcomm, MPI_STATUS_IGNORE);

    LOG(4, 0, 0, "ISLAND %d of %d RECV END: %d inbound channels from topology %d\r\n", isle.id, t.world_size, send_size, t.id);

    // initialize the set of process send channels

    isle.senders.clear();
    isle.senders.resize(send_size);

    // complete the MPI_Send of island ids to which this island will send ...
    
    LOG(6, 0, 0, "ISLAND %d of %d RECV INIT: inbound island ids from topology %d\r\n", isle.id, t.world_size, t.id);
    
    MPI_Recv(&isle.senders[0], send_size, MPI_INT, 0, (isle.id*10)+2, isle.tcomm, MPI_STATUS_IGNORE);

    LOG(4, 0, 0, "ISLAND %d of %d RECV END: %lu inbound island ids from topology %d\r\n", isle.id, t.world_size, isle.senders.size(), t.id);
    
    if(isle.senders.empty()) {
        LOG(6, 0, 0, "[X]\r\n");
    } else {
        for(int i=0; i<isle.senders.size(); i++) { LOG(6, 0, 0, "[%d]", isle.senders[i]); }
        LOG(6, 0, 0, "\r\n");
    }

    // complete the MPI_Send of queue size for island import
    
    LOG(6, 0, 0, "ISLAND %d of %d RECV INIT: outbound channel size from topology %d\r\n", isle.id, t.world_size, t.id);

    MPI_Recv(&rec_size, 1, MPI_INT, 0, (isle.id*10)+3, isle.tcomm, MPI_STATUS_IGNORE);

    LOG(4, 0, 0, "ISLAND %d of %d RECV END: %lu outbound channels from topology %d\r\n", isle.id, t.world_size, isle.senders.size(), t.id);

    // initialize the set of process receive channels

    isle.receivers.clear();
    isle.receivers.resize(rec_size);

    // complete the MPI_Send of island ids from which this island will receive ...

    LOG(6, 0, 0, "ISLAND %d of %d RECV INIT: outbound island ids from topology %d\r\n", isle.id, t.world_size, t.id);
    
    MPI_Recv(&isle.receivers[0], rec_size, MPI_INT, 0, (isle.id*10)+4, isle.tcomm, MPI_STATUS_IGNORE);
    
    LOG(4, 0, 0, "ISLAND %d of %d RECV END: %lu outbound island ids from topology %d\r\n", isle.id, t.world_size, isle.senders.size(), t.id);

    LOG(7, 0, 0, "topology %d: island %d got %lu exports to ->   ", t.id, isle.id, isle.receivers.size());
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
    
    if(isle.id != 0) { return; }
    
    // issue an MPI_Send of the required i\o sockets for the provided corresponding process (island) ...
    
    for(int i=0; i<this->world_size; i++) {
        
        int send_size = (int)this->channels[i].senders.size();
        int rec_size = (int)this->channels[i].receivers.size();
        
        LOG(6, 0, 0, "ISLAND %d SEND INIT: sending %d send channels to %d of %d from topology %d\r\n", isle.id, send_size, i, this->world_size, this->id);
        
        MPI_Send(&send_size, 1, MPI_INT, i, (i*10)+1, isle.tcomm); // needed to size target senders vector
        MPI_Send(&this->channels[i].senders[0], send_size, MPI_INT, i, (i*10)+2, isle.tcomm);
        
        LOG(4, 0, 0, "ISLAND %d SEND END: sent %d send channels to %d of %d from topology %d\r\n", isle.id, send_size, i, this->world_size, this->id);
        
        LOG(6, 0, 0, "ISLAND %d SEND INIT: sending %d receive channels to %d of %d from topology %d\r\n", isle.id, rec_size, i, this->world_size, this->id);
        
        MPI_Send(&rec_size, 1, MPI_INT, i, (i*10)+3, isle.tcomm); // needed to size target receivers vector
        MPI_Send(&this->channels[i].receivers[0], rec_size, MPI_INT, i, (i*10)+4, isle.tcomm);
        
        LOG(4, 0, 0, "ISLAND %d SEND END: sent %d receive channels to %d of %d from topology %d\r\n", isle.id, rec_size, i, this->world_size, this->id);
        
    }
    
}

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
        LOG(3, 0, 0, "VALIDATION: topology %d failed: island %d sender count (%lu) does not match source (%lu).\r\n", this->id, isle.id, isle.senders.size(), s.size());
    }
    
    // for island imports, check to ensure there is a matching receive operation on the destination island ...
    
    for(int i=0; i<isle.senders.size(); i++) {

        int recv_from = isle.senders[i];

        if(std::find(s.begin(), s.end(), isle.id) != s.end()) {
            result = false;
            LOG(3, 0, 0, "VALIDATION: topology %d failed: island %d no matching send op found from %d\r\n", this->id, isle.id, recv_from);
        } else {
            LOG(5, 0, 0, "topology %d matching send op found [%d] -> [%d]\r\n", this->id, recv_from, isle.id)
        }
    }

    // check size of this island's send operations against what we got from the topology source ...
    
    if(r.size() != isle.receivers.size()) {
        result = false;
        LOG(3, 0, 0, "VALIDATION: topology %d failed: island %d receiver count (%lu) does not match source (%lu).\r\n", this->id, isle.id, isle.receivers.size(), r.size());
    }

    // for island exports, check to ensure there is a matching receive operation on the destination island ...
    
    for(int i=0; i<isle.receivers.size(); i++) {

        int send_to = isle.receivers[i];

        if(std::find(r.begin(), r.end(), isle.id) != r.end()) {
            result = false;
            LOG(3, 0, 0, "VALIDATION: topology %d failed: island %d no matching recv op found to %d\r\n", this->id, isle.id, send_to);
        } else {
            LOG(5, 0, 0, "topology %d matching rcv op found [%d] -> [%d]\r\n", this->id, isle.id, send_to)
        }

    }
    
    // return success if we passed all checks ...
    
    if(result == true) {
        LOG(3, 0, 0, " *** topology %d: island %d send receive operations VALIDATED ***\r\n", this->id, isle.id);
    }
    
}

template<> template<typename e> void objective<topology>::begin(objective_run &run, e &meta) {
    
    LOG(3, 0, 0, "BEGIN ISLAND %d META objective<topology> %d RUN %d\r\n", meta.variant.isle.id, this->id, this->run.id);
    
    this->run.begin();
    
    this->population.clear();
    this->population.resize(this->mu);
    
    this->log_begin(run, meta);
    
}

template<> template<typename e> void objective<topology>::begin(objective_eval &eval, e &meta) {
    
    LOG(3, 0, 0, "BEGIN ISLAND %d META objective<topology> %d EVAL %d -> ", meta.variant.isle.id, this->id, this->run.eval.id);
    
    this->run.eval.begin();
    this->log_begin(eval, meta);
    
}

template<> template<typename e> void objective<topology>::end(objective_run &run, e &meta) {
        
    LOG(3, 0, 0, "END ISLAND %d META objective<topology> %d RUN %d -> ", meta.variant.isle.id, this->id, this->run.id);
    
    this->run.end();
    
}

template<> template<typename e> void objective<topology>::end(objective_eval &eval, e &meta) {
    
    LOG(3, 0, 0, "END ISLAND %d META objective<topology> %d EVAL %d -> ", meta.variant.isle.id, this->id, this->run.eval.id);
    
    this->run.eval.end();
    
}

template<> template<typename e> void objective<topology>::log_begin(objective_run &run, e &meta) {
    
    if(meta.variant.isle.id != 0 || meta.topologies.run.id == 0) { return; }
    
    LOG(3, 0, 0, "LOG META OBJECTIVE %d RUN %d BEGIN\r\n", this->id, this->run.id);
    
}

template<> template<typename e> void objective<topology>::log_end(objective_run &run, e &meta) {

    if(meta.variant.isle.id != 0 || meta.topologies.run.id == 0) { return; }

    LOG(2, meta.variant.isle.id, 0, "LOG META OBJECTIVE %d RUN %d END\r\n", this-id, this->run.id);

    fprintf(config::topo_run_stats_out, "average_topo_fitness, global_best_topo_id, global_best_topo_rounds, global_best_topo_channels, global_best_topo_round_fitness, global_best_topo_fitness1, local_best_topo_fitness, global_best_topo_fitness2, average_local_best_topo_fitness, average_global_best_topo_fitness, t_id, t_rounds, t_channels, t_fitness\r\n");

    std::fprintf(config::topo_run_stats_out, "%d,%f,%f,%f,%f,%f,%d", meta.topologies.run.id, meta.topologies.run.stats.run_duration, meta.run.eval.stats.average_local_best_topo_fitness, meta.run.eval.stats.average_global_best_topo_fitness, meta.run.eval.stats.global_best_topo_fitness, meta.run.eval.stats.total_migrate_time, meta.run.stats.total_channels);

    fflush(config::topo_run_stats_out);

}

template<> template<typename e> void objective<topology>::log_begin(objective_eval &eval, e &meta) {
    
    if(meta.variant.isle.id != 0 || meta.topologies.run.eval.id == 0) { return; }
    
    LOG(3, 0, 0, "LOG META OBJECTIVE %d EVAL %d END\r\n", this->id, this->run.eval.id);
    
}

template<> template<typename e> void objective<topology>::log_end(objective_eval &eval, e &meta) {
    
    if(meta.variant.isle.id != 0 || meta.topologies.run.eval.id == 0) { return; }
    
    LOG(3, 0, 0, "LOG META OBJECTIVE %d EVAL %d END\r\n", this->id, this->run.eval.id);
    
}
    
template<> template<typename e, typename m, typename g> void objective<topology>::log_stats(objective_eval &eval, e &solver, m &meta, g &genome) {
    
    LOG(3, 0, 0, "STATS META OBJECTIVE %d EVAL %d\r\n", this->id, this->run.eval.id);
    
    if(meta.variant.isle.id != 0) { return; }
    
    log_fn_topology_stats(solver, meta, genome);
    
}

//template<> template<typename e> void objective<topology>::end(e &solver) {
//    
//    ea_end(solver, *this);
//    
//}

#endif /* topology_h */

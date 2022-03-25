#ifndef SIMPLIFICATION_AUX_STRUCTURES_H
#define SIMPLIFICATION_AUX_STRUCTURES_H


#include "basic_types/edge.h"
#include "basic_types/basic_wrappers.h"
const double Zero=1e-13;
struct Geom_Edge{
    ivect edge;
    double val;
   inline Geom_Edge(ivect &edge, double val){this->edge=edge; this->val=val;};
};

struct CompareEdge{

 

    bool operator()(Geom_Edge* e1,Geom_Edge* e2){

      if(fabs(e1->val-e2->val)>Zero)
      return (e1->val-e2->val)>Zero;
      else if(e1->edge[0]!=e2->edge[0])
      return e1->edge[0]>e2->edge[0];
      else
      return e1->edge[1]>e2->edge[1];
    }
};
typedef std::priority_queue<Geom_Edge*, std::vector<Geom_Edge*>, CompareEdge> edge_queue;

class contraction_parameters
{
public: 
  contraction_parameters(){checked_edges = 0;maximum_limit=0;simplification_counter=0;sum_edge_queue_sizes=0;QEM_based=false;parallel=false;all_edges=false;calculate_stats=false;}
    inline void increment_contracted_edges_counter() { checked_edges++; }
    inline int get_contracted_edges_num() { return checked_edges; }
    inline void set_maximum_limit(double l){this->maximum_limit=l;}
    inline double get_maximum_limit() { return maximum_limit; }
    inline void increment_counter(){simplification_counter++;}
    inline int get_counter(){return simplification_counter;}
    inline void add_edge_queue_size(int size){sum_edge_queue_sizes+=size;}
    inline int get_sum_edge_queue_sizes(){return sum_edge_queue_sizes;}
    inline void queue_criterion_QEM(){QEM_based=true;}
    inline void queue_criterion_length(){QEM_based=false;}
    inline bool is_QEM(){return QEM_based;}
    inline bool is_parallel(){return parallel;}
    inline bool no_limit(){return all_edges;}
    inline void parallel_compute(){parallel=true;}
    inline void sequential_compute(){parallel=false;}
    inline void contract_all_possible_edges(){all_edges=true;}
    inline void print_simplification_counters()
    {
        if(checked_edges >0)
            cerr<<"[STATS] contracted edges "<<checked_edges<<endl;
    }
    inline void calc_stats(){calculate_stats=true;}
protected:
  int checked_edges;
  double maximum_limit;
  int simplification_counter; //number of deleted triangles
  int sum_edge_queue_sizes;
  bool QEM_based;
  bool parallel;
  bool all_edges;
  bool calculate_stats;
};

#endif // SIMPLIFICATION_AUX_STRUCTURES_H
/*
    This file is part of the Terrain Trees library.

    Author(s): Riccardo Fellegara (riccardo.fellegara@gmail.com)

    This project has been supported by the Italian Ministry of Education and
    Research under the PRIN 2009 program, and by the National Science Foundation
    under grant number IIS-1116747.

    The Terrain Trees library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The Terrain Trees library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the Terrain Trees library.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef FORMANGRADIENTVECTOR_AUX_STRUCTURE_H
#define FORMANGRADIENTVECTOR_AUX_STRUCTURE_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <boost/function.hpp>

#include "basic_types/basic_wrappers.h"
#include "basic_types/edge.h"
#include "basic_types/triangle.h"
#include "triangle_gradient.h"
#include "ig.h"
#include "utilities/cli_parameters.h"

using namespace std;

namespace BitTwiddle {

// 32 bit version of popcount -- returns the number of bits in V that are set
template<typename T>
unsigned int popcount(T v)
{
    // from bit twiddling hacks: http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
    static unsigned int const masks[] = {
        0x55555555        // 0 -- 0b 01010101 01010101 01010101 01010101
        , 0x33333333        // 1 -- 0b 00110011 00110011 00110011 00110011
        , 0x0F0F0F0F        // 2 -- 0b 00001111 00001111 00001111 00001111
        , 0x01010101        // 3 -- 0b 00000001 00000001 00000001 00000001
    };

    static unsigned int const shifts[] = {
        1               // 0
        ,	2               // 1
        ,	4               // 2
        ,	8               // 3
    };

    //typename boost::make_unsigned<T>::type  v = t;

    v = v - ((v >> shifts[0]) & masks[0]);				                // reuse input as temporary
    v = (v & masks[1]) + ((v >> shifts[1]) & masks[1]);		            // temp
    return  ((v + (v >> shifts[2]) & masks[2]) * masks[3]) >> 24;           	// count
}
}

//typedef unsigned int uint;
#define CASES 512
#define OUTSIDECASES 1000

typedef enum { EDGE, TRIANGLE } TopType;

typedef map<itype,itype> map_2paths; //il primo int rappresenta il triangolo dove e' stato fermato il cammino, il secondo e' la label del massimo di partenza
typedef map<utype,map_2paths> leaves_2_desc_map;

typedef map<pair<itype,itype>,itype> map_paths_pair; //prima coppia è l'edge, l'int la label
typedef map<utype,map_paths_pair> leaves_1_desc_map;
typedef map<utype,map_paths_pair> leaves_2_asc_map;

class triple
{
protected:
    pair<itype,short> coppia;
    itype label;

public:
    triple(pair<itype,short> c, itype l)
    {
        coppia = c;
        label = l;
    }

    inline itype get_first() const { return coppia.first; }
    inline short get_second() const { return coppia.second; }
    inline pair<itype,short> get_pair() const { return coppia; }
    inline itype get_label() const { return label; }

    inline friend bool operator== (const triple &p, const triple &q)
    {
        if(p.coppia!=q.coppia) return false;
        if(p.label!=q.label) return false;
        return true;
    }
    inline friend bool operator!= (const triple &p, const triple &q) { return !(p==q); }
    inline bool operator<(const triple& s) const
    {
        return ((this->label < s.label) || (this->label == s.label && this->coppia < s.coppia));
    }
    inline bool operator>(const triple& s) const
    {
        return ((this->label > s.label) || (this->label == s.label && this->coppia > s.coppia));
    }
};

class quadruple : public triple
{
private:
    iNode* node;
public:
    quadruple(pair<itype,short> c, itype l, iNode* n) : triple(c,l) { node = n; }

    inline iNode* get_node() const { return node; }

    inline friend bool operator== (const quadruple &p, const quadruple &q)
    {
        if(p.coppia!=q.coppia) return false;
        if(p.label!=q.label) return false;
        if(p.node!=q.node) return false;
        return true;
    }
    inline friend bool operator!= (const quadruple &p, const quadruple &q) { return !(p==q); }

    inline bool operator<(const quadruple& s) const
    {
        return ((this->label < s.label) ||
                (this->label == s.label && this->coppia < s.coppia) ||
                (this->label == s.label && this->coppia == s.coppia && this->node < s.node));
    }
    inline bool operator>(const quadruple& s) const
    {
        return ((this->label > s.label) ||
                (this->label == s.label && this->coppia > s.coppia) ||
                (this->label == s.label && this->coppia == s.coppia && this->node > s.node));
    }
};

class desc1_mig_quadruple
{
private:
    Edge *e;
    itype last_v;
    short start_v_id;
    iNode *node;

public:
    desc1_mig_quadruple(Edge *edge, itype last, short l, iNode* n)
    {
        e = edge;
        last_v = last;
        start_v_id = l;
        node = n;
    }

    inline Edge* get_edge() const { return e; }
    inline itype get_last_v() const { return last_v; }
    inline short get_start_v() const { return start_v_id; }
    inline iNode* get_node() const { return node; }

    inline friend bool operator== (const desc1_mig_quadruple &p, const desc1_mig_quadruple &q)
    {
        if(*p.e!=*q.e) return false;
        if(p.last_v!=q.last_v) return false;
        if(p.start_v_id!=q.start_v_id) return false;
        if(p.node!=q.node) return false;
        return true;
    }
    inline friend bool operator!= (const desc1_mig_quadruple &p, const desc1_mig_quadruple &q) { return !(p==q); }
    inline bool operator<(const desc1_mig_quadruple& s) const
    {
        return ((this->last_v < s.last_v) ||
                (this->last_v == s.last_v && this->start_v_id < s.start_v_id) ||
                (this->last_v == s.last_v && this->start_v_id == s.start_v_id && this->node < s.node) ||
                (this->last_v == s.last_v && this->start_v_id == s.start_v_id && this->node == s.node && *(this->e) < *(s.e)) );
    }
    inline bool operator>(const desc1_mig_quadruple& s) const
    {
        return ((this->last_v > s.last_v) ||
                (this->last_v == s.last_v && this->start_v_id > s.start_v_id) ||
                (this->last_v == s.last_v && this->start_v_id == s.start_v_id && this->node > s.node) ||
                (this->last_v == s.last_v && this->start_v_id == s.start_v_id && this->node == s.node && *(this->e) > *(s.e)) );
    }
};

typedef set<triple> triple_set;
typedef map<utype,triple_set> leaves_asc_map;

typedef set<quadruple> set_asc_quadruple;
typedef map<utype,set_asc_quadruple> leaves_1_asc_mig_map;

typedef map<ivect, itype> simplices_map;

/// used by MIG -- start -- ///
typedef set<desc1_mig_quadruple> set_1paths;
typedef map<utype,set_1paths> leaves_1_desc_mig_map;
/// used by MIG -- end -- ///

struct desc1rels
{
    leaf_ETstar etstars;
    leaf_VTstar vtstars;

    inline void init(utype num_items)
    {
        vtstars.assign(num_items,-1);
    }
    inline void set_vtstar(utype v_pos, itype t) { vtstars.at(v_pos) = t; }
    inline utype get_vtstar(utype v_pos) { return vtstars.at(v_pos); }

    inline leaf_ETstar::iterator begin_etstars() { return etstars.begin(); }
    inline leaf_ETstar::iterator end_etstars() { return etstars.end(); }
};

struct asc2rels
{
    leaf_VV vvs;
    leaf_VTstar vtstars;

    inline void init(itype num_items)
    {
        vtstars.assign(num_items,-1);
        iset tmp_vv;
        vvs.assign(num_items,tmp_vv);
    }

    inline utype get_vtstar(utype v_pos) { return vtstars.at(v_pos); }

    inline void addVV(utype v_pos, utype v) { vvs.at(v_pos).insert(v); }
};

class local_VTstar_ET
{
protected:
    leaf_VTstar vtstars;
    leaf_ET ets;

public:
    inline void init(utype num_items)
    {
        vtstars.assign(num_items,-1);
    }

    inline leaf_VTstar& get_VTstars() { return vtstars; }
    inline leaf_ET& get_ETs() { return ets; }

    inline void add_VTstar(utype v_pos, utype t) { vtstars.at(v_pos) = t; }
    inline void add_ET(ivect &e, ET &et) { ets.insert(make_pair(e,et)); }

    inline utype get_VTstar(utype v_pos) { return vtstars.at(v_pos); }
    inline ET& get_ET(const ivect &e) { return ets.at(e); }

    inline leaf_ET::iterator find_ET(ivect &e) { return ets.find(e); }
    inline leaf_ET::iterator begin_ETs() { return ets.begin(); }
    inline leaf_ET::iterator end_ETs() { return ets.end(); }
};

namespace forman_aux_structures {

//questa struct permette di memorizzare in una sola lista tutti i simplessi
//con un overhead (posizioni dell'array inutilizzate) che varia a seconda del tipo di simplesso memorizzato
//il tipo di simplesso si ottiene gratis dal numero di posizioni utilizzate
struct simplex
{
    ivect simpl_id;
    vector<simplex*> coboundary;

    simplex(ivect &ids)
    {
        coboundary = vector<simplex*>();
        simpl_id.assign(ids.begin(),ids.end());
    }

    ~simplex()
    {
        simpl_id.clear();
        coboundary.clear();
    }

    inline void init(ivect &ids)
    {
        simpl_id.assign(ids.begin(),ids.end());

    }

    //for debug
    inline void print()
    {
        cout<<"simplex indices: ";
        for(utype i=0; i<simpl_id.size(); i++)
            cout<<simpl_id.at(i)-1<<" ";
        cout<<endl;
        if(coboundary.size() > 0)
        {
            cout<<"co boundary"<<endl;
            for(utype i=0; i<coboundary.size(); i++){
                for(utype j=0; j<coboundary.at(i)->simpl_id.size(); j++){
                    cout << coboundary.at(i)->simpl_id.at(j)-1 << " ";
                }
                cout << "-";
            }
            cout<<endl;
        }
    }
    inline void print_simplex()
    {
        cout<<"simplex indices: ";
        for(utype i=0; i<simpl_id.size(); i++)
            cout<<simpl_id[i]-1<<" ";
        cout<<endl;
    }

    //per avere la chiave della map
    inline string to_string()
    {
        //assumo gli indici gia' ordinati per field discendente
        stringstream ss;
        for(utype i=0; i<simpl_id.size(); i++)
        {
            if(i>0)
                ss << "_";
            ss << simpl_id[i];
        }

        return ss.str();
    }

};

typedef map<ivect, simplex *> lower_star_map;
typedef set<ivect, boost::function<bool(const ivect &, const ivect &)>> lower_star_set;

typedef ivect lower_star_set_SF;

struct gradient
{   
    simplices_map vertices_vector;
    simplices_map edges_vector;

    gradient()
    {
        vertices_vector = map<ivect, itype>();
        edges_vector = map<ivect, itype>();
    }

    ~gradient()
    {
        clear();
    }

    inline pair<simplices_map::iterator,bool> insert(const ivect &key, itype v)
    {
        if(key.size() == 1)
            return vertices_vector.insert(make_pair(key,v));
        else if(key.size() == 2)
            return edges_vector.insert(make_pair(key,v));

        cout<<"[gradient - insert] non dovrei mai arrivarci"<<endl;
        int a; cin>>a;
        return edges_vector.insert(make_pair(key,v));
    }

    inline void erase(const ivect &key, simplices_map::iterator it)
    {
        if(key.size() == 1)
            vertices_vector.erase(it);
        else if(key.size() == 2)
            edges_vector.erase(it);
    }

    inline simplices_map::iterator find(const ivect &key)
    {
        if(key.size() == 1)
            return vertices_vector.find(key);
        else if(key.size() == 2)
            return edges_vector.find(key);

        cout<<"[gradient - find] non dovrei mai arrivarci"<<endl;
        int a; cin>>a;
        return edges_vector.find(key);
    }

    inline simplices_map::iterator begin(const int key) //mi serve il vector per ritornare l'iteratore corretto
    {
        if(key == 1)
            return vertices_vector.begin();
        else if(key == 2)
            return edges_vector.begin();

        cout<<"[gradient - begin] non dovrei mai arrivarci"<<endl;
        int a; cin>>a;
        return edges_vector.begin();
    }

    inline simplices_map::iterator end(const ivect &key) //mi serve il vector per ritornare l'iteratore corretto
    {
        if(key.size() == 1)
            return vertices_vector.end();
        else if(key.size() == 2)
            return edges_vector.end();

        cout<<"[gradient - end] non dovrei mai arrivarci"<<endl;
        int a; cin>>a;
        return edges_vector.end();
    }

    inline simplices_map::iterator end(const int key) //mi serve il vector per ritornare l'iteratore corretto
    {
        if(key == 1)
            return vertices_vector.end();
        else if(key == 2)
            return edges_vector.end();

        cout<<"[gradient - end] non dovrei mai arrivarci"<<endl;
        int a; cin>>a;
        return edges_vector.end();
    }

    inline void add_local_gradients(gradient &vecs)
    {
        vertices_vector.insert(vecs.vertices_vector.begin(),vecs.vertices_vector.end());
        edges_vector.insert(vecs.edges_vector.begin(),vecs.edges_vector.end());
    }

    //mi serve durante l'estrazione dei cammini discendenti per salvare il gradiente dei simplessi usati nell'estrazione
    //il tutto e' usato in una fase di post-processing che mi calcola i cammini che sono ancora incompleti dopo una visita dell'indice pr-star
    inline void add_one_local_gradient(gradient &vecs, int size)
    {
        if(size==1)
            vertices_vector.insert(vecs.vertices_vector.begin(),vecs.vertices_vector.end());
        else if(size==2)
            edges_vector.insert(vecs.edges_vector.begin(),vecs.edges_vector.end());
    }

    inline void size()
    {
        cout<<"gradient entity numbers"<<endl;
        cout<<"vertices "<<this->vertices_vector.size()<<endl;
        cout<<"edge "<<this->edges_vector.size()<<endl;
    }

    inline void print()
    {
        cout<<"vertices"<<endl;
        for(simplices_map::iterator it=vertices_vector.begin(); it!=vertices_vector.end(); ++it)
        {
            for(itype i=0; i<it->first.size(); i++)
                cout<<it->first.at(i)<<" ";
            cout<<" - ";
            cout<<it->second;
            cout<<endl;
        }
        int a; cin>>a;
        cout<<"edges"<<endl;
        for(simplices_map::iterator it=edges_vector.begin(); it!=edges_vector.end(); ++it)
        {
            for(itype i=0; i<it->first.size(); i++)
                cout<<it->first.at(i)<<" ";
            cout<<" - ";
            cout<<it->second;
            cout<<endl;
        }
        cin>>a;
    }

    inline void clear()
    {
        simplices_map::iterator it= vertices_vector.begin();
        while(it!=vertices_vector.end())
        {
            vertices_vector.erase(it++);
        }
        simplices_map::iterator it2= edges_vector.begin();
        while(it2!=edges_vector.end())
        {
            edges_vector.erase(it2++);
        }
    }

    inline bool is_empty()
    {
        return ((vertices_vector.size() == 0) && (edges_vector.size() == 0) /*&& (faces_vector.size() == 0)*/);
    }

    inline void swap(gradient& gr)
    {
        this->vertices_vector.swap(gr.vertices_vector);
        this->edges_vector.swap(gr.edges_vector);
    }
};

struct forman_structure
{
    simplices_map critical_simplexes;
    gradient gradient_vectors;

    forman_structure()
    {
        critical_simplexes = simplices_map();
        gradient_vectors = gradient();
    }

    ~forman_structure()
    {
        simplices_map::iterator it= critical_simplexes.begin();
        while(it!=critical_simplexes.end())
            critical_simplexes.erase(it++);
        gradient_vectors.clear();
    }

    inline void clear()
    {
        simplices_map::iterator it= critical_simplexes.begin();
        while(it!=critical_simplexes.end())
            critical_simplexes.erase(it++);

        gradient_vectors.clear();
    }
};

class ig_paths
{
private:
    leaves_1_asc_mig_map paths_asc1cells;
    leaves_1_desc_mig_map paths_1cells;

public:
    ig_paths()
    {
        paths_asc1cells = leaves_1_asc_mig_map();
        paths_1cells = leaves_1_desc_mig_map();
    }

    ~ig_paths()
    {
        paths_asc1cells.clear();
        paths_1cells.clear();
    }

    inline bool visited_all()
    {
        return ((paths_asc1cells.size() == 0) &&
                (paths_1cells.size() == 0));
    }

    inline leaves_1_asc_mig_map& get_asc1cells_path() { return this->paths_asc1cells; }
    inline leaves_1_desc_mig_map& get_desc1cells_path() { return this->paths_1cells; }

    inline leaves_1_asc_mig_map::iterator find_asc1cells(utype key) { return this->paths_asc1cells.find(key); }
    inline leaves_1_desc_mig_map::iterator find_desc1cell(utype key) { return this->paths_1cells.find(key); }

    inline leaves_1_asc_mig_map::iterator begin_asc1cells() { return this->paths_asc1cells.begin(); }
    inline leaves_1_desc_mig_map::iterator begin_desc1cell() { return this->paths_1cells.begin(); }

    inline leaves_1_asc_mig_map::iterator end_asc1cells() { return this->paths_asc1cells.end(); }
    inline leaves_1_desc_mig_map::iterator end_desc1cell() { return this->paths_1cells.end(); }

    inline void erase(leaves_1_asc_mig_map::iterator it) { this->paths_asc1cells.erase(it); }
    inline void erase(leaves_1_desc_mig_map::iterator it) { this->paths_1cells.erase(it); }

    inline utype size_asc1cells() { return this->paths_asc1cells.size(); }
    inline utype size_desc1cell() { return this->paths_1cells.size(); }
};

class critical_clusters
{
private:
    uvect v_labels;
    uvect t_labels;
    vector<uset> vPerL;
    vector<uset> tPerL;
    utype ind, clusterized_min, clusterized_max;

public:

    critical_clusters(utype num_v, utype num_t, utype min_num, utype max_num)
    {
        ind = 1;
        v_labels = uvect(num_v,0);
        t_labels = uvect(num_t,0);
        vPerL = vector<uset>(min_num+1,uset());
        tPerL = vector<uset>(max_num+1,uset());

        clusterized_max = clusterized_min = 0;
    }
    inline void set_v_label(utype pos, utype l) { v_labels[pos-1] = l; }
    inline void set_t_label(utype pos, utype l) { t_labels[pos-1] = l; }
    inline void add_vPerL(utype pos, utype l) { vPerL[pos].insert(l); }
    inline void add_tPerL(utype pos, utype l) { tPerL[pos].insert(l); }
    inline void add_vPerL(utype pos, uset &s) { vPerL[pos].insert(s.begin(),s.end()); }
    inline void add_tPerL(utype pos, uset &s) { tPerL[pos].insert(s.begin(),s.end()); }

    inline utype get_v_label(utype pos) { return v_labels[pos-1]; }
    inline utype get_t_label(utype pos) { return t_labels[pos-1]; }

    inline utype get_vPerL_size(utype pos) { return vPerL[pos].size(); }
    inline utype get_tPerL_size(utype pos) { return tPerL[pos].size(); }

    inline void clear_vPerL(utype pos) { vPerL[pos].clear(); }
    inline void clear_tPerL(utype pos) { tPerL[pos].clear(); }

    inline uset& get_vPerL(utype pos) { return vPerL[pos]; }
    inline uset& get_tPerL(utype pos) { return tPerL[pos]; }

    inline void reset_counter() { ind = 1; }
    inline utype get_counter_value() { return ind; }
    inline void increment_counter() { ind++; }

    inline void increment_cmin_counter() { clusterized_min++; }
    inline void increment_cmax_counter() { clusterized_max++; }
    inline utype get_cmin_counter_value() { return clusterized_min; }
    inline utype get_cmax_counter_value() { return clusterized_max; }
};

typedef pair<int, int> simplex_id; //highest_vertex_index-dim

}

#endif // FORMANGRADIENTVECTOR_AUX_STRUCTURE_H

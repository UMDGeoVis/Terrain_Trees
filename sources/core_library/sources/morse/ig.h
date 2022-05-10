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

#ifndef IG_H
#define IG_H

#include <set>
#include <vector>
#include <map>

#include "basic_types/basic_wrappers.h"

using namespace std;

class IG_Node{
protected:
    itype critical_index;

public:
    inline itype get_critical_index() const {return critical_index;}
    virtual ~IG_Node() { }

    inline friend std::ostream& operator<< (std::ostream &out, IG_Node &q)
    {
        out << "("<<q.critical_index<< ")";
        return out;
    }
};

class Arc{

    IG_Node* node_i;
    IG_Node* node_j; //j=i+1;

    itype simplex_to_node_i;
    itype simplex_to_node_j;

    int label;

public:
    inline Arc() {}
    virtual ~Arc() { }
    inline Arc(IG_Node* node1, itype simplex_i, IG_Node* node2, itype simplex_j){node_i=node1; node_j=node2;label=1;
                                                                             simplex_to_node_i=simplex_i; simplex_to_node_j=simplex_j;}

    inline IG_Node* getNode_i() const {return node_i;}
    inline IG_Node* getNode_j() const {return node_j;}
    inline int getLabel(){return label;}
    inline void setLabel(int lab){label=lab;}
    inline itype getSimplexi(){return simplex_to_node_i;}
    inline itype getSimplexj(){return simplex_to_node_j;}

    inline friend std::ostream& operator<< (std::ostream &out, Arc &q)
    {
        out << "E["<<*q.node_i<< " " << *q.node_j << "] ";
        out << "I["<<q.simplex_to_node_i<< " " << q.simplex_to_node_j << "] ";
        out << "L["<<q.label<< "]";
        return out;
    }
};

class nNode : public IG_Node
{
private:
  set<Arc*> arcs;

public:
  inline nNode(itype ci){critical_index=ci;}
  virtual ~nNode() { arcs.clear(); }

  inline bool addArc(Arc*  arc){ return arcs.insert(arc).second; }
  inline set<Arc*>& getArcs(){ return arcs; }
  inline vector<Arc*> get_vector_Arcs(){ return vector<Arc*>(arcs.begin(), arcs.end()); }
  inline void removeArc(Arc* arc_to_remove){ arcs.erase(arc_to_remove); }
  inline void removeArc(set<Arc*>::iterator arc_to_remove){ arcs.erase(arc_to_remove); }
  inline void clear_arcs(){ arcs.clear(); arcs=set<Arc*>(); }

  inline set<Arc*>::iterator begin() { return arcs.begin(); }
  inline set<Arc*>::iterator end() { return arcs.end(); }
  inline int size() { return arcs.size(); }

};

class iNode : public IG_Node
{
    set<Arc*> arcs_up;
    set<Arc*> arcs_down;

    pair<itype,itype> edge_id;

public:
    inline iNode(itype edge){ critical_index = edge; edge_id = make_pair(-1,-1); }

    /// gli archi rimanenti vengono rimossi dalle selle
    virtual ~iNode()
    {
        for(set<Arc*>::iterator it=arcs_up.begin(); it!=arcs_up.end(); ++it)
            delete *it;
        for(set<Arc*>::iterator it=arcs_down.begin(); it!=arcs_down.end(); ++it)
            delete *it;
        arcs_up.clear();
        arcs_down.clear();
    }

    inline void add_edge_id(itype t1, itype t2){edge_id=pair<itype,itype>(t1,t2);}
    inline pair<itype,itype> get_edge_id(){return edge_id;}
    inline bool addArc(bool up, Arc* arc){
        if(up) return arcs_up.insert(arc).second;
        else return arcs_down.insert(arc).second;
    }
    inline set<Arc*>& getArcs(bool up)
    {
        if(up) return arcs_up;
        else return arcs_down;
    }
    inline vector<Arc*> get_vector_Arcs(bool up)
    {
        if(up) return vector<Arc*>(arcs_up.begin(), arcs_up.end());
        else return vector<Arc*>(arcs_down.begin(), arcs_down.end());
    }
    inline void removeArc(bool up, Arc* arc_to_remove){
        if(up) arcs_up.erase(arc_to_remove);
        else arcs_down.erase(arc_to_remove);
    }
    inline void removeArc(bool up, set<Arc*>::iterator arc_to_remove){
        if(up) arcs_up.erase(arc_to_remove);
        else arcs_down.erase(arc_to_remove);
    }
    inline void clear_arcs(){ arcs_up.clear(); arcs_down.clear(); }

    inline friend std::ostream& operator<< (std::ostream &out, iNode &q)
    {
        out << "c_index("<<q.critical_index<< ") et_id("<<q.edge_id.first<<" "<<q.edge_id.second<<")";
        return out;
    }
};


class IG
{

private:
    map<itype, nNode*> minima; /// key is v_id
    map<pair<itype,itype>, iNode*> saddle; /// key is <max_et_id,min_et_id>
    map<itype, nNode*> maxima; /// key is t_id

    vector<set<Arc*> > arcs;

public:
    inline IG(){arcs =vector<set<Arc*> >(2,set<Arc*>());}

    inline map<itype, nNode*>& getMinima(){return minima;}
    inline map<pair<itype,itype>, iNode*>& getSaddles(){return saddle;}
    inline map<itype, nNode*>& getMaxima(){return maxima;}

    inline utype minima_number() { return minima.size(); }
    inline utype saddle_number() { return saddle.size(); }
    inline utype maxima_number() { return maxima.size(); }
    inline utype critica_number() { return minima.size() + saddle.size() + maxima.size(); }

    inline set<Arc*>& getLevelArcs(int lvl){return arcs[lvl];}
    inline utype get_arcs_number() { return (arcs[0].size() + arcs[1].size()); }

    inline bool add_minimum(itype key, nNode* node) {return minima.insert(make_pair(key,node)).second; }
    inline bool add_saddle(pair<itype,itype> key, iNode* node) {return saddle.insert(make_pair(key,node)).second; }
    inline bool add_maximum(itype key, nNode* node) {return maxima.insert(make_pair(key,node)).second; }

    inline nNode* find_minimum(itype& key)
    {
        map<itype, nNode*>::iterator it = minima.find(key);
        if(it == minima.end())
            return NULL;
        else
            return it->second;
    }

    inline iNode* find_saddle(ET& key)
    {
        map<pair<itype,itype>, iNode*>::iterator it = saddle.find(key);
        if(it == saddle.end())
            return NULL;
        else
            return it->second;
    }

    inline nNode* find_maximum(itype& key)
    {
        map<itype, nNode*>::iterator it = maxima.find(key);
        if(it == maxima.end())
            return NULL;
        else
            return it->second;
    }

    inline void removeArc(int lvl, Arc* arc_to_remove){
        arcs[lvl].erase(arc_to_remove);
    }

    inline void remove_minimum(int key, nNode* node) { minima.erase(key); delete node; }
    inline void remove_saddle(pair<itype,itype> key, iNode* node) { saddle.erase(key); delete node;}
    inline void remove_maximum(int key, nNode* node) { maxima.erase(key); delete node; }

    inline Arc* addArc(IG_Node* n1, int s_n1, IG_Node* n2, int s_n2, int lvl)
    {
        Arc* new_arc = new Arc(n1,s_n1,n2,s_n2);
        arcs[lvl].insert(new_arc);
        if(lvl == 0)
        {
            ((nNode*)n1)->addArc(new_arc);
            ((iNode*)n2)->addArc(true, new_arc);
        }
        else
        {
            //assert(lvl==1);
            ((iNode*)n1)->addArc(false, new_arc);
            ((nNode*)n2)->addArc(new_arc);
        }

        return new_arc;
    }

    inline Arc* already_connected(nNode* extrema, iNode* saddle)
    {
        for(set<Arc*>::const_iterator it=extrema->begin(); it!=extrema->end(); ++it)
            if((((*it)->getNode_i() == extrema && (*it)->getNode_j() == saddle) ||
               ((*it)->getNode_j() == extrema && (*it)->getNode_i() == saddle))
                    && (*it)->getLabel() > 0 )
                return *it;

        return NULL;
    }

    /// init again the arcs vector
    inline void init() { arcs =vector<set<Arc*> >(2,set<Arc*>()); }

    inline void clear()
    {
        for(map<itype,nNode*>::iterator it=minima.begin(); it!=minima.end(); )
        {
            delete it->second;
            minima.erase(it);
            it = minima.begin();

        }
        minima.clear();
//        cout<<"minima cleared"<<endl;

        for(map<itype,nNode*>::iterator it=maxima.begin(); it!=maxima.end(); )
        {
            delete it->second;
            maxima.erase(it);
            it = maxima.begin();
        }
        maxima.clear();
//        cout<<"maxima cleared"<<endl;

//        cout<<"saddle number: "<<saddle.size()<<endl;

        for(map<pair<itype,itype>,iNode*>::iterator it=saddle.begin(); it!=saddle.end(); )
        {
//            cout<<"edge-id: "<<it->first.first<<" "<<it->first.second<<" --- ";
//            cout<<*(it->second)<<endl;
            delete it->second;
//            cout<<"B"<<endl;
            saddle.erase(it);
//            cout<<"C"<<endl;
            it = saddle.begin();
//            cout<<"inner--> saddle number: "<<saddle.size()<<endl;
        }
        saddle.clear();
//        cout<<"saddle cleared"<<endl;

        /// gli archi rimanenti vengono deallocati quando rimuovo le selle
        arcs.clear();
//        cout<<"arcs cleared"<<endl;
    }

    inline void print_saddles()
    {
        for(map<pair<itype,itype>,iNode*>::iterator it=saddle.begin(); it!=saddle.end(); ++it)
        {
            cout<<"S("<<it->first.first<<" "<<it->first.second<<") -> "<<*(it->second)<<endl;
        }
    }

    inline void print_stats(bool verbose)
    {
        cerr<<"-- MIG statistics --"<<endl;
        cerr<<"critical: "<<minima_number()<<" "<<saddle_number()<<" "<<maxima_number()<<endl;
        cerr<<"tot critical points: "<<this->critica_number()<<endl;
        cerr<<"arcs number: "<<get_arcs_number()<<endl;

        if(verbose)
        {
            cerr<<"minimum-saddle arcs: "<<arcs[0].size()<<endl;
            cerr<<"saddle-maximum arcs: "<<arcs[1].size()<<endl;
            int arcs_up = 0, arcs_down = 0;
            int num_saddle_with_two_min = 0, num_saddle_with_two_max = 0;
            for(map<pair<itype,itype>,iNode*>::iterator it=saddle.begin(); it!=saddle.end(); ++it)
            {
                arcs_up += it->second->getArcs(true).size();
                if(it->second->getArcs(true).size() == 2)
                    num_saddle_with_two_min++;
                arcs_down += it->second->getArcs(false).size();
                if(it->second->getArcs(false).size() == 2)
                    num_saddle_with_two_max++;
            }
            cerr<<"a_u: "<<arcs_up<<" a_d: "<<arcs_down<<" tot: "<<(arcs_up+arcs_down)<<endl;
            cerr<<"saddle_with_two_min: "<<num_saddle_with_two_min<<endl;
            cerr<<"saddle_with_two_max: "<<num_saddle_with_two_max<<endl;
        }
    }
};

class Topo_Sempl
{

public:
    Arc* arc;
    coord_type val;
    int filt_s0,filt_s1;
    ivect filt_ex;
    int lvl;

    Topo_Sempl() { arc = NULL; val = -1; lvl = -1;filt_s0=filt_s1=-1;filt_ex=ivect(); }
    Topo_Sempl(Arc* arc, coord_type val, int lvl,int filt_s0,int filt_s1,ivect filt_ex){ this->arc = arc; this->val=val; this->lvl=lvl;this->filt_s0=filt_s0;this->filt_s1=filt_s1;this->filt_ex=filt_ex;}

    inline friend std::ostream& operator<< (std::ostream &out, Topo_Sempl &q)
    {
        out << "A("<<*q.arc<< ") V["<<q.val<<"] lvl["<<q.lvl<<"]";
        return out;
    }
};

struct sort_arcs_topo{
    bool operator()(Topo_Sempl &s1, Topo_Sempl &s2)
    {   
        if(s1.val!=s2.val)
           return s1.val > s2.val;
        else if(s1.lvl!=s2.lvl)
            return s1.lvl>s2.lvl;
      /*  else if(s1.arc->getNode_i()->get_critical_index()!=s2.arc->getNode_i()->get_critical_index())
            return s1.arc->getNode_i()->get_critical_index()>s2.arc->getNode_i()->get_critical_index();
        else
            return s1.arc->getNode_j()->get_critical_index()>s2.arc->getNode_j()->get_critical_index();
      */
        else if(s1.filt_s0!=s2.filt_s0)
            return s1.filt_s0>s2.filt_s0;
        else if(s1.filt_s1!=s2.filt_s1)
        {
           // cout<<"All the same"<<s1.arc->getLabel()<<" AND "<<s2.arc->getLabel()<<endl;
            return s1.filt_s1>s2.filt_s1;
        }
        else
        {
            for (int i=0;i<s1.filt_ex.size();i++)
            {
                   if( s1.filt_ex[i]!=s2.filt_ex[i])
                    return s1.filt_ex[i] > s2.filt_ex[i];
            }
            return s1.filt_ex[s1.filt_ex.size()-1] > s2.filt_ex[s1.filt_ex.size()-1];
        }
           

    }
};

#endif // IG_H

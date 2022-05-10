/*
    This file is part of the Terrain Trees library.

    Author(s): Riccardo Fellegara (riccardo.fellegara@gmail.com)
               Federico Iuricich (federico.iuricich@gmail.com)

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

#include "forman_gradient_computation.h"
#include "utilities/container_utilities.h"

void Forman_Gradient_Computation::compute_gradient_vector(Forman_Gradient &gradient, Node_V &n, Mesh &mesh,
                                                                       Spatial_Subdivision &division)
{
   // this->initial_filerting(mesh);
        _min=saddle=_max=0;
    this->compute_local_gradient_vector(gradient,n,mesh,division);
    cerr <<"minima("<< this->_min
        << ") saddles(" << this->saddle
        << ") maxima(" << this->_max<< ")" << endl;
    
   
}

void Forman_Gradient_Computation::compute_gradient_vector(Forman_Gradient &gradient, Node_T &n, Box &n_dom, Mesh &mesh,
                                                                       Spatial_Subdivision &division)
{
   // this->initial_filerting(mesh);
    _min=saddle=_max=0;
    this->compute_local_gradient_vector(gradient,n,n_dom,0,mesh,division);
    cerr <<"minima("<< this->_min
        << ") saddles(" << this->saddle
        << ") maxima(" << this->_max<< ")" << endl;
     
}

void Forman_Gradient_Computation::compute_local_gradient_vector(Forman_Gradient &gradient, Node_V &n, Mesh &mesh,
                                                                             Spatial_Subdivision &division)
{
    if (n.is_leaf())
    {
        compute_local_gradient_vector_leaf(gradient,n,mesh);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                this->compute_local_gradient_vector(gradient,*n.get_son(i),mesh,division);
            }
        }
    }
}

void Forman_Gradient_Computation::compute_local_gradient_vector(Forman_Gradient &gradient, Node_T &n, Box &n_dom, int level, Mesh &mesh,
                                                                             Spatial_Subdivision &division)
{
    if (n.is_leaf())
    {
        compute_local_gradient_vector_leaf(gradient,n,n_dom,mesh);
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(n_dom,level,i);
            int son_level = level +1;
            if(n.get_son(i)!=NULL)
            {
                this->compute_local_gradient_vector(gradient,*n.get_son(i),son_dom,son_level,mesh,division);
            }
        }
    }
}

void Forman_Gradient_Computation::compute_local_gradient_vector_leaf(Forman_Gradient &fg, Node_V &n, Mesh &mesh)
{
    /// if there are no vertices in the leaf we have nothing to do..
    if(!n.indexes_vertices())
        return;

    // auto foo = bind(&Forman_Gradient_Computation::cmp_filtered_simplices, this,_1,_2);

    itype v_start = n.get_v_start();
    itype v_end = n.get_v_end();

    leaf_VT vts = leaf_VT();
    vts.assign(v_end-v_start, ivect());

    lower_star_map tmp = lower_star_map();
    vector< lower_star_map > lower_stars = vector<lower_star_map>();
    lower_stars.assign(v_end-v_start,tmp);
    leaf_ET ets;

    extract_topological_relations(n,mesh,lower_stars,vts,ets);
    compute_gradient(fg,v_start,vts,ets,lower_stars,mesh);
}

void Forman_Gradient_Computation::compute_local_gradient_vector_leaf(Forman_Gradient &fg, Node_T &n, Box &n_dom, Mesh &mesh)
{
    itype v_start;
    itype v_end;

    n.get_v_range(v_start,v_end,n_dom,mesh); // we need to gather the vertices range..

    if(v_start == v_end) //no internal vertices..
        return;

   // auto foo = bind(&Forman_Gradient_Computation::cmp_filtered_simplices, this,_1,_2);

    leaf_VT vts = leaf_VT();
    vts.assign(v_end-v_start, ivect());

    lower_star_map tmp = lower_star_map();
    vector< lower_star_map > lower_stars = vector<lower_star_map>();
    lower_stars.assign(v_end-v_start,tmp);

    leaf_ET ets;

    extract_topological_relations(n,v_start,v_end,mesh,lower_stars,vts,ets);
    compute_gradient(fg,v_start,vts,ets,lower_stars,mesh);
}

void Forman_Gradient_Computation::compute_gradient(Forman_Gradient &fg, itype v_start, leaf_VT &vts, leaf_ET &ets,
                                                                vector<lower_star_map> &lower_stars, Mesh &mesh)
{
//  vector<lower_star_map> lower_star_group;
    ivect v =ivect();
    for(uint i=0; i< lower_stars.size(); i++)
    {
        itype current_v=v_start+i;
        v.clear();
        v.push_back(current_v);

          if(lower_stars[i].size() == 0)
        {
            _min++;
        }
        else{
            //          cout<<"step1"<<endl;
            forman_structure local_forman;         

    
            compute_gradient_in_a_vertex(v, lower_stars[i], vts[i], local_forman.critical_simplexes, local_forman.gradient_vectors, mesh);
             
            initialize_forman_gradient(fg,local_forman.gradient_vectors,ets,mesh);

        }     
        
    }

           for(vector< lower_star_map >::iterator iter = lower_stars.begin(); iter!=lower_stars.end(); ++iter)
    {
        //when I computed the gradient in n, I have to free the memory used by the lower stars
        for(lower_star_map::iterator it = iter->begin(); it != iter->end(); it++)
        {
            simplex * simpl = it->second;
            simpl->coboundary.clear();
            delete simpl;
        }
        iter->clear();
    }

}

void Forman_Gradient_Computation::compute_gradient_in_a_vertex(ivect &v, lower_star_map &lower_star, VT &vts,
                                                              simplices_map &critical_points, gradient &local_vectors, Mesh &mesh)
{
     vector<int> minimal_edge = vector<int>(0); //extract the smallest edge (it exists as the lower star is not empty)
    vector<vector<int> > other_edges = vector<vector<int> >(); //save the other edges here for later..

    for(map< vector<int>, simplex *>::iterator it = lower_star.begin(); it != lower_star.end(); it++)
    {
        if(it->first.size() == 2)
        {
            bool same=false;

            if(minimal_edge.size()==0)
            {
           
                minimal_edge = it->first;
            }
            else if(is_lower(it->first, minimal_edge, &same))
            {
                
                other_edges.push_back(minimal_edge);
                minimal_edge = it->first;
            }
            else
            {
                 
                other_edges.push_back(it->first);
            }
        }
    }

    //initialize the two "priority queues" using lists as I have to remove elements from pq_zero..
    //do not use std::priority_queue as the removal is much more inefficient
    list<vector<int> > pq_one = list<vector<int> >();
    list<vector<int> > pq_zero = list<vector<int> >();
    //cout<<v[0]<<" and "<<minimal_edge[1]<<endl;
  //  cout<<v[0]<<"'s other edges: "<<other_edges.size()<<endl;
    local_vectors.insert(v,minimal_edge[1]);  //add the first edge
    for(uint i=0; i<other_edges.size(); i++)
    {
        push_ordered(pq_zero, other_edges[i]);
    }
    
    vector<simplex*> &coboundary_faces = lower_star[minimal_edge]->coboundary;

    for(uint i=0; i<coboundary_faces.size(); i++)
    {
        if(num_unpared_faces(coboundary_faces[i]->simpl_id, local_vectors, critical_points,mesh) == 1)
        {
            push_ordered(pq_one, coboundary_faces[i]->simpl_id);
        }
    }

    other_edges.clear();

    vector<int> popped = vector<int>();
    vector<int> pair = vector<int>();
    vector<int> critico = vector<int>();

    while(pq_one.size() != 0 || pq_zero.size() != 0)
    {
  //      cerr<<v[0]<<"'s pq_zero size (before):"<<pq_zero.size()<<", pg_one size: "<<pq_one.size()<<endl;
        while(pq_one.size() != 0)
        {
            popped = pq_one.front();
            pq_one.pop_front();

            if(num_unpared_faces(popped, local_vectors, critical_points,mesh) == 0)
            {
                push_ordered(pq_zero, popped);
        //         cerr<<"checkpoint A"<<endl;
            }
            else
            {                
                //init a new part of the gradient and remove it frompq_zero
                pair = unique_pairable_face(popped, local_vectors, critical_points,mesh);
  //             cerr<<"checkpoint B"<<endl;
                if(pair.size()==2) // edge-case
                {
                    int t_id = get_t_id(popped,vts,mesh);
                    local_vectors.insert(pair,t_id);
   //                 cerr<<"Case A;";
                }
                else //otherwise I add only the linked simplex
                {
                    local_vectors.insert(pair,extract_missing_id(pair,popped));
     //               cerr<<"Case B;";
                }
        //        cerr<<endl;

  //              cerr<<endl;
               // cerr<<pq_zero.front()[0]<<"..."<<pq_zero.front()[1]<<endl;
                pq_zero.remove(pair);

                //add the co-boundary faces of simplex "popped" only if the have a single face to be paired with
                vector<simplex*> &coboundary_faces_popped = lower_star[popped]->coboundary;

                for(uint i=0; i<coboundary_faces_popped.size(); i++)
                {
                    if(num_unpared_faces(coboundary_faces_popped[i]->simpl_id, local_vectors, critical_points,mesh) == 1)
                    {
                        push_ordered(pq_one, coboundary_faces_popped[i]->simpl_id);
                    }
                }

                //add the co-boundary faces of simplex "pair" only if the have a single face to be paired with
                vector<simplex*> &coboundary_faces_pair= lower_star[pair]->coboundary;

                for(uint i=0; i<coboundary_faces_pair.size(); i++)
                {
                    if(num_unpared_faces(coboundary_faces_pair[i]->simpl_id, local_vectors, critical_points,mesh) == 1)
                    {
                        push_ordered(pq_one, coboundary_faces_pair[i]->simpl_id);
                    }
                }
            }
        }

        // cerr<<v[0]<<"'s pq_zero size (after):"<<pq_zero.size()<<", pg_one size: "<<pq_one.size()<<endl;
        if(pq_zero.size() != 0)
        {
            critico = pq_zero.front();
            pq_zero.pop_front();

            //the critical simplex is identified by the vertex with lower position index
            if(critico.size() == 2)
            {
                critical_points[critico] = critico.front();
//                cout<<"ONESADDLE: "<<critico.front()-1<<endl;
                saddle++;
            }
            if(critico.size() == 3)
            {
                critical_points[critico] = critico.front();
//                cout<<"MAXIMUM: "<<critico.front()-1<<endl;
//                cout<<critico[0]-1<<" "<<critico[1]-1<<" "<<critico[2]-1<<endl;
                _max++;
            }

            //add all simplices of the co-boundary in pq_one
            vector<simplex*> &coboundary_faces_critico= lower_star[critico]->coboundary;

            for(uint i=0; i<coboundary_faces_critico.size(); i++)
            {
                if(num_unpared_faces(coboundary_faces_critico[i]->simpl_id, local_vectors, critical_points,mesh) == 1)
                {
                    push_ordered(pq_one, coboundary_faces_critico[i]->simpl_id);
                }
            }

        }
    }

}


void Forman_Gradient_Computation::initial_filtering(Mesh &mesh)
{
    utype num_v = mesh.get_vertices_num();
   

    cout<<"vertices number: "<<num_v<<endl;
   

    filtration = uvect(num_v,0); //final filtration
    //componentBasedFiltration = vector<uvect>(num_fields,uvect(num_v,0)); //simulation of simplicity for each component
    vector<pair<coord_type,itype> > injectiveF(num_v);
 for(utype j=0;j<num_v; j++)
        {
            coord_type val=mesh.get_vertex(j+1).get_field(0);
            injectiveF[j]=pair<coord_type,itype>(val,j);
        }

        sort(injectiveF.begin(),injectiveF.end(),bind(&Forman_Gradient_Computation:: filtration_cmp, this,_1,_2));
        utype ind=0;
        for(auto p : injectiveF)
        {
            filtration[p.second]=ind++;
        }
       


}

void Forman_Gradient_Computation::initial_filtering_IA(Mesh &mesh)
{
    utype num_v = mesh.get_vertices_num();
   map<double, vector<int> > vert;
      cout<<"vertices number: "<<num_v<<endl;
   

     filtration = uvect(num_v,0); //final filtration
        for(int i=0; i< num_v; i++){
            vert[mesh.get_vertex(i+1).get_z()].push_back(i);
        }

        int count=0;
        for(auto m : vert){
            for(auto v : m.second){
                filtration[v]=count++;
            }
        }

 
    


    cout<<filtration[100]<<"; "<<filtration[59]<<endl;
}

void Forman_Gradient_Computation::reset_filtering(Mesh &mesh,ivect& original_vertex_indices, bool new_to_old){
    utype num_v = mesh.get_vertices_num();
   

    cout<<"vertices number: "<<num_v<<endl;
    uvect tmp = vector<utype>(num_v, 0);
    //final filtration

    if(new_to_old){
        cout<<"after reindexing"<<endl;
        for(utype j=0;j<num_v;j++){
            tmp[j]=filtration[original_vertex_indices[j]-1];
    
        }
    }
    else{
        cout<<"after simplification"<<endl;
        cout<<filtration.size()<<endl;
        cout<<original_vertex_indices.size()<<endl;
        for(utype j=0;j<original_vertex_indices.size();j++){
            if(original_vertex_indices[j]!=-1)
                tmp[original_vertex_indices[j]-1]=filtration[j];
        }
    }
    filtration = tmp;

}




void Forman_Gradient_Computation::extract_topological_relations(Node_V &n, Mesh &mesh, vector<lower_star_map> &lower_stars,
                                                                             leaf_VT &vts, leaf_ET &ets)
{
    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        this->extract_topological_relation(n,*t_id,mesh.get_triangle(*t_id),lower_stars,vts,ets,mesh);
    }
}

void Forman_Gradient_Computation::extract_topological_relations(Node_T &n, itype v_start, itype v_end, Mesh &mesh,
                                                                             vector<lower_star_map> &lower_stars, leaf_VT &vts, leaf_ET &ets)
{
    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        this->extract_topological_relation(n,v_start,v_end,*t_id,mesh.get_triangle(*t_id),lower_stars,vts,ets,mesh);
    }
}

void Forman_Gradient_Computation::extract_topological_relation(Node_V &n, itype t_id, Triangle &t, vector<lower_star_map> &lower_stars,
                                                                            leaf_VT &vts, leaf_ET &ets,Mesh &mesh)
{
    for(int v=0; v<t.vertices_num(); v++)
    {
        ///init della ET
        ivect e;
        t.TE(v,e);
        /// HERE I HAVE TO GET THE MAXIMUM FIELD VERTEX
        itype max_field_id = this->sort_filtered_vertices(e[0],e[1]) ? e[0] : e[1];
        this->add_ET_entry(e,max_field_id,t_id,n,ets);
        ///init VT
        itype v_id = t.TV(v);
        this->add_VT_entry(v_id,t_id,n,vts);
        ///init LOWER_STAR
        if(n.indexes_vertex(v_id))
           this->add_lower_star_entry(v_id,t,lower_stars[v_id-n.get_v_start()],mesh);
    }
}

void Forman_Gradient_Computation::extract_topological_relation(Node_T &n, itype v_start, itype v_end, itype t_id, Triangle &t,
                                                                            vector<lower_star_map> &lower_stars, leaf_VT &vts, leaf_ET &ets,Mesh &mesh)
{
    for(int v=0; v<t.vertices_num(); v++)
    {
        ///init della ET
        ivect e;
        t.TE(v,e);
        /// HERE I HAVE TO GET THE MAXIMUM FIELD VERTEX
        itype max_field_id = this->sort_filtered_vertices(e[0],e[1]) ? e[0] : e[1];
        if(n.indexes_vertex(v_start,v_end,max_field_id))
        {
            leaf_ET::iterator it = ets.find(e);
            if(it != ets.end())
            {
                ET &inside = it->second;
                inside.second = t_id;
            }
            else
            {
                ET new_entry = make_pair(t_id,-1);
                ets.insert(make_pair(e,new_entry));
            }
        }

        itype v_id = t.TV(v);
        if(n.indexes_vertex(v_start,v_end,v_id))
        {
            itype v_pos = v_id - v_start;
            ///init VT
            vts[v_pos].push_back(t_id);
            ///init LOWER_STAR
            this->add_lower_star_entry(v_id,t,lower_stars[v_pos],mesh);
        }
    }
}

bool Forman_Gradient_Computation::cmp_filtered_simplices(const ivect& lhs, const ivect& rhs)
{
    if(lhs.size() == rhs.size())
    {
        itype fValuesL;
        itype fValuesR;
        //itype similar_pair=-1;
        //cout<<"size of lower star"<<lhs.size()<<endl;
        for(utype i=0; i<lhs.size();i++)
        {
            fValuesL=filtration[lhs[lhs.size()-i-1]-1];
            fValuesR=filtration[rhs[lhs.size()-i-1]-1];
            if(fValuesL==fValuesR)
                continue;
            else
                return fValuesL<fValuesR;
        }

    }
    
    return lhs.size() < rhs.size();
}

utype Forman_Gradient_Computation::simplex_filtration(const ivect &simpl)
{
    utype filtr=0;

        for(auto v : simpl)
        {
            if(filtration[v-1] > filtr)
                filtr=filtration[v-1];
        }
    
    return filtr;
}


void Forman_Gradient_Computation::add_lower_star_entry(itype v_id, Triangle &t, lower_star_map &lower_stars,Mesh &mesh)
{
      vector<int> simpl =  vector<int>();
      vector<int> sub_simplex =  vector<int>();
     //  cout<<"start";
      simpl.push_back(v_id);
      int v_pos=t.vertex_index(v_id);

      for(int i=1;i<t.vertices_num();i++){
          int v2_id=abs(t.TV((v_pos+i)%t.vertices_num()));
          if(filtration[v_id-1]>filtration[v2_id-1])
          {
              simpl.push_back(v2_id);
             
          }
      }
  
      if(simpl.size() == 1) {
          //cout<<"v_id"<<v_id<<endl;
          return;}

     //  sort(simpl.begin(),simpl.end(),bind(&Forman_Gradient_SF::sort_filtered_vertices,this,_1,_2));
       sort_simplex(simpl,mesh);


      simplex* s = new simplex(simpl);
            if(lower_stars.find(simpl) == lower_stars.end())
            {
        
                lower_stars[simpl] = s; //metto in relazione il simplesso top con la sua rappresentazione
            }
   
            if(simpl.size()>2)
            {
                for(int i=1; i<simpl.size(); i++)
                {
                    //qua giro nel simplesso ed elimino i vertici uno ad uno per costruire i sottosimplessi (non elimino mai il primo vertice)
                    sub_simplex = simpl;
                    sub_simplex.erase(sub_simplex.begin()+i); // cancello la posizione i-esima

                    lower_star_map::iterator it = lower_stars.find(sub_simplex);
                    simplex *sopra;

                    if(it == lower_stars.end())
                    {
                        sopra = new simplex(sub_simplex);
                        lower_stars[sub_simplex] = sopra;
                    }
                    else
                    {
                        sopra = it->second;
                    }
   
                    push_in_coboundary(sopra, lower_stars[simpl]);
                                            
                }
            }
        
}

int Forman_Gradient_Computation::numPairableLowerStar(const ivect &next, const lower_star_map& sset, ivect& pair)
{
    int num=0;
    ivect sub_simplex =  ivect();
    for(utype i=0; i<next.size(); i++)
    {
        sub_simplex = next;
        sub_simplex.erase(sub_simplex.begin()+i);
        if(sset.find(sub_simplex) != sset.end())
        {
            num++;
            pair=sub_simplex;
        }
    }
    return num;
}

void Forman_Gradient_Computation::initialize_forman_gradient(Forman_Gradient &fg, gradient &grad,leaf_ET &ets, Mesh &mesh)
{
  //  cout<<"set_VE"<<endl;
    for(simplices_map::iterator it=grad.begin(1);it!=grad.end(1);it++){
        fg.set_VE(it->first.front(),it->second,ets,mesh);
    }
 //  cout<<"set_ET"<<endl;
    for(simplices_map::iterator it = grad.begin(2); it != grad.end(2); it++)
    {
        fg.set_ET(it->second, it->first, mesh);
    }
}

void Forman_Gradient_Computation::push_in_coboundary(simplex *sface, simplex *co_sface)
{
    uint i=0;
    for(; i<sface->coboundary.size(); i++){

        if(sface->coboundary[i]->simpl_id == co_sface->simpl_id)
        {
            break;
        }
    }

    if(i==sface->coboundary.size())
    {
        sface->coboundary.push_back(co_sface);
    }
}

itype Forman_Gradient_Computation::extract_missing_id(ivect &smaller, ivect &bigger)
{
    itype ret=-1;

    for(ivect::iterator it=bigger.begin(); it!=bigger.end(); ++it)
    {
        bool found=false;
        for(ivect::iterator it2=smaller.begin(); it2!=smaller.end(); ++it2)
        {
            if(*it==*it2)
            {
                found = true;
                break;
            }
        }
        if(!found)
        {
            ret = *it;
            break;
        }
    }
    return ret;
}

void Forman_Gradient_Computation::add_VT_entry(itype v_id, itype t_id, Node_V &n, leaf_VT &vts)
{
    if(n.indexes_vertex(v_id))
    {
        itype v_pos = v_id - n.get_v_start();
        vts[v_pos].push_back(t_id);
    }
}

void Forman_Gradient_Computation::add_ET_entry(ivect &e, itype max_field_id, itype t_id, Node_V &n, leaf_ET &ets)
{
    if(n.indexes_vertex(max_field_id))
    {
        leaf_ET::iterator it = ets.find(e);
        if(it != ets.end())
        {
            ET &inside = it->second;
            inside.second = t_id;
        }
        else
        {
            ET new_entry = make_pair(t_id,-1);
            ets.insert(make_pair(e,new_entry));
        }
    }
}




bool Forman_Gradient_Computation::is_lower(const vector<int> &simplex1, const vector<int> &simplex2, bool* same){
    
      int similar_pair=-1;
       if(simplex1.size() < simplex2.size()) return true;
        if(simplex1.size() > simplex2.size()) return false;

    //  if(simplex1.size()==simplex2.size()){
          for(uint i=1; i<simplex1.size(); i++) /// new parto da 1 per alla posizione 0 ho sempre il vertice al centro del lower star
          {
              if(filtration[simplex1[i]-1] ==filtration[simplex2[i]-1])
              {
                  if(simplex1[i]!=simplex2[i])similar_pair=i;
                  continue;
              }
              else if(filtration[simplex1[i]-1] >filtration[simplex2[i]-1]) {
                  return false;
              }
              else
              {
                  return true;
              }
          }

          if(simplex1.back() == simplex2.back() && similar_pair == -1)
              *same=true;
          else
              *same=false;

    //  }
   //   cout<<simplex1[similar_pair]<<"____" <<simplex2[similar_pair];
//      
          return (simplex1[similar_pair] < simplex2[similar_pair]);
}

int Forman_Gradient_Computation::get_t_id(const vector<int> &popped, vector<int> &vt, Mesh &mesh)
{
    for(vector<int>::iterator it=vt.begin(); it!=vt.end(); ++it)
    {
        Triangle& top = mesh.get_triangle(*it);
//        cout<<top<<endl;
        int count_v = 0;

        for(int v=0; v<top.vertices_num(); v++)
        {
            for(vector<int>::const_iterator it2=popped.begin(); it2!=popped.end(); ++it2)
            {
                if(abs(top.TV(v))==*it2)
                {
                    count_v++;
                    break;
                }
            }

            if(count_v != v + 1)
                break;
        }

        if(count_v == top.vertices_num())
            return *it;
    }

    //non ci dovrebbe mai arrivare se i due vector sono inizializzati correttamente
    cout<<"[get_t_id] non dovrei arrivarci"<<endl;
    int a; cin>>a;
    return -1;
}

vector<int> Forman_Gradient_Computation::unique_pairable_face(vector<int> &s_face, gradient &local_vectors, simplices_map &critical_points, Mesh &mesh)
{
    vector<int> simpl_face = vector<int>();
    vector<int> simpl_face_second;
    for(unsigned int i=1; i<s_face.size(); i++)
    {
        simpl_face = s_face;
        simpl_face.erase(simpl_face.begin()+i);

        if(local_vectors.find(simpl_face) != local_vectors.end(simpl_face)) continue;
        if(critical_points.find(simpl_face) != critical_points.end()) continue;

        int j=-1;
        if(simpl_face.size() > 1){
           // cout<<simpl_face.size()<<endl;
            for(j=1; j<(int)simpl_face.size();j++){
                simpl_face_second = simpl_face;

                simpl_face_second.erase(simpl_face_second.begin()+j);

                simplices_map::iterator it = local_vectors.find(simpl_face_second);
                if(it != local_vectors.end(simpl_face_second)){
                    vector<int> simpl;

                    if(simpl_face_second.size() == 1){
                        simpl = simpl_face_second;
                        simpl.push_back(it->second);
                    }
                    else{
                        for(int pos=0; pos<3; pos++){
                            simpl.push_back(mesh.get_triangle(it->second).TV(pos));
                        }
                    }
                   //  sort(simpl.begin(),simpl.end(),bind(&Forman_Gradient_SF::sort_filtered_vertices,this,_1,_2));
                        sort_simplex(simpl,mesh);

                    if(simpl == simpl_face)
                        break;
                }
            }
        }
        if(j==(int)simpl_face.size())
        {
            return simpl_face;
        }
    }

    cout<<"[unique_pairable_face] should not print this.. execution paused"<<endl;
    int a; cin>>a;
    return vector<int>();
}


int Forman_Gradient_Computation::num_unpared_faces(vector<int> &co_face, gradient &local_vectors, simplices_map &critical_points, Mesh &mesh)
{
    // give the simplex co_face
    // first extract its faces
    // then, for each of such faces:
    // 1. check if they are paired as "head"
    // 2. check that they are not paired as "tail", searching, as "head" one of its co-faces.

    //NOTICE: a i-simplex can be paired as "head" of a (i-1)-simplex or as "tail" of a (i+1)-simplex

    int num_unpaired=0;
    vector<int> simpl_face = vector<int>();
    vector<int> simpl_face_second = vector<int>();

    for(uint i=1; i<co_face.size(); i++)
    {
        simpl_face = co_face;
        simpl_face.erase(simpl_face.begin()+i);

        if(local_vectors.find(simpl_face) != local_vectors.end(simpl_face))
        {
          //  cerr<<"Case 1.1";
            continue; //if it is already present in the gradient field... nothing to do
        }

        if(critical_points.find(simpl_face) != critical_points.end())
        {
            //cerr<<"Case 1.2";
            continue; //if it is critical... nothing to do
        }
        int j=-1;
       // cout<<"size:"<<simpl_face.size()<<endl;
        if(simpl_face.size() > 1)
        {
            for(j=1; j<(int)simpl_face.size();j++)
            {
                simpl_face_second = simpl_face;

                simpl_face_second.erase(simpl_face_second.begin()+j);

                simplices_map::iterator it = local_vectors.find(simpl_face_second);
                if(it != local_vectors.end(simpl_face_second))
                {
                    //compare two integers.. if I find a identical pair I keep going..
                    vector<int> simpl;

                    if(simpl_face_second.size() == 1){
                        simpl = simpl_face_second;
                        simpl.push_back(it->second);
                    }
                    else{
                        for(int pos=0; pos<3; pos++){
                            simpl.push_back(mesh.get_triangle(it->second).TV(pos));
                        }
                    }
           // sort(simpl.begin(),simpl.end(),bind(&Forman_Gradient_SF::sort_filtered_vertices,this,_1,_2));
                sort_simplex(simpl,mesh);

                    if(simpl == simpl_face)
                        break;
                }
            }
        }
        if(j==(int)simpl_face.size())
        {
            num_unpaired++;
        }
    }
    return num_unpaired;
}

void Forman_Gradient_Computation::sort_simplex(vector<int> &simplex, Mesh &mesh){
        list<int> ordered = list<int>();
       
        for(uint i=0; i<simplex.size(); i++)
        {
            if(i==0)
                ordered.push_back(simplex[i]);
            else
            {
                list<int>::iterator it = ordered.begin();
                for(; it!=ordered.end(); ++it)
                {
                    if(filtration[simplex[i]-1] > filtration[*it-1])
                        break;
                }

                if(it == ordered.begin()) ordered.push_front(simplex[i]);
                else if(it == ordered.end()) ordered.push_back(simplex[i]);
                else ordered.insert(it, simplex[i]);
            }
        }
   
        simplex.clear();
        simplex.assign(ordered.begin(), ordered.end());
   
}
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Roughness.cpp
 * Author: ytsong
 * 
 * Created on January 6, 2019, 2:33 PM
 */

#include "Roughness.h"
#include "timer.h"





void Roughness::compute(Node_V &n, Mesh &mesh, Spatial_Subdivision &division)
 {
    /// if there are no vertices in the leaf we have nothing to do..

      if (n.is_leaf())
    {
        this->roughness_leaf(n,mesh);
        // 
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            if(n.get_son(i)!=NULL)
            {
                this->compute(*n.get_son(i),mesh,division);
            }
        }
       
    }

 }

void Roughness::compute(Node_T &n, Box &n_dom, int level, Mesh &mesh, Spatial_Subdivision &division)
{
    if (n.is_leaf())
    {
        this->roughness_leaf(n,n_dom,mesh);
      //   cout<<"Relation time:"<<time_for_relations<<" Compute time:"<<time_for_compute<<endl;
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(n_dom,level,i);
            int son_level = level +1;
            if(n.get_son(i)!=NULL)
            {
                this->compute(*n.get_son(i),son_dom,son_level,mesh,division);
            }
        }
       
    }
    
   
}

    
void Roughness::roughness_leaf(Node_V& n, Mesh& mesh)
{  
    
    if(!n.indexes_vertices())
        return;

    itype v_start = n.get_v_start();
    itype v_end = n.get_v_end();
    itype v_range = v_end - v_start;
   

    //vvs(v_range,VV_vec());

   
    n.get_VV_vector(vvs,mesh);



    
    for(unsigned i=0;i<v_range;i++)
    {

        itype real_v_id=v_start+i;
      
        VV_vec &vv = vvs[i];
        coord_type zDistSum = 0;
        coord_type zSum=0;
        Vertex &v = mesh.get_vertex(real_v_id);
        zSum+=2*v.get_z();
           for(auto vid:vv)
        { 
            Vertex &v2=mesh.get_vertex(vid);
            zSum += v2.get_z();
           
        }
        coord_type zAVG=zSum/(vv.size()+2);
        zDistSum +=2* (v.get_z()-zAVG) * (v.get_z()-zAVG);
        for(auto vid:vv)
        {
            Vertex &v2=mesh.get_vertex(vid);
            zDistSum+= (v2.get_z()-zAVG)*(v2.get_z()-zAVG);
                    
        }
        
        //roughness[real_v_id]
        coord_type roughness=sqrt(zDistSum/(vv.size()+2));
        mesh.get_vertex(i+v_start).add_field(roughness);
    }
    vvs.clear();

 
}

void Roughness::roughness_leaf(Node_T &n, Box &n_dom, Mesh &mesh)
{
    itype v_start;
    itype v_end;

    n.get_v_range(v_start,v_end,n_dom,mesh); // we need to gather the vertices range..
    
    if(v_start == v_end) //no internal vertices..
        return;
   
   // leaf_VV_vec vvs(v_end-v_start,VV_vec());
    //dvect v_aux(v_range,0.0);
 
       //    
    n.get_VV_vector(vvs,v_start,v_end,mesh);

    for(unsigned i=0;i<vvs.size();i++)
    {
        itype real_v_id=v_start+i;
       
        VV_vec &vv = vvs[i];
        coord_type zDistSum = 0.0;
        coord_type zSum=0.0;
        Vertex &v = mesh.get_vertex(real_v_id);
        zSum+=2*v.get_z();
        for(auto vid:vv)
        {
            Vertex &v2=mesh.get_vertex(vid);
            zSum += v2.get_z();
        
        }
        coord_type zAVG=zSum/(vv.size()+2);
        zDistSum += 2*(v.get_z()-zAVG) * (v.get_z()-zAVG);
        for(auto vid:vv)
        {
            Vertex &v2=mesh.get_vertex(vid);
            zDistSum+= (v2.get_z()-zAVG)*(v2.get_z()-zAVG);
                    
        }
        coord_type roughness=sqrt(zDistSum/vv.size()+2);
        mesh.get_vertex(i+v_start).add_field(roughness);
    }

    vvs.clear();

}


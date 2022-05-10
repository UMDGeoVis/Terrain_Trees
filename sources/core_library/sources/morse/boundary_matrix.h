#ifndef BOUNDARY_MATRIX_H
#define BOUNDARY_MATRIX_H

#include <sys/time.h>
#include <sys/resource.h>

#include <list>
#include <sstream>
#include <queue>
#include <cstring>

#include "basic_types/mesh.h"
#include "utilities/timer.h"

using namespace std;

template <class dataType>
class Boundary_Matrix
{
public:
    inline void getMatrix(vector<vector<int> >& origMatrix){
            origMatrix = vector<vector<int> >(matrix.size(), vector<int>());
            for(int i=0; i<matrix.size(); i++){
                vector<int> vec;
                for(auto v : matrix[i]){
                    vec.push_back(filtration[v]);
                }
                origMatrix[filtration[i]] = vec;
            }
        }

    inline void addValue(int simplex, int boundary, dataType filtr){
        
        if(matrix.size() <= simplex){
            for(int i=matrix.size(); i<=simplex; i++){
                matrix.push_back(vector<int>());
                values.push_back(0);
            }
        }
        values[simplex] = filtr;

        if(boundary > -1){
            auto lb = std::lower_bound(matrix[simplex].begin(), matrix[simplex].end(), boundary);
            if (lb != matrix[simplex].end() && *lb == boundary) {
                auto ub = std::upper_bound(lb, matrix[simplex].end(), boundary);
                matrix[simplex].erase(lb, ub);
            }
            else{
                matrix[simplex].insert(lb,boundary);
            }
        }
    }

    inline void sort(){

        //organize the simplices according to the filtration value
        map<int, vector<int> > organize_values;
        for(int i=0; i<values.size(); i++){

            if(organize_values.find(values[i]) == organize_values.end())
                organize_values[values[i]] = vector<int>();
            organize_values[values[i]].push_back(i);
            
        }

        //compute the filtration value they should have in the matrix
        filtration = vector<int>(values.size(), -1);
        int fI = 0;
        for(auto levelset : organize_values){
            vector<int> simplices = levelset.second;
            //here i should sort these guys

            for(int i=0; i < simplices.size(); i++){
                int s = simplices[i];

                if(matrix[s].size() == 0){
                    filtration[s]=fI++;
                }
                else{
                    if( boundary_done(s, filtration)){
                        filtration[s]=fI++;
                    }
                    else{
                        int j=0;
                        for(j=i+1; j<simplices.size(); j++){

                            if(boundary_done(simplices[j],filtration)){
                                simplices[i] = simplices[j];
                                simplices[j] = s;
                                break;
                            }

                            
                        }
                        i--;
                    }
                }
            }
        }
        organize_values.clear();

        //sort the structure accordingly
        vector<int> newfiltration = filtration;
        for(int i=0; i<matrix.size(); i++){
            if(i == newfiltration[i]){
                update_indices(matrix[i], filtration);
            }
            else{
                vector<int> temp1 = matrix[i];
                matrix[i] = matrix[newfiltration[i]];
                matrix[newfiltration[i]] = temp1;
                
                dataType val = values[i];
                values[i] = values[newfiltration[i]];
                values[newfiltration[i]] = val;

                int index = newfiltration[newfiltration[i]];
                newfiltration[newfiltration[i]] = newfiltration[i];
                newfiltration[i] = index;
                i--;
            }
        }

        for(int i=0; i<filtration.size(); i++){
            newfiltration[filtration[i]] = i;
        }
        filtration=newfiltration;

    }

    inline void print(){

        int index=0;    
        for(auto vect : matrix){
            cout << values[index] << "| " << index++ << ": ";
            for(auto v : vect){
                cout << v << " ";
            }
            cout << endl;
        }
    }

    inline void reduce(){
        cout << "Reducing a matrix with " << matrix.size() << " columns" << endl;
        allpivots = vector<int> (matrix.size(), -1);

        for(int i=0; i<matrix.size(); i++){
            
            if(pivot(i) != -1){
                
                while(allpivots[pivot(i)] != -1){
                    addcolumn(i,allpivots[pivot(i)]);
                    if(pivot(i) == -1)
                        break;
                }

                if(matrix.size() != 0 && pivot(i) != -1)
                    allpivots[pivot(i)]=i;
            } 
        }

    }

    inline void getPairs(map<int,int>& pairs){

        pairs = map<int,int>();
        for(int i=0; i<allpivots.size(); i++){
            if(allpivots[i] != -1){
                    pairs[filtration[i]]=filtration[allpivots[i]];
            }
        }
    }

    inline void getHomology(vector<int>& ret){

        // pivots = vector<int>(allpivots.size(), -1);
        ret = vector<int>();

        set<int> pivots;
        for(int i=0; i<allpivots.size(); i++){
            pivots.insert(allpivots[i]);
        }


        for(int i=0; i<allpivots.size(); i++){
            if(allpivots[i] == -1 && pivots.find(i) == pivots.end())
                ret.push_back(filtration[i]);
        }
    }

private:

    vector<dataType> values; //contain the filtration value as it was passed
    vector<int> filtration; //contain the index to the original simplex as it was passed
    
    vector<vector<int> > matrix;
    vector<int> allpivots;


    inline bool boundary_done(int simpl, vector<int>& filtration){
        for(auto bd : matrix[simpl]){
            if(filtration[bd] == -1)
                return false;
        }
        return true;
    }

    inline void update_indices(vector<int>& vec, vector<int>& filtration){
        for(int i=0; i<vec.size(); i++)
            vec[i] = filtration[vec[i]];
        std::sort(vec.begin(), vec.end());
    }

    inline int pivot(int i){
        if(matrix[i].size() == 0)
            return -1;
        return matrix[i].back();
    }

    inline void addcolumn(int i, int j){

        for(auto v : matrix[j]){

            auto pr = equal_range(matrix[i].begin(), matrix[i].end(), v);

            if(pr.first == pr.second){
                matrix[i].insert(pr.second,v);
            }
            else{
                matrix[i].erase(pr.first, pr.second);
            }
            
        }
    }

};

#endif
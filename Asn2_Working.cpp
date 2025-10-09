#include <bits/stdc++.h>
using namespace std;

int main(){
    vector<vector<int>> dna_coords = {
    {0,3,2}, {1,2,3}, {2,1,4}, {1,0,5}, {0,1,4},
    {1,2,5}, {2,3,4}, {3,4,5}, {2,5,4}, {1,4,3},
    {0,5,2}, {1,4,1}, {2,5,0}, {3,4,1}, {4,3,2},
    {3,4,3}, {4,5,4}, {5,4,5}, {4,3,4}, {3,2,5}
};
    vector<vector<int>> pro_coords = {
    {4,5,2}, {5,4,1}, {4,3,2}, {3,4,3}, {4,5,4},
    {3,4,3}, {4,5,4}, {2,3,4}, {1,2,5}, {2,3,4},
    {1,2,5}, {1,2,5}, {0,1,4}, {1,2,5}, {2,3,4},
    {2,1,4}, {1,2,3}, {2,1,4}, {1,2,3}, {0,3,2},
    {1,2,1}, {0,3,0}, {1,2,1}, {2,1,0}, {1,2,1},
    {0,3,0}, {0,5,2}
};
    set<vector<int>> st;
    for(auto coords : dna_coords)st.insert(coords);
    int prev = -1, cur = -1;
    int steps_1d = 0, steps_3d = 0, steps_dis = 0, steps_ass = 0;
    for(auto coords : pro_coords){
        if(st.find(coords) == st.end()){
            cur = 0;
        }else{
            cur = 1;
        }
        if(prev != -1){
            if(cur == 1){
                if(prev == 1){
                    steps_1d++;
                }else{
                    steps_ass++;
                }
            }else{
                if(prev == 1){
                    steps_dis++;
                }else{
                    steps_3d++;
                }   
            }
        }
        swap(prev,cur);
    }

    cout << "1D steps : " << steps_1d << "\n" << "3D steps : " << steps_3d << "\n" << "Association steps : " << steps_ass << "\n" << "Dissociation steps : " << steps_dis << "\n";

    return 0;
}

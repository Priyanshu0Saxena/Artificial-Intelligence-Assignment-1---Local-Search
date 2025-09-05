#include "solver.h"
#include <iostream>
#include <chrono>
#include <random>

using namespace std;

struct State{
    vector<vector<int>> helicopterVillage;
    State(){}
    State(int n): helicopterVillage(n){};
};

State Random_State(const ProblemData& problem){
    random_device rd;
    mt19937 rng(rd());
    int H=problem.helicopters.size();
    int V=problem.villages.size();
    State state;
    state.helicopterVillage.resize(H);
    int target_Helicopter;
    uniform_int_distribution<int>dist(0,H-1);
    for(int i=0;i<V;i++){
        target_Helicopter=dist(rng);
        state.helicopterVillage[target_Helicopter].push_back(problem.villages[i].id);
    }
    return state;
}

Solution solve(const ProblemData& problem) {
    cout << "Starting solver..." << endl;

    Solution solution;
    

    cout << "Solver finished." << endl;
    return solution;
}

#include "solver.h"
#include <iostream>
#include <chrono>

using namespace std;

struct State{
    vector<vector<Point>> helicopterVillage;
    State(){}
    State(int n): helicopterVillage(n){};
};

Solution solve(const ProblemData& problem) {
    cout << "Starting solver..." << endl;

    Solution solution;
    

    cout << "Solver finished." << endl;
    return solution;
}

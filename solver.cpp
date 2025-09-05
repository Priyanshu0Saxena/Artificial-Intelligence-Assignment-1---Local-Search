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

double calculateTourDistance(int hid, const vector<int>& villages, ProblemData& problem)
{
    if(villages.size()==0)
    return 0.0;

    const auto& helicopter = problem.helicopters[hid];
    Point home = problem.cities[helicopter.home_city_id - 1];
    Point curr = home;
    double totalDistance = 0.0;

    for(int vill: villages)
    {
        Point des = problem.villages[vill].coords;
        totalDistance += distance(curr, des);
        curr = des;
    }
    totalDistance += distance(curr, home);
    return totalDistance;
}
double calculateVillageValue(int village_id, int dry_food, int perishable_food, int other_supplies, double current_food_delivered, double current_other_delivered)
{

}


double calculateObjectiveValue(const Solution& solution, ProblemData& problem)
{
    double totalValue, totalCost;
    totalCost=totalValue=0;

    int H = (int)problem.helicopters.size();
    int V = (int)problem.villages.size();

    vector<double> foodDelivered(V+1, 0.0);
    vector<double> otherDelivered(V+1, 0.0);
    vector<double> helicopterTotalDistance(H, 0.0);

    for(int h = 0;h<H; h++)
    {
        const auto& plan = solution[h];
        const auto& helicopter = problem.helicopters[plan.helicopter_id-1];

        if(plan.trips.empty())
        continue;

        totalCost += helicopter.fixed_cost;

        for(const auto& trip: plan.trips)
        {
            if(trip.drops.empty())
            continue;

            double tripWeight = trip.dry_food_pickup*problem.packages[0].weight;
            tripWeight += trip.perishable_food_pickup*problem.packages[1].weight;
            tripWeight += trip.other_supplies_pickup*problem.packages[2].weight;

            if(tripWeight > helicopter.weight_capacity)
            return -1e9;
            int totalDryDrops = 0, totalPerishableDrops = 0, totalOtherDrops = 0;

            for(const auto& drop : trip.drops)
            {
                totalDryDrops += drop.dry_food;
                totalOtherDrops += drop.other_supplies;
                totalPerishableDrops += drop.perishable_food;
            }

            if (totalDryDrops > trip.dry_food_pickup ||
                totalPerishableDrops > trip.perishable_food_pickup ||
                totalOtherDrops > trip.other_supplies_pickup)
                return -1e9;

            vector<int> route;
            for(const auto& drop: trip.drops)
            route.push_back(drop.village_id);

            double tripDistance = calculateTourDistance(h, route, problem);
            if(tripDistance > helicopter.distance_capacity)
            return -1e9;

            if(helicopterTotalDistance[h] > problem.d_max - tripDistance)
            return -1e9;

            helicopterTotalDistance[h] += tripDistance;

            for(const auto& drop: trip.drops)
            {
                totalValue += calculateVillageValue(drop.village_id, drop.dry_food, drop.perishable_food, drop.other_supplies, foodDelivered[drop.village_id], otherDelivered[drop.village_id]);
                foodDelivered[drop.village_id] += drop.dry_food + drop.perishable_food;
                otherDelivered[drop.village_id] += drop.other_supplies;
            }
            
            totalCost += helicopter.alpha*tripDistance;
        }

    }
    return totalValue - totalCost;
}

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

State neighbor(State state)
{
    int source = -1;
    int destination = -1;
    int h = state.helicopterVillage.size();
    if (h <= 1)
    {
        return state;
    }
    int limit;
    random_device rd;
    mt19937 rng(rd());
    int choice = uniform_int_distribution<int>(1, 3)(rng);
    uniform_int_distribution<int> dist(0, h - 1);
    source=dist(rng);
    destination=dist(rng);
    switch (choice)
    {
    case 1:
        cout<<"Case 1"<<endl;
        while (source == destination)
        {
            source = dist(rng);
            destination = dist(rng);
        }
        swap(state.helicopterVillage[source], state.helicopterVillage[destination]);
        break;
    case 2:
        cout<<"case 2 "<<endl;
        limit=100;
        while(state.helicopterVillage[source].size()<=1 && limit){
            source=dist(rng);
            limit--;
        }
        if(state.helicopterVillage[source].size()>1)
        shuffle(state.helicopterVillage[source].begin(),state.helicopterVillage[source].end(),rng);
        break;
    case 3:
    cout<<"case 3"<<endl;
    limit =h+10;
        while((state.helicopterVillage[source].size()==0 || source==destination) && limit){
            source=dist(rng);
            destination=dist(rng);
            limit--;
        }
        if(state.helicopterVillage[source].size()==0){
            return state;
        }
        int randomvillag;
        uniform_int_distribution<int> dist1(0,state.helicopterVillage[source].size()-1);
        randomvillag=dist1(rng);
        state.helicopterVillage[destination].push_back(state.helicopterVillage[source][randomvillag]);
        state.helicopterVillage[source].erase(state.helicopterVillage[source].begin()+randomvillag);
        break;
    }
    return state; 
}

Solution solve(const ProblemData& problem) {
    cout << "Starting solver..." << endl;

    Solution solution;
    

    cout << "Solver finished." << endl;
    return solution;
}

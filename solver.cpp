#include "solver.h"
#include <iostream>
#include <random>
#include <algorithm>
using namespace std;

struct State
{
    vector<vector<int>> helicopterVillage;
    State(int n) : helicopterVillage(n){};
};

State randomState(const ProblemData& problem)
{
    random_device rd;
    mt19937 rng(rd());
    int H = problem.helicopters.size();
    int V = problem.villages.size();
    State state(H);

    for (int i = 0; i < V; i++) {
        int target_Helicopter = rand()%H;
        state.helicopterVillage[target_Helicopter].push_back(problem.villages[i].id);
    }
    return state;
}

State getNeighbour(State state)
{
    int h = state.helicopterVillage.size();
    if (h <= 1) return state;
    random_device rd;
    mt19937 rng(rd());
    int choice = uniform_int_distribution<int>(1, 3)(rng);
    uniform_int_distribution<int> dist(0, h - 1);
    int src = dist(rng);
    int dest = dist(rng);
    int limit;
    switch (choice)
    {
        case 1:
            while (src == dest)
            {
                src = dist(rng);
                dest = dist(rng);
            }
            swap(state.helicopterVillage[src], state.helicopterVillage[dest]);
            break;
        case 2:
            limit = 100;
            while (state.helicopterVillage[src].size() <= 1 && limit)
            {
                src = dist(rng);
                limit--;
            }
            if (state.helicopterVillage[src].size() > 1)
                shuffle(state.helicopterVillage[src].begin(), state.helicopterVillage[src].end(), rng);
            break;
        case 3:
        {
            limit = h + 10;
            while ((state.helicopterVillage[src].size() == 0 || src == dest) && limit)
            {
                src = dist(rng);
                dest = dist(rng);
                limit--;
            }
            if (state.helicopterVillage[src].size() == 0) return state;
            uniform_int_distribution<int> dist1(0, state.helicopterVillage[src].size() - 1);
            int randomVillage = dist1(rng);
            state.helicopterVillage[dest].push_back(state.helicopterVillage[src][randomVillage]);
            state.helicopterVillage[src].erase(state.helicopterVillage[src].begin() + randomVillage);
            break;
        }
    }
    return state;
}

int pickup(double dryOther, double wcap)
{
    int dryOtherCount = wcap / dryOther;
    random_device rd;
    mt19937 rng(rd());
    return uniform_int_distribution<int>(0, dryOtherCount)(rng);
}

double getTripDistance(int hid, const vector<int>& villages, ProblemData& problem)
{
    if (villages.empty())
    return 0.0;

    const auto& helicopter = problem.helicopters[hid];
    Point homeCity = problem.cities[helicopter.home_city_id - 1];

    if (villages.size() == 1)
    {
        Point village_pos = problem.villages[villages[0] - 1].coords;
        return 2.0 * distance(homeCity, village_pos);
    }
    vector<Point> points;
    points.push_back(homeCity);

    for (int vid : villages)
    points.push_back(problem.villages[vid - 1].coords);

    int n = (int)points.size();
    vector<bool> visited(n, 0);
    visited[0] = true;

    int current = 0;
    double totalDis = 0.0;

    for (int i = 1; i < n; i++)
    {
        double bestDis = 1e9;
        int next = -1;
        for (int j = 1; j < n; j++)
        {
            if (!visited[j])
            {
                double d = distance(points[current], points[j]);
                if (d < bestDis) {
                    bestDis = d;
                    next = j;
                }
            }
        }
        if (next != -1) {
            visited[next] = true;
            totalDis += bestDis;
            current = next;
        }
    }
    totalDis += distance(points[current], homeCity);
    return totalDis;
}

Solution getSolution(const State& state, const ProblemData& problem)
{
    const double EPS = 1e-9;
    int V = (int)problem.villages.size();
    int H = (int)problem.helicopters.size();
    vector<int> peopleRem(V + 1, 0);

    for (const auto &v : problem.villages)
    peopleRem[v.id] = v.population;

    double dryComboWeight = 9.0 * problem.packages[0].weight + problem.packages[2].weight;
    double wetComboWeight = 9.0 * problem.packages[1].weight + problem.packages[2].weight;

    Solution solution;
    solution.resize(H);

    for (int h = 0; h < H; ++h)
    {
        const auto &heli = problem.helicopters[h];
        Point home = problem.cities[heli.home_city_id - 1];
        double heliTotalDistance = 0.0;
        HelicopterPlan plan;
        plan.helicopter_id = heli.id;
        const vector<int> &assigned = state.helicopterVillage[h];
        unsigned int index = 0;
        while (index < assigned.size())
        {
            int vid = assigned[index];
            if (vid < 1 || vid > V) { ++index; continue; }
            if (peopleRem[vid] <= 0) { ++index; continue; }

            int dryComboTaken = pickup(dryComboWeight, heli.weight_capacity);
            double remainingWeightCap = heli.weight_capacity - dryComboTaken * dryComboWeight;
            int wetComboTaken = 0;

            if (remainingWeightCap > EPS)
            wetComboTaken = (int)floor(remainingWeightCap / wetComboWeight);

            int dryMealTaken = dryComboTaken * 9;
            int perishableMealTaken = wetComboTaken * 9;
            int otherTaken = dryComboTaken + wetComboTaken;

            Trip trip;
            trip.dry_food_pickup = dryMealTaken;
            trip.perishable_food_pickup = perishableMealTaken;
            trip.other_supplies_pickup = otherTaken;
            Point currentPos = home;

            double tripDistance = 0.0;
            bool haveDrops = false;

            while (index < assigned.size())
            {
                int vid2 = assigned[index];
                if (vid2 < 1 || vid2 > V) { ++index; continue; }
                if (peopleRem[vid2] <= 0) { ++index; continue; }
                Point villagePos = problem.villages[vid2 - 1].coords;
                double distanceToVillage = distance(currentPos, villagePos);
                double distanceVillageToHome = distance(villagePos, home);
                double distanceIfVisit = tripDistance + distanceToVillage + distanceVillageToHome;
                if (distanceIfVisit > heli.distance_capacity + EPS)
                {
                    if (!haveDrops)
                    {
                        double roundtrip_from_home = distance(home, villagePos) + distance(villagePos, home);
                        if (roundtrip_from_home > heli.distance_capacity + EPS) { ++index; continue; }
                        else break;
                    }
                    else 
                    {
                        tripDistance += distance(currentPos, home);
                        heliTotalDistance += tripDistance;
                        if (heliTotalDistance > problem.d_max + EPS) goto changeHeli;
                        if (!trip.drops.empty()) plan.trips.push_back(trip);
                        break;
                    }
                }

                int totalFoodInTrip = dryMealTaken + perishableMealTaken;
                int numberOfPeopleSatisfiedUsingOnlyFood = totalFoodInTrip / 9;
                int numberOfPeopleSatisfiedUsingOnlyOthers = otherTaken;
                int maxPeopleSatisfiedUsingBoth = min(numberOfPeopleSatisfiedUsingOnlyFood, numberOfPeopleSatisfiedUsingOnlyOthers);
                if (maxPeopleSatisfiedUsingBoth <= 0)
                {
                    if (!haveDrops)
                    break;
                    else
                    {
                        tripDistance += distance(currentPos, home);
                        heliTotalDistance += tripDistance;
                        if (heliTotalDistance > problem.d_max + EPS) goto changeHeli;
                        if (!trip.drops.empty()) plan.trips.push_back(trip);
                        break;
                    }
                }
                int people = peopleRem[vid2];
                int serveCount = min(people, maxPeopleSatisfiedUsingBoth);
                int foodReq = serveCount * 9;
                int PerishableUse = min(perishableMealTaken, foodReq);
                foodReq -= PerishableUse;
                int dryUse = min(dryMealTaken, foodReq);
                foodReq -= dryUse;

                if (foodReq > 0)
                {
                    serveCount -= (int)ceil((double)foodReq / 9.0);
                    if (serveCount <= 0) break;
                    foodReq = serveCount * 9;
                    PerishableUse = min(perishableMealTaken, foodReq);
                    foodReq -= PerishableUse;
                    dryUse = min(dryMealTaken, foodReq);
                    foodReq -= dryUse;
                }

                int otherUse = min(otherTaken, serveCount);
                if (otherUse < serveCount)
                {
                    int newServe = otherUse;
                    if (newServe <= 0) {
                        if (!haveDrops) break;
                        else {
                            tripDistance += distance(currentPos, home);
                            heliTotalDistance += tripDistance;
                            if (heliTotalDistance > problem.d_max + EPS) goto changeHeli;
                            if (!trip.drops.empty()) plan.trips.push_back(trip);
                            break;
                        }
                    }
                    serveCount = newServe;
                    foodReq = serveCount * 9;
                    PerishableUse = min(perishableMealTaken, foodReq);
                    foodReq -= PerishableUse;
                    dryUse = min(dryMealTaken, foodReq);
                    foodReq -= dryUse;
                }

                Drop drop;
                drop.village_id = vid2;
                drop.perishable_food = PerishableUse;
                drop.dry_food = dryUse;
                drop.other_supplies = serveCount;
                trip.drops.push_back(drop);
                haveDrops = true;
                perishableMealTaken -= PerishableUse;
                dryMealTaken -= dryUse;
                otherTaken -= serveCount;
                peopleRem[vid2] -= serveCount;
                tripDistance += distanceToVillage;
                currentPos = villagePos;

                if (peopleRem[vid2] > 0)
                {
                    tripDistance += distanceVillageToHome;
                    heliTotalDistance += tripDistance;
                    if (heliTotalDistance > problem.d_max + EPS) goto changeHeli;
                    if (!trip.drops.empty()) plan.trips.push_back(trip);
                    break;
                }
                else
                ++index;
            }

            if (haveDrops && (trip.drops.size() > 0))
            {
                if (!(currentPos.x == home.x && currentPos.y == home.y))
                {
                    tripDistance += distance(currentPos, home);
                    heliTotalDistance += tripDistance;
                    if (heliTotalDistance <= problem.d_max + EPS) plan.trips.push_back(trip);
                }
                else
                    plan.trips.push_back(trip);
            }
            if (heliTotalDistance > problem.d_max + EPS)
            break;
        }
        changeHeli:
        solution[h] = plan;
    }
    return solution;
}


double calculateVillageValue(int vid, int dryFood, int perishableFood, int others, double foodToAdd, double othersToAdd, ProblemData& problem)
{
    const auto& village = problem.villages[vid - 1];
    double maxFoodReq = village.population * 9.0;
    double maxOtherReq = village.population * 1.0;

    double foodLeft = max(0.0, maxFoodReq - foodToAdd);
    double foodInDrop = dryFood + perishableFood;

    double finalFoodInDrop = min(foodInDrop, foodLeft);
    double effPerishable = min((double)perishableFood, finalFoodInDrop);
    double remFinalFoodInDrop = finalFoodInDrop - effPerishable;
    double effDry = min((double)dryFood, remFinalFoodInDrop);
    double otherLeft = max(0.0, maxOtherReq - othersToAdd);
    double effOther = min((double)others, otherLeft);
    return effPerishable * problem.packages[1].value + effDry * problem.packages[0].value + effOther * problem.packages[2].value;
}


double getObjVal(Solution solution, ProblemData problem)
{
    double totalValue = 0.0;
    double totalCost = 0.0;
    int H = (int)problem.helicopters.size();
    int V = (int)problem.villages.size();

    vector<double> foodDelivered(V + 1, 0.0);
    vector<double> otherDelivered(V + 1, 0.0);
    vector<double> heliTotalDistance(H, 0.0);

    for (int h = 0; h < H; ++h) {
        const auto& plan = solution[h];
        const auto& helicopter = problem.helicopters[plan.helicopter_id-1];
        if (plan.trips.empty()) continue;
        totalCost += helicopter.fixed_cost;
        for (const auto& trip : plan.trips)
        {
            if (trip.drops.empty()) continue;
            double tripWeight = trip.dry_food_pickup * problem.packages[0].weight + trip.perishable_food_pickup * problem.packages[1].weight + trip.other_supplies_pickup * problem.packages[2].weight;
            if (tripWeight > helicopter.weight_capacity + 1e-6) return -1e9;
            int totalDryDrops = 0, totalPerishableDrops = 0, totalOtherDrops = 0;
            for (const auto& drop : trip.drops) {
                totalDryDrops += drop.dry_food;
                totalPerishableDrops += drop.perishable_food;
                totalOtherDrops += drop.other_supplies;
            }
            if (totalDryDrops > trip.dry_food_pickup || totalPerishableDrops > trip.perishable_food_pickup || totalOtherDrops > trip.other_supplies_pickup)
            return -1e9;

            vector<int> path;
            path.reserve(trip.drops.size());

            for (const auto& drop : trip.drops)
            path.push_back(drop.village_id);

            double tripDistance = getTripDistance(plan.helicopter_id - 1, path, problem);

            if (tripDistance > helicopter.distance_capacity + 1e-6)
            return -1e9;

            heliTotalDistance[h] += tripDistance;
            
            if (heliTotalDistance[h] > problem.d_max + 1e-6)
            return -1e9;

            for (const auto& drop : trip.drops) {
                totalValue += calculateVillageValue(drop.village_id, drop.dry_food, drop.perishable_food, drop.other_supplies, foodDelivered[drop.village_id], otherDelivered[drop.village_id], problem);
                foodDelivered[drop.village_id] += drop.dry_food + drop.perishable_food;
                otherDelivered[drop.village_id] += drop.other_supplies;
            }
            totalCost += helicopter.alpha * tripDistance;
        }
    }
    return totalValue - totalCost;
}

State tweakIt(const State& current, const ProblemData& problem) {
    State neighbour = current;
    random_device rd;
    mt19937 rng(rd());
    int choice = uniform_int_distribution<int>(1, 2)(rng);
    switch (choice)
    {
        case 1:
            neighbour = randomState(problem);
            break;
        case 2:
            neighbour = getNeighbour(neighbour);
            break;
    }
    return neighbour;
}

Solution helper(const ProblemData& problem)
{
    Solution bestSol;
    double bestScore = 0;
    auto currentTime = chrono::steady_clock::now();
    auto timeLimit = (long long)(problem.time_limit_minutes*60.0*1000*0.95);

    int restartCount = 0;
    while(chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - currentTime).count() < timeLimit)
    {
        restartCount++;
        State current = randomState(problem);
        Solution currentSolution = getSolution(current, problem);
        double currentObjValue = getObjVal(currentSolution, problem);
        int noImprovCount = 0, maxItrCount = 1000;
        while (noImprovCount < maxItrCount && chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - currentTime).count() < timeLimit)
        {
            State neighbour = tweakIt(current, problem);
            Solution neighbourSolution = getSolution(neighbour, problem);
            double neighbourObjValue = getObjVal(neighbourSolution, problem);
            if (neighbourObjValue > currentObjValue) {
                current = move(neighbour);
                currentSolution = move(neighbourSolution);
                currentObjValue = neighbourObjValue;
                noImprovCount = 0;
            }
            else
            noImprovCount++;
        }
        if (currentObjValue > bestScore)
        {
            bestScore = currentObjValue;
            bestSol = currentSolution;
        }
    }
    return bestSol;
}

Solution solve(const ProblemData& problem)
{
    Solution solution = helper(problem);
    return solution;
}

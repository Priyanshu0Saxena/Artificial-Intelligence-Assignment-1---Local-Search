#include "solver.h"
#include <iostream>
#include <chrono>
#include <random>
#include <algorithm>
#include <climits>
#include <cmath>
#include <map>
#include <vector>
#include <set>
#include<utility>

using namespace std;

struct State {
    vector<vector<int>> helicopter_villages;
    State(int num_helicopters = 0) : helicopter_villages(num_helicopters) {}
};

class LocalSearchSolver {
private:
    const ProblemData& problem;
    mt19937 rng;

    // Calculate TSP-like tour distance for multiple villages
    double calculateTourDistance(int heli_idx, const vector<int>& villages) const {
        if (villages.empty()) return 0.0;
        
        const auto& helicopter = problem.helicopters[heli_idx];
        Point home = problem.cities[helicopter.home_city_id - 1];
        
        if (villages.size() == 1) {
            Point village_pos = problem.villages[villages[0] - 1].coords;
            return 2.0 * distance(home, village_pos);
        }
        
        // For multiple villages, find a reasonable tour
        vector<Point> points;
        points.push_back(home);
        for (int vid : villages) {
            points.push_back(problem.villages[vid - 1].coords);
        }
        
        // Simple nearest neighbor TSP approximation
        vector<bool> visited(points.size(), false);
        visited[0] = true;
        int current = 0;
        double total_dist = 0.0;
        
        for (int i = 1; i < points.size(); ++i) {
            double min_dist = 1e9;
            int next = -1;
            for (int j = 1; j < points.size(); ++j) {
                if (!visited[j]) {
                    double d = distance(points[current], points[j]);
                    if (d < min_dist) {
                        min_dist = d;
                        next = j;
                    }
                }
            }
            if (next != -1) {
                visited[next] = true;
                total_dist += min_dist;
                current = next;
            }
        }
        
        // Return to home
        total_dist += distance(points[current], home);
        return total_dist;
    }

    double calculateVillageValue(int village_id, int dry_food, int perishable_food, int other_supplies,
                                 double current_food_delivered, double current_other_delivered) const {
        const auto& village = problem.villages[village_id - 1];
        double max_food_needed = village.population * 9.0;
        double max_other_needed = village.population * 1.0;

        double food_room_left = max(0.0, max_food_needed - current_food_delivered);
        double food_in_this_drop = dry_food + perishable_food;
        double effective_food_this_drop = min(food_in_this_drop, food_room_left);

        double effective_perishable = min((double)perishable_food, effective_food_this_drop);
        double value_from_perishable = effective_perishable * problem.packages[1].value;

        double remaining_effective_food = effective_food_this_drop - effective_perishable;
        double effective_dry = min((double)dry_food, remaining_effective_food);
        double value_from_dry = effective_dry * problem.packages[0].value;

        double other_room_left = max(0.0, max_other_needed - current_other_delivered);
        double effective_other = min((double)other_supplies, other_room_left);
        double value_from_other = effective_other * problem.packages[2].value;

        return value_from_perishable + value_from_dry + value_from_other;
    }


State Random_State(const ProblemData& problem){
    random_device rd;
    mt19937 rng(rd());
    int H=problem.helicopters.size();
    int V=problem.villages.size();
    State state;
    state.helicopter_villages.resize(H);
    int target_Helicopter;
    uniform_int_distribution<int>dist(0,H-1);
    for(int i=0;i<V;i++){
        target_Helicopter=dist(rng);
        state.helicopter_villages[target_Helicopter].push_back(problem.villages[i].id);
    }
    return state;
}

State neighborGen(State state)
{
    int source = -1;
    int destination = -1;
    int h = state.helicopter_villages.size();
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
        while (source == destination)
        {
            source = dist(rng);
            destination = dist(rng);
        }
        swap(state.helicopter_villages[source], state.helicopter_villages[destination]);
        break;
    case 2:
        limit=100;
        while(state.helicopter_villages[source].size()<=1 && limit){
            source=dist(rng);
            limit--;
        }
        if(state.helicopter_villages[source].size()>1)
        shuffle(state.helicopter_villages[source].begin(),state.helicopter_villages[source].end(),rng);
        break;
    case 3:
    limit =h+10;
        while((state.helicopter_villages[source].size()==0 || source==destination) && limit){
            source=dist(rng);
            destination=dist(rng);
            limit--;
        }
        if(state.helicopter_villages[source].size()==0){
            return state;
        }
        int randomvillag;
        uniform_int_distribution<int> dist1(0,state.helicopter_villages[source].size()-1);
        randomvillag=dist1(rng);
        state.helicopter_villages[destination].push_back(state.helicopter_villages[source][randomvillag]);
        state.helicopter_villages[source].erase(state.helicopter_villages[source].begin()+randomvillag);
        break;
    }
    return state; 
}
    State perturbState(const State& current) {
        State neighbor = current;
        int operation = uniform_int_distribution<int>(1, 2)(rng);
        
        switch(operation) {
            case 1: neighbor = Random_State(problem); break;
            case 2: neighbor = neighborGen(neighbor); break;
        }
        
        return neighbor;
    }

public:
    LocalSearchSolver(const ProblemData& problem_data)
        : problem(problem_data),
          rng((uint32_t)chrono::steady_clock::now().time_since_epoch().count()) {}
int randint(int i, int j)
{
    random_device rd;
    mt19937 rng(rd());
    return uniform_int_distribution<int>(i, j)(rng);
}

int pickup(double dry_other, double wcap)
{
    int dry_other_count = wcap / dry_other;
    return randint(0, dry_other_count);
}

    Solution stateToSolution(const State& state, const ProblemData& problem) {
    const double EPS = 1e-9;
    int numVillages = (int)problem.villages.size();
    int numHelis = (int)problem.helicopters.size();

    // remaining people to serve per village (1-based indexing)
    vector<int> people_remaining(numVillages + 1, 0);
    for (const auto &v : problem.villages) {
        people_remaining[v.id] = v.population;
    }

    // combo weights: one combo = 9 food packets + 1 other supply
    double dry_combo_weight = 9.0 * problem.packages[0].weight + problem.packages[2].weight;
    double wet_combo_weight  = 9.0 * problem.packages[1].weight + problem.packages[2].weight;

    Solution solution;
    solution.resize(numHelis);

    // Process each helicopter independently (state.helicopterVillage[h] lists village ids assigned to h)
    for (int h = 0; h < numHelis; ++h) {
        const auto &heli = problem.helicopters[h];
        Point home = problem.cities[heli.home_city_id - 1]; // convert to 0-based
        double heli_total_distance = 0.0;                    // DMax across trips
        HelicopterPlan plan;
        plan.helicopter_id = heli.id;

        // Iterate assigned villages in order, serving one village fully before next
        const vector<int> &assigned = state.helicopter_villages[h];
        size_t idx = 0;

    heli_outer_loop:
        while (idx < assigned.size()) {
            int vid = assigned[idx];
            if (vid < 1 || vid > numVillages) { ++idx; continue; } // defensive

            if (people_remaining[vid] <= 0) { ++idx; continue; } // already served

            // Start a new trip (load random combos according to your rule)
            int dry_combos = pickup(dry_combo_weight, heli.weight_capacity);
            double rem_w = heli.weight_capacity - dry_combos * dry_combo_weight;
            int wet_combos = 0;
            if (rem_w > EPS) wet_combos = (int)std::floor(rem_w / wet_combo_weight);

            // Available packets on this trip (packets, not combos)
            int dry_food_avail = dry_combos * 9;          // number of dry packets
            int perishable_food_avail = wet_combos * 9;   // number of perishable packets
            int other_avail = dry_combos + wet_combos;   // number of other supplies

            // Record what was picked up on this trip (used later by objective checker to compare)
            Trip trip;
            trip.dry_food_pickup = dry_food_avail;
            trip.perishable_food_pickup = perishable_food_avail;
            trip.other_supplies_pickup = other_avail;

            // trip-local state
            Point currentPos = home;
            double trip_distance = 0.0;
            bool trip_has_drops = false;

            // Serve villages as long as we have load and distance allowance
            while (idx < assigned.size()) {
                int vid2 = assigned[idx];
                if (vid2 < 1 || vid2 > numVillages) { ++idx; continue; }
                if (people_remaining[vid2] <= 0) { ++idx; continue; }

                Point villagePos = problem.villages[vid2 - 1].coords;

                // Check if this trip can physically visit this village and return to base within per-trip capacity
                double dist_to_village = distance(currentPos, villagePos);
                double dist_village_to_home = distance(villagePos, home);
                double projected_trip_if_visit = trip_distance + dist_to_village + dist_village_to_home;

                if (projected_trip_if_visit > heli.distance_capacity + EPS) {
                    // Can't include this village in CURRENT trip.
                    // If trip is empty (no drops added yet) then this helicopter cannot reach this village even from home.
                    if (!trip_has_drops) {
                        // If even a fresh trip from home can't reach it, stop trying with this helicopter.
                        double roundtrip_from_home = distance(home, villagePos) + distance(villagePos, home);
                        if (roundtrip_from_home > heli.distance_capacity + EPS) {
                            // unreachable by this helicopter; skip this village (it was assigned in state, but not feasible)
                            ++idx; // skip to next assigned village
                            continue;
                        } else {
                            // Strange case: we computed projected_trip_if_visit > dcap but trip had no drops;
                            // treat as cannot serve in this trip and break to attempt new trip (though new trip would be same)
                            break;
                        }
                    } else {
                        // finish this trip by returning to base and commit it
                        trip_distance += distance(currentPos, home);
                        heli_total_distance += trip_distance;
                        if (heli_total_distance > problem.d_max + EPS) {
                            // helicopter exceeded total distance budget -> stop all work
                            goto finish_helicopter;
                        }
                        if (!trip.drops.empty()) plan.trips.push_back(trip);

                        // start a fresh trip (outer while will do so)
                        break;
                    }
                }

                // How many people can current load serve?
                int total_food_avail = dry_food_avail + perishable_food_avail; // packets
                int max_people_by_food = total_food_avail / 9;                // floor
                int max_people_by_other = other_avail;
                int max_people_served_now = std::min(max_people_by_food, max_people_by_other);

                if (max_people_served_now <= 0) {
                    // No capacity left on this trip. Finish trip (return to base) and start new one.
                    if (!trip_has_drops) {
                        // Trip had no capacity at all -> nothing to do; start fresh but if still 0 combos, we cannot proceed
                        // try to start a new trip once (the pickup() might produce 0 again); we'll break to restart above
                        break;
                    } else {
                        // end current trip
                        trip_distance += distance(currentPos, home);
                        heli_total_distance += trip_distance;
                        if (heli_total_distance > problem.d_max + EPS) goto finish_helicopter;
                        if (!trip.drops.empty()) plan.trips.push_back(trip);
                        break; // outer loop will start new trip
                    }
                }

                // Determine how many people to serve at this village in this trip
                int people_here = people_remaining[vid2];
                int serve_count = std::min(people_here, max_people_served_now);

                // Allocate food packets: prefer perishable first (higher value)
                int need_food_packets = serve_count * 9;
                int use_perishable = std::min(perishable_food_avail, need_food_packets);
                need_food_packets -= use_perishable;
                int use_dry = std::min(dry_food_avail, need_food_packets);
                need_food_packets -= use_dry;

                // Sanity check (should be zero now)
                if (need_food_packets > 0) {
                    // Unexpected: not enough food (shouldn't happen because serve_count limited by food)
                    serve_count -= (int)std::ceil((double)need_food_packets / 9.0);
                    if (serve_count <= 0) {
                        // can't serve anyone
                        break;
                    }
                    // recompute allocations for reduced serve_count
                    need_food_packets = serve_count * 9;
                    use_perishable = std::min(perishable_food_avail, need_food_packets);
                    need_food_packets -= use_perishable;
                    use_dry = std::min(dry_food_avail, need_food_packets);
                    need_food_packets -= use_dry;
                }

                // Other supplies allocation (one per person)
                int use_other = std::min(other_avail, serve_count);
                if (use_other < serve_count) {
                    // can't fully serve serve_count people because of other limits -> reduce
                    int new_serve = use_other;
                    if (new_serve <= 0) {
                        // no other supplies left -> finish trip
                        if (!trip_has_drops) {
                            // nothing possible this trip
                            break;
                        } else {
                            trip_distance += distance(currentPos, home);
                            heli_total_distance += trip_distance;
                            if (heli_total_distance > problem.d_max + EPS) goto finish_helicopter;
                            if (!trip.drops.empty()) plan.trips.push_back(trip);
                            break;
                        }
                    }
                    serve_count = new_serve;
                    // recompute food allocation for the reduced serve_count
                    need_food_packets = serve_count * 9;
                    use_perishable = std::min(perishable_food_avail, need_food_packets);
                    need_food_packets -= use_perishable;
                    use_dry = std::min(dry_food_avail, need_food_packets);
                    need_food_packets -= use_dry;
                }

                // Commit allocation as a Drop
                Drop drop;
                drop.village_id = vid2;
                drop.perishable_food = use_perishable;
                drop.dry_food = use_dry;
                drop.other_supplies = serve_count; // exactly one per person covered

                trip.drops.push_back(drop);
                trip_has_drops = true;

                // Update remaining supplies on the trip
                perishable_food_avail -= use_perishable;
                dry_food_avail -= use_dry;
                other_avail -= serve_count;

                // Update people remaining at village
                people_remaining[vid2] -= serve_count; // may become zero or positive

                // Move to village: update trip distance and current position
                trip_distance += dist_to_village;
                currentPos = villagePos;

                // If village still has people, we need to return to base (trip ended) to reload
                if (people_remaining[vid2] > 0) {
                    // finish trip by returning home
                    trip_distance += dist_village_to_home;
                    heli_total_distance += trip_distance;
                    if (heli_total_distance > problem.d_max + EPS) goto finish_helicopter;
                    if (!trip.drops.empty()) plan.trips.push_back(trip);
                    // start a new trip (outer loop)
                    break;
                } else {
                    // village fully served, move to next assigned village (stay on same trip if capacity remains)
                    ++idx;
                    // continue serving next village in same trip (loop)
                }
            } // end while serving this trip

            // if we finished inner trip loop without pushing trip (because we exited due to capacity or completed villages),
            // and if trip has drops but wasn't yet committed, we must commit it and add return distance.
            if (trip_has_drops && (trip.drops.size() > 0)) {
                // If the trip hasn't been committed above (we might have committed inside), ensure we return home and count distance
                // But avoid double-adding distance if already added above (we added to heli_total_distance when we pushed).
                // Here we check if last commit has already accounted for this trip by comparing distances; simplest approach:
                // If currentPos != home, return and add distance.
                if (!(currentPos.x == home.x && currentPos.y == home.y)) {
                    trip_distance += distance(currentPos, home);
                    heli_total_distance += trip_distance;
                    if (heli_total_distance > problem.d_max + EPS) {
                        // trip exceeded overall budget -> discard this trip to keep feasibility
                        // (alternatively we could remove trip, but here we'll not push it)
                    } else {
                        plan.trips.push_back(trip);
                    }
                } else {
                    // already at home (rare), but still push if not already pushed
                    plan.trips.push_back(trip);
                }
            }

            // if heli_total_distance is over global per-heli budget, stop
            if (heli_total_distance > problem.d_max + EPS) {
                break;
            }
        } // end outer while over assigned villages

    finish_helicopter:
        solution[h] = plan;
    } // end for each helicopter

    return solution;
}

double calculateObjectiveValue(const Solution& solution) const {
    double total_value = 0.0;
    double total_cost = 0.0;

    int num_helicopters = (int)problem.helicopters.size();
    int num_villages = (int)problem.villages.size();

    vector<double> food_delivered(num_villages + 1, 0.0);
    vector<double> other_delivered(num_villages + 1, 0.0);
    vector<double> helicopter_total_distance(num_helicopters, 0.0);

    for (int h = 0; h < num_helicopters; ++h) {
        const auto& plan = solution[h];
        const auto& helicopter = problem.helicopters[plan.helicopter_id-1];

        if (plan.trips.empty()) continue; // unused heli → no cost

        // fixed cost charged once per helicopter
        total_cost += helicopter.fixed_cost;

        for (const auto& trip : plan.trips) {
            if (trip.drops.empty()) continue;

            // --- Weight check ---
            double trip_weight =
                trip.dry_food_pickup * problem.packages[0].weight +
                trip.perishable_food_pickup * problem.packages[1].weight +
                trip.other_supplies_pickup * problem.packages[2].weight;

            if (trip_weight > helicopter.weight_capacity + 1e-6) return -1e9;

            // --- Package conservation check ---
            int total_dry_drops = 0, total_perishable_drops = 0, total_other_drops = 0;
            for (const auto& drop : trip.drops) {
                total_dry_drops += drop.dry_food;
                total_perishable_drops += drop.perishable_food;
                total_other_drops += drop.other_supplies;
            }

            if (total_dry_drops > trip.dry_food_pickup ||
                total_perishable_drops > trip.perishable_food_pickup ||
                total_other_drops > trip.other_supplies_pickup) {
                return -1e9;
            }

            // --- Distance check ---
            vector<int> route;
            route.reserve(trip.drops.size());
            for (const auto& drop : trip.drops) route.push_back(drop.village_id);

            double trip_distance = calculateTourDistance(plan.helicopter_id - 1, route);
            if (trip_distance > helicopter.distance_capacity + 1e-6) return -1e9;

            helicopter_total_distance[h] += trip_distance;
            if (helicopter_total_distance[h] > problem.d_max + 1e-6) return -1e9;

            // --- Value accumulation ---
            for (const auto& drop : trip.drops) {
                total_value += calculateVillageValue(
                    drop.village_id, drop.dry_food, drop.perishable_food, drop.other_supplies,
                    food_delivered[drop.village_id], other_delivered[drop.village_id]
                );
                food_delivered[drop.village_id] += drop.dry_food + drop.perishable_food;
                other_delivered[drop.village_id] += drop.other_supplies;
            }

            // --- Variable cost (distance) ---
            total_cost += helicopter.alpha * trip_distance;
        }
    }

    return total_value - total_cost;
}



    Solution solveProblem() {
        Solution best_solution;
        double best_score = -1e18;

        auto start_time = chrono::steady_clock::now();
        auto time_limit_ms = (long long)(problem.time_limit_minutes * 60.0 * 1000.0 * 0.9);
        if (time_limit_ms < 1000) time_limit_ms = 1000;

        int restart_count = 0;
        
        while (chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - start_time).count() < time_limit_ms) {
            restart_count++;
            
            State current_state = Random_State(problem);
            Solution current_solution = stateToSolution(current_state, problem);
            double current_score = calculateObjectiveValue(current_solution);

            int iterations_without_improvement = 0;
            const int max_iterations = 1000;

            while (iterations_without_improvement < max_iterations &&
                   chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - start_time).count() < time_limit_ms) {
                
                State neighbor_state = perturbState(current_state);
                Solution neighbor_solution = stateToSolution(neighbor_state, problem);
                double neighbor_score = calculateObjectiveValue(neighbor_solution);

                if (neighbor_score > current_score) {
                    current_state = move(neighbor_state);
                    current_solution = move(neighbor_solution);
                    current_score = neighbor_score;
                    iterations_without_improvement = 0;
                } else {
                    iterations_without_improvement++;
                }
            }

            if (current_score > best_score) {
                best_score = current_score;
                best_solution = current_solution;
                cerr << "Restart " << restart_count << ", new best score = " << best_score << endl;
            }
        }

        cerr << "Total restarts: " << restart_count << ", final best score = " << best_score << endl;
        return best_solution;
    }
};

Solution solve(const ProblemData& problem) {
    LocalSearchSolver solver(problem);
    Solution solution = solver.solveProblem();
    return solution;
}

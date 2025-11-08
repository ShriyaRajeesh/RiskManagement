#include <iostream>
#include <vector>
#include <random>
#include <omp.h>
#include <fstream>
#include <iomanip>
using namespace std;

// Function to simulate portfolio value given a base value, a shock, and iterations (days)
double simulate_portfolio(double base_value, double shock_percent, int iterations)
{
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<> noise(0, 0.02); // small random daily fluctuation

    double value = base_value;
    for (int i = 0; i < iterations; i++)
    {
        double daily_return = shock_percent + noise(gen);
        value *= (1 + daily_return);
    }
    return value;
}

int main()
{
    // ---------------- Configuration ----------------
    const int NUM_SCENARIOS = 8;
    const int DAYS = 252;              // trading days in a year
    const double BASE_PORTFOLIO = 1e6; // starting portfolio value

    // Different stress test shocks (as daily return offsets)
    vector<double> shocks = {-0.10, -0.20, -0.30, 0.05, -0.05, 0.15, -0.25, 0.10};
    vector<double> results(NUM_SCENARIOS);

    cout << "Running parallel stress test simulation...\n";
    double start_time = omp_get_wtime();

// Run scenarios in parallel
#pragma omp parallel for
    for (int i = 0; i < NUM_SCENARIOS; i++)
    {
        results[i] = simulate_portfolio(BASE_PORTFOLIO, shocks[i], DAYS);
    }

    double end_time = omp_get_wtime();

    // ---------------- Output to CSV ----------------
    ofstream file("results.csv");
    if (!file.is_open())
    {
        cerr << "Error: Unable to open results.csv for writing.\n";
        return 1;
    }

    file << "Scenario,Shock(%),FinalValue\n";
    for (int i = 0; i < NUM_SCENARIOS; i++)
    {
        file << (i + 1) << ","
             << fixed << setprecision(2) << shocks[i] * 100 << ","
             << setprecision(4) << results[i] << "\n";
    }
    file.close();

    cout << "Stress test complete!\n";
    cout << "Results saved to results.csv\n";
    cout << "Total time: " << fixed << setprecision(3)
         << (end_time - start_time) << " seconds\n";

    return 0;
}

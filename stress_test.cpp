// // stress_test.cpp
// // Full Monte Carlo stress-test engine with OpenMP, correlations (Cholesky), per-thread buffers,
// // CSV output for sample returns and scenario summary.
// // Compile: g++ -O3 -std=c++17 -fopenmp stress_test.cpp -o stress_test
// // Run:    OMP_NUM_THREADS=8 ./stress_test

// #include <bits/stdc++.h>
// #include <omp.h>
// using namespace std;

// // ----------------- Utility RNG wrapper (per-thread) -----------------
// struct RNGWrapper {
//     std::mt19937_64 gen;
//     std::normal_distribution<double> nd;
//     RNGWrapper(unsigned long long seed=1ULL) : gen((uint64_t)seed), nd(0.0,1.0) {}
//     double normal() { return nd(gen); }
// };

// // ----------------- Simple Cholesky for positive-definite matrices -----------------
// // Compute lower-triangular L such that A = L * L^T
// // A is n x n, stored row-major as vector<vector<double>>
// bool cholesky_decompose(const vector<vector<double>>& A, vector<vector<double>>& L) {
//     int n = (int)A.size();
//     L.assign(n, vector<double>(n, 0.0));
//     for (int i = 0; i < n; ++i) {
//         for (int j = 0; j <= i; ++j) {
//             double s = 0.0;
//             for (int k = 0; k < j; ++k) s += L[i][k] * L[j][k];
//             if (i == j) {
//                 double val = A[i][i] - s;
//                 if (val <= 0.0) return false;
//                 L[i][j] = sqrt(val);
//             } else {
//                 L[i][j] = (1.0 / L[j][j]) * (A[i][j] - s);
//             }
//         }
//     }
//     return true;
// }

// // ----------------- Scenario struct -----------------
// struct Scenario {
//     string name;
//     vector<double> mu_shift; // additive shift to mu per asset
//     double scale;            // multiplicative scale on portfolio return
//     Scenario(const string &n, const vector<double> &shift, double s=1.0) : name(n), mu_shift(shift), scale(s) {}
// };

// // ----------------- CSV helpers -----------------
// void write_csv_sample(const string &fname, const vector<double> &vals) {
//     ofstream ofs(fname);
//     if (!ofs.is_open()) {
//         cerr << "Cannot open file for writing: " << fname << "\n";
//         return;
//     }
//     ofs << "return\n";
//     for (double v : vals) ofs << std::setprecision(12) << v << "\n";
//     ofs.close();
// }

// // Append summary line to summary CSV
// void append_summary_csv(const string &fname, const string &scenario_name, int n_sim, double mean, double sd, double worst, double best, double VaR, double CVaR) {
//     bool exists = ifstream(fname).good();
//     ofstream ofs(fname, ios::app);
//     if (!ofs.is_open()) {
//         cerr << "Cannot open summary file: " << fname << "\n";
//         return;
//     }
//     if (!exists) {
//         ofs << "scenario,n_sim,mean,sd,worst,best,VaR95,CVaR95\n";
//     }
//     ofs << "\"" << scenario_name << "\"," << n_sim << "," << setprecision(12) << mean << "," << sd << "," << worst << "," << best << "," << VaR << "," << CVaR << "\n";
//     ofs.close();
// }

// // ----------------- Main Monte Carlo -----------------
// int main() {
//     ios::sync_with_stdio(false);
//     cin.tie(nullptr);

//     // ---------------- CONFIG (tweak these) ----------------
//     int n_assets = 5;
//     vector<double> weights = {0.3, 0.25, 0.2, 0.15, 0.1}; // sum to 1
//     vector<double> mu      = {0.06, 0.05, 0.04, 0.03, 0.02}; // expected returns (annual)
//     vector<double> sigma   = {0.15, 0.18, 0.12, 0.20, 0.10}; // volatilities (annual)

//     // Optional correlation matrix (n_assets x n_assets) - 1 on diagonal
//     // Example: low-to-moderate correlations
//     vector<vector<double>> rho = {
//         {1.0, 0.3, 0.2, 0.25, 0.15},
//         {0.3, 1.0, 0.25, 0.3, 0.2},
//         {0.2, 0.25, 1.0, 0.15, 0.1},
//         {0.25, 0.3, 0.15, 1.0, 0.2},
//         {0.15, 0.2, 0.1, 0.2, 1.0}
//     };

//     // Number of Monte Carlo samples (per scenario)
//     int n_sim = 500000; // change as desired (be aware of memory)
//     // If you want less memory usage, reduce n_sim or store only a subset.
//     // Number of sample returns to keep for plotting per scenario (cap)
//     int keep_sample = 100000;

//     // Create scenarios
//     vector<Scenario> scenarios;
//     scenarios.emplace_back("Baseline", vector<double>(n_assets, 0.0), 1.0);
//     scenarios.emplace_back("Rate-Hike+200bp", vector<double>{-0.02, -0.018, -0.015, -0.01, -0.005}, 1.0);
//     scenarios.emplace_back("Market-Crash-30pct", vector<double>(n_assets, 0.0), 0.7);
//     scenarios.emplace_back("Rate+Crash", vector<double>{-0.02, -0.018, -0.015, -0.01, -0.005}, 0.7);

//     // Output files
//     string summary_fname = "scenario_summary.csv";

//     // Precompute covariance matrix: Cov[i][j] = sigma[i] * sigma[j] * rho[i][j]
//     vector<vector<double>> cov(n_assets, vector<double>(n_assets, 0.0));
//     for (int i = 0; i < n_assets; ++i) {
//         for (int j = 0; j < n_assets; ++j) {
//             cov[i][j] = sigma[i] * sigma[j] * rho[i][j];
//         }
//     }

//     // Compute Cholesky L (lower triangular) of covariance
//     vector<vector<double>> L;
//     bool ok = cholesky_decompose(cov, L);
//     if (!ok) {
//         cerr << "Cholesky decomposition failed (covariance not positive definite). Exiting.\n";
//         return 1;
//     }

//     cout << "Starting Monte Carlo stress test\n";
//     cout << "Assets: " << n_assets << "  sims per scenario: " << n_sim << "  threads(omp): " << omp_get_max_threads() << "\n";

//     // Iterate scenarios
//     for (auto &sc : scenarios) {
//         cout << "Running scenario: " << sc.name << " ...\n";
//         // Containers for results: per-thread local storage then merge
//         int n_threads = omp_get_max_threads();
//         vector<vector<double>> thread_returns(n_threads);

//         // Reserve approximate capacity per thread
//         int per_thread_est = (n_sim + n_threads - 1) / n_threads;
//         for (int t = 0; t < n_threads; ++t) thread_returns[t].reserve(per_thread_est);

//         // Timer
//         double t0 = omp_get_wtime();

//         // Base seed (use steady clock)
//         unsigned long long base_seed = (unsigned long long)std::chrono::high_resolution_clock::now().time_since_epoch().count();

//         // Parallel region
//         #pragma omp parallel
//         {
//             int tid = omp_get_thread_num();
//             unsigned long long seed = base_seed + (unsigned long long)(tid * 11400714819323198485ULL); // large stride
//             RNGWrapper rng(seed);

//             // Each thread does its portion
//             #pragma omp for schedule(static)
//             for (int i = 0; i < n_sim; ++i) {
//                 // generate independent normals z (size n_assets)
//                 vector<double> z(n_assets);
//                 for (int a = 0; a < n_assets; ++a) z[a] = rng.normal();

//                 // correlate: z_corr = L * z (L is lower-triangular)
//                 vector<double> zcorr(n_assets, 0.0);
//                 for (int r = 0; r < n_assets; ++r) {
//                     double s = 0.0;
//                     for (int c = 0; c <= r; ++c) s += L[r][c] * z[c];
//                     zcorr[r] = s;
//                 }

//                 // compute portfolio return : r_a = mu[a] + mu_shift[a] + zcorr[a]
//                 double port_ret = 0.0;
//                 for (int a = 0; a < n_assets; ++a) {
//                     double mu_eff = mu[a] + sc.mu_shift[a];
//                     double ra = mu_eff + zcorr[a];
//                     port_ret += weights[a] * ra;
//                 }
//                 port_ret *= sc.scale;

//                 thread_returns[tid].push_back(port_ret);
//             } // end for
//         } // end parallel

//         double t1 = omp_get_wtime();
//         cout << "  simulation done in " << (t1 - t0) << " s\n";

//         // Merge all thread_returns into a single vector all_returns
//         vector<double> all_returns;
//         size_t total = 0;
//         for (int t = 0; t < n_threads; ++t) total += thread_returns[t].size();
//         all_returns.reserve(total);
//         for (int t = 0; t < n_threads; ++t) {
//             all_returns.insert(all_returns.end(), thread_returns[t].begin(), thread_returns[t].end());
//             // optionally free thread vector memory
//             vector<double>().swap(thread_returns[t]);
//         }

//         if ((int)all_returns.size() != n_sim) {
//             cerr << "Warning: merged returns size (" << all_returns.size() << ") != n_sim (" << n_sim << ")\n";
//         }

//         // Sort to compute VaR/CVaR
//         sort(all_returns.begin(), all_returns.end()); // ascending: worst (lowest returns) first

//         // Basic stats
//         double mean = 0.0;
//         for (double v : all_returns) mean += v;
//         mean /= all_returns.size();

//         double sd = 0.0;
//         for (double v : all_returns) sd += (v - mean) * (v - mean);
//         sd = sqrt(sd / max(1.0, (double)all_returns.size() - 1.0));

//         double worst = all_returns.front();
//         double best  = all_returns.back();

//         // VaR (95% -> 5th percentile): index floor(p * n) - 1 (0-based)
//         double p = 0.05;
//         int idx = max(0, (int)floor(p * all_returns.size()) - 1);
//         if (idx < 0) idx = 0;
//         double VaR = all_returns[idx];

//         // CVaR (Expected Shortfall): average of worst p fraction
//         int k = max(1, (int)ceil(p * all_returns.size()));
//         double sum_tail = 0.0;
//         for (int i = 0; i < k; ++i) sum_tail += all_returns[i];
//         double CVaR = sum_tail / k;

//         cout << "  mean=" << mean << "  sd=" << sd << "  worst=" << worst << "  best=" << best << "\n";
//         cout << "  VaR(95%)=" << VaR << "  CVaR(95%)=" << CVaR << "\n";

//         // Store sample subset to CSV for plotting
//         vector<double> sample_for_plot;
//         sample_for_plot.reserve(min((size_t)keep_sample, all_returns.size()));
//         // pick equispaced samples from sorted returns for variety, or just copy first keep_sample random entries
//         // We'll pick a random subset uniformly across all_returns to show distribution (not just tails)
//         std::mt19937_64 pickgen((uint64_t)base_seed ^ 0x5bf03635);
//         std::uniform_int_distribution<size_t> idxdist(0, all_returns.size()-1);
//         size_t to_keep = min((size_t)keep_sample, all_returns.size());
//         for (size_t i = 0; i < to_keep; ++i) {
//             size_t j = idxdist(pickgen);
//             sample_for_plot.push_back(all_returns[j]);
//         }

//         string safe_name = sc.name;
//         // sanitize file name
//         for (char &c : safe_name) if (!(isalnum((unsigned char)c) || c=='-' || c=='_' || c==' ')) c = '_';
//         // Replace spaces with underscores for filenames
//         for (char &c : safe_name) if (c == ' ') c = '_';

//         string sample_fname = "returns_" + safe_name + ".csv";
//         write_csv_sample(sample_fname, sample_for_plot);
//         cout << "  wrote sample returns to " << sample_fname << " (" << sample_for_plot.size() << " rows)\n";

//         // Append to summary CSV
//         append_summary_csv(summary_fname, sc.name, (int)all_returns.size(), mean, sd, worst, best, VaR, CVaR);
//         cout << "  appended summary to " << summary_fname << "\n\n";

//         // Optional: free memory
//         vector<double>().swap(all_returns);
//         vector<double>().swap(sample_for_plot);
//     } // end scenario loop

//     cout << "All scenarios complete.\n";
//     return 0;
// }
#include <iostream>
#include <vector>
#include <random>
#include <omp.h>
#include <fstream>
using namespace std;

double simulate_portfolio(double base_value, double shock_percent, int iterations) {
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<> noise(0, 0.02); // small randomness for realism

    double value = base_value;
    for (int i = 0; i < iterations; i++) {
        double daily_return = shock_percent + noise(gen);
        value *= (1 + daily_return);
    }
    return value;
}

int main() {
    const int NUM_SCENARIOS = 8;
    const int DAYS = 252;
    double base_portfolio = 1e6;

    vector<double> shocks = {-0.10, -0.20, -0.30, 0.05, -0.05, 0.15, -0.25, 0.10};
    vector<double> results(NUM_SCENARIOS);

    double start_time = omp_get_wtime();

    #pragma omp parallel for
    for (int i = 0; i < NUM_SCENARIOS; i++) {
        results[i] = simulate_portfolio(base_portfolio, shocks[i], DAYS);
    }

    double end_time = omp_get_wtime();

    // Save results to CSV file
    ofstream file("results.csv");
    file << "Scenario,Shock(%),FinalValue\n";
    for (int i = 0; i < NUM_SCENARIOS; i++) {
        file << (i + 1) << "," << shocks[i] * 100 << "," << results[i] << "\n";
    }
    file.close();

    cout << "Stress test complete! Results saved to results.csv\n";
    cout << "Time taken: " << end_time - start_time << " seconds\n";
    return 0;
}

// loa_beale.cpp
// LOA (exact equations from the paper) applied to the Beale function
// Compile: g++ -std=c++17 -O2 loa_beale.cpp -o loa_beale
// Run: ./loa_beale

#include <bits/stdc++.h>
using namespace std;

// ---------------- Random Generator ----------------
mt19937 rng((unsigned)time(NULL));
double rnd(double a, double b) {
    uniform_real_distribution<double> dist(a, b);
    return dist(rng);
}

// ---------------- Beale Function ----------------
// f(x,y) = (1.5 - x + x*y)^2 + (2.25 - x + x*y^2)^2 + (2.625 - x + x*y^3)^2
double beale(const vector<double>& v) {
    double x = v[0];
    double y = v[1];
    double t1 = 1.5 - x + x*y;
    double t2 = 2.25 - x + x*y*y;
    double t3 = 2.625 - x + x*y*y*y;
    return t1*t1 + t2*t2 + t3*t3;
}

// ---------------- LOA PARAMETERS ----------------
int POP = 30;       // population size
int DIM = 2;        // Beale is 2D
int MAX_IT = 6000;  // iterations
double LB = -4.5;   // lower bound
double UB = 4.5;    // upper bound

// ---------------- Initialization (Eq 2) ----------------
vector<vector<double>> init_pop() {
    vector<vector<double>> P(POP, vector<double>(DIM));
    for (int i = 0; i < POP; ++i)
        for (int j = 0; j < DIM; ++j)
            P[i][j] = rnd(LB, UB);
    return P;
}

// ---------------- Escape Update (Eq 6) ----------------
// xP1_i,j = x_i,j + r_{i,j} * (SSA_{i,j} - I_{i,j} * x_i,j)
vector<double> escape_exact(const vector<double>& Xi, const vector<double>& SSA) {
    vector<double> newX = Xi;
    for (int j = 0; j < DIM; ++j) {
        double r = rnd(0.0, 1.0);
        int I = (rng() % 2) + 1; // 1 or 2
        newX[j] = Xi[j] + r * (SSA[j] - I * Xi[j]);
        // bound correction
        if (newX[j] < LB) newX[j] = LB;
        if (newX[j] > UB) newX[j] = UB;
    }
    return newX;
}

// ---------------- Hide Update (Eq 8) ----------------
// xP2_i,j = x_i,j + (1 - 2*r_{i,j}) * (ub_j - lb_j) / t
vector<double> hide_exact(const vector<double>& Xi, int t) {
    vector<double> newX = Xi;
    for (int j = 0; j < DIM; ++j) {
        double r = rnd(0.0, 1.0);
        newX[j] = Xi[j] + (1.0 - 2.0 * r) * ((UB - LB) / (double)t);
        // bound correction
        if (newX[j] < LB) newX[j] = LB;
        if (newX[j] > UB) newX[j] = UB;
    }
    return newX;
}

// ---------------- MAIN LOA ----------------
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    auto pop = init_pop();
    vector<double> fit(POP);

    // Evaluate initial population
    for (int i = 0; i < POP; ++i)
        fit[i] = beale(pop[i]);

    // Track best solution
    double bestFit = numeric_limits<double>::infinity();
    vector<double> bestSol(DIM);

    for (int i = 0; i < POP; ++i) {
        if (fit[i] < bestFit) {
            bestFit = fit[i];
            bestSol = pop[i];
        }
    }

    // Main LOA loop (Eq 4,5,6,7,8,9)
    for (int t = 1; t <= MAX_IT; ++t) {
        for (int i = 0; i < POP; ++i) {

            // Safe areas SA_i (Eq 5): indices k with F_k < F_i
            vector<int> safe;
            safe.reserve(POP);
            for (int k = 0; k < POP; ++k)
                if (fit[k] < fit[i])
                    safe.push_back(k);

            vector<double> candidate;
            double rp = rnd(0.0, 1.0);

            // Update decision (Eq 4)
            if (rp <= 0.5 && !safe.empty()) {
                // Phase 1: Escape (Eq 6)
                int idx = safe[rng() % safe.size()];
                candidate = escape_exact(pop[i], pop[idx]);
            } else {
                // Phase 2: Hide (Eq 8)
                // ensure t > 0 (t starts at 1)
                candidate = hide_exact(pop[i], max(1, t));
            }

            double newFit = beale(candidate);

            // Accept if improved (Eq 7 & Eq 9)
            if (newFit <= fit[i]) {
                pop[i] = candidate;
                fit[i] = newFit;
                if (newFit < bestFit) {
                    bestFit = newFit;
                    bestSol = candidate;
                }
            }
        }

        // optional progress print
        if (t % 500 == 0 || t == 1) {
            cout << "Iter " << t << " | Best = " << scientific << setprecision(8) << bestFit << '\n';
        }
    }

    // Final results
    cout << fixed << setprecision(10);
    cout << "\nBest Solution Found:\n";
    cout << "x = " << bestSol[0] << "\n";
    cout << "y = " << bestSol[1] << "\n";
    cout << "\nBeale Value = " << bestFit << "\n";

    return 0;
}

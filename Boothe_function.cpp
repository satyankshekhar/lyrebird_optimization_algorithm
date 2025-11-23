#include <bits/stdc++.h>
using namespace std;

// ---------------- Random Generator ----------------
mt19937 rng(time(NULL));
double rnd(double a, double b) {
    uniform_real_distribution<double> dist(a, b);
    return dist(rng);
}

// ---------------- Booth Function ----------------
// f(x,y) = (x + 2y - 7)^2 + (2x + y - 5)^2
double booth(const vector<double>& x) {
    double X = x[0];
    double Y = x[1];

    double t1 = X + 2.0 * Y - 7.0;
    double t2 = 2.0 * X + Y - 5.0;

    return t1*t1 + t2*t2;
}

// ---------------- LOA PARAMETERS ----------------
int POP = 30;       // population size
int DIM = 2;        // Booth is 2D
int MAX_IT = 6000;  // iterations
double LB = -10;    // lower bound
double UB = 10;     // upper bound

// ---------------- Initialization (Eq 2) ----------------
vector<vector<double>> init_pop() {
    vector<vector<double>> P(POP, vector<double>(DIM));
    for (int i = 0; i < POP; i++)
        for (int j = 0; j < DIM; j++)
            P[i][j] = rnd(LB, UB);
    return P;
}

// ---------------- Escape Update (Eq 6) ----------------
// xP1_i,j = x_i,j + r (SSA_j - I * x_i,j)
vector<double> escape_exact(const vector<double> &Xi, const vector<double> &SSA) {
    vector<double> newX = Xi;

    for (int j = 0; j < DIM; j++) {
        double r = rnd(0, 1);
        int I = (rng() % 2) + 1;  // I = 1 or 2

        newX[j] = Xi[j] + r * (SSA[j] - I * Xi[j]);

        // Bound correction
        newX[j] = max(LB, min(UB, newX[j]));
    }
    return newX;
}

// ---------------- Hide Update (Eq 8) ----------------
// xP2_i,j = x_i,j + (1 - 2*r) * (UB - LB) / t
vector<double> hide_exact(const vector<double> &Xi, int t) {
    vector<double> newX = Xi;

    for (int j = 0; j < DIM; j++) {
        double r = rnd(0, 1);

        newX[j] = Xi[j] + (1 - 2*r) * ((UB - LB) / (double)t);

        // Bound correction
        newX[j] = max(LB, min(UB, newX[j]));
    }
    return newX;
}

// ---------------- MAIN LOA ----------------
int main() {

    auto pop = init_pop();
    vector<double> fit(POP);

    // Evaluate initial population
    for (int i = 0; i < POP; i++)
        fit[i] = booth(pop[i]);

    // Track best solution
    double bestFit = 1e18;
    vector<double> bestSol(DIM);

    for (int i = 0; i < POP; i++)
        if (fit[i] < bestFit) {
            bestFit = fit[i];
            bestSol = pop[i];
        }

    // ---------------- LOA Loop ----------------
    for (int t = 1; t <= MAX_IT; t++) {

        for (int i = 0; i < POP; i++) {

            // Safe areas (Eq 5)
            vector<int> safe;
            for (int k = 0; k < POP; k++)
                if (fit[k] < fit[i])
                    safe.push_back(k);

            vector<double> candidate;
            double rp = rnd(0, 1);

            // Decide escape or hide (Eq 4)
            if (rp <= 0.5 && !safe.empty()) {
                // Escape (Eq 6)
                int idx = safe[rng() % safe.size()];
                candidate = escape_exact(pop[i], pop[idx]);
            }
            else {
                // Hide (Eq 8)
                candidate = hide_exact(pop[i], t);
            }

            double newFit = booth(candidate);

            // Update if improved (Eq 7 and Eq 9)
            if (newFit <= fit[i]) {
                pop[i] = candidate;
                fit[i] = newFit;

                if (newFit < bestFit) {
                    bestFit = newFit;
                    bestSol = candidate;
                }
            }
        }

        if (t % 500 == 0)
            cout << "Iter " << t << " | Best = " << bestFit << endl;
    }

    // Final Best Output
    cout << "\nBest Solution Found:\n";
    cout << "x = " << bestSol[0] << endl;
    cout << "y = " << bestSol[1] << endl;

    cout << "\nBooth Value = " << bestFit << endl;

    return 0;
}

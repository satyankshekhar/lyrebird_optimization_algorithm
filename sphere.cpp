#include <bits/stdc++.h>
using namespace std;

// Random generator
mt19937 rng(time(NULL));
double randF(double a, double b)
{
    uniform_real_distribution<double> dist(a, b);
    return dist(rng);
}

// ------------------------------
// Sphere Function
// ------------------------------
double sphere(const vector<double> &x)
{
    double s = 0;
    for (double v : x)
        s += v * v;
    return s;
}

// ------------------------------
// LOA PARAMETERS
// ------------------------------
int POP = 20;        // population size
int DIM = 5;         // dimension of sphere function
int MAX_IT = 200; // number of iterations
double LB = -100;    // lower bound
double UB = 100;     // upper bound

// ------------------------------
// Initialize population
// ------------------------------
vector<vector<double>> init_population()
{
    vector<vector<double>> pop(POP, vector<double>(DIM));
    for (int i = 0; i < POP; i++)
        for (int d = 0; d < DIM; d++)
            pop[i][d] = randF(LB, UB);
    return pop;
}

// ------------------------------
// Escape move (global search)
// ------------------------------
vector<double> escape(const vector<double> &x,
                      const vector<double> &SSA)
{
    int DIM = x.size();
    vector<double> newX(DIM);

    for (int j = 0; j < DIM; j++)
    {
        double r = randF(0, 1);
        int I = (rng() % 2) + 1; // I = 1 or 2

        newX[j] = x[j] + r * (SSA[j] - I * x[j]);

        // Bound check
        if (newX[j] < LB)
            newX[j] = LB;
        if (newX[j] > UB)
            newX[j] = UB;
    }

    return newX;
}

// ------------------------------
// Hide move (local search)
// ------------------------------
vector<double> hide(const vector<double> &Xi, int t) {
    vector<double> newX = Xi;

    for (int j = 0; j < DIM; j++) {
        double r = randF(0, 1);

        newX[j] = Xi[j] + (1 - 2*r) * (UB - LB) / t;

        // bounds
        if (newX[j] < LB) newX[j] = LB;
        if (newX[j] > UB) newX[j] = UB;
    }
    return newX;
}

// ------------------------------
// MAIN LOA ALGORITHM
// ------------------------------
int main()
{
    auto pop = init_population(); // random solution create kia idhar
    vector<double> fitness(POP);  // sabhi lyrebird ki fitness store krne k liye

    // Evaluate initial fitness
    for (int i = 0; i < POP; i++)
        fitness[i] = sphere(pop[i]); // initial fitness calculate kr liya

    // Track best
    double bestFit = 1e18; // initial best fitness ko bahut bada rkh dia

    vector<double> bestSol(DIM);

    for (int i = 0; i < POP; i++)
    {
        if (fitness[i] < bestFit)
        {
            bestFit = fitness[i];
            bestSol = pop[i];
        }
    }

    // -------------------------------------
    // LOA main loop
    // -------------------------------------
    for (int it = 1; it <= MAX_IT; it++)
    {

        for (int i = 0; i < POP; i++)
        {

            // Pick random better candidate
            int betterIdx = -1;

            vector<int> betterList;
            for (int j = 0; j < POP; j++)
                if (fitness[j] < fitness[i])
                    betterList.push_back(j);

            if (!betterList.empty())
                betterIdx = betterList[rng() % betterList.size()];

            // Choose escape or hide
            vector<double> candidate;
            double rp = randF(0, 1);

            if (rp < 0.5 && betterIdx != -1)
                candidate = escape(pop[i], pop[betterIdx]);
            else
                candidate = hide(pop[i], it);

            double f = sphere(candidate);

            // Accept if better
            if (f < fitness[i])
            {
                pop[i] = candidate;
                fitness[i] = f;

                if (f < bestFit)
                {
                    bestFit = f;
                    bestSol = candidate;
                }
            }
        }

        // OPTIONAL: Print progress
        cout << "Iteration " << it
             << " | Best = " << bestFit << endl;
    }

    // -------------------------------------
    // Final Output
    // -------------------------------------
    cout << "\nBest Solution Found:\n";
    for (int d = 0; d < DIM; d++)
        cout << "x" << d + 1 << " = " << bestSol[d] << endl;

    cout << "\nBest Sphere Value = " << bestFit << endl;

    return 0;
}

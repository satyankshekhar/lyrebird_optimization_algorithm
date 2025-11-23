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
// Rosenbrock Function
// ------------------------------
double rosenbrock(const vector<double> &x)
{
    int n = x.size();
    double s = 0.0;
    for (int i = 0; i < n - 1; ++i)
        s += 100.0 * pow(x[i + 1] - x[i] * x[i], 2) + pow(x[i] - 1.0, 2);
    return s;
}

// ------------------------------
// LOA PARAMETERS
// ------------------------------
int POP = 20;        // population size
int DIM = 5;         // dimension
int MAX_IT = 200;    // number of iterations
double LB = -30.0;   // lower bound
double UB = 30.0;    // upper bound

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
    auto pop = init_population();
    vector<double> fitness(POP);

    for (int i = 0; i < POP; i++)
        fitness[i] = rosenbrock(pop[i]);

    double bestFit = 1e18;
    vector<double> bestSol(DIM);

    for (int i = 0; i < POP; i++)
    {
        if (fitness[i] < bestFit)
        {
            bestFit = fitness[i];
            bestSol = pop[i];
        }
    }

    for (int it = 1; it <= MAX_IT; it++)
    {
        for (int i = 0; i < POP; i++)
        {
            int betterIdx = -1;

            vector<int> betterList;
            for (int j = 0; j < POP; j++)
                if (fitness[j] < fitness[i])
                    betterList.push_back(j);

            if (!betterList.empty())
                betterIdx = betterList[rng() % betterList.size()];

            vector<double> candidate;
            double rp = randF(0, 1);

            if (rp < 0.5 && betterIdx != -1)
                candidate = escape(pop[i], pop[betterIdx]);
            else
                candidate = hide(pop[i], it);

            double f = rosenbrock(candidate);

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

        cout << "Iteration " << it << " | Best = " << bestFit << endl;
    }

    cout << "\nBest Solution Found:\n";
    for (int d = 0; d < DIM; d++)
        cout << "x" << d + 1 << " = " << bestSol[d] << endl;

    cout << "\nBest Rosenbrock Value = " << bestFit << endl;

    return 0;
}

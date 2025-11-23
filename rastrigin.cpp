#include <bits/stdc++.h>
using namespace std;

mt19937 rng(time(NULL));
double randF(double a, double b){ uniform_real_distribution<double> d(a,b); return d(rng); }

double rastrigin(const vector<double>& x){
    int n = x.size();
    const double PI = acos(-1.0);
    double s = 10.0 * n;
    for (int i = 0; i < n; ++i)
        s += x[i]*x[i] - 10.0 * cos(2.0*PI*x[i]);
    return s;
}

int POP = 20, DIM = 5, MAX_IT = 200;
double LB = -5.12, UB = 5.12;

vector<vector<double>> init_population(){
    vector<vector<double>> pop(POP, vector<double>(DIM));
    for(int i=0;i<POP;i++) for(int d=0;d<DIM;d++) pop[i][d]=randF(LB,UB);
    return pop;
}

vector<double> escape(const vector<double>& x,const vector<double>& SSA){
    int dim = x.size(); vector<double> newX(dim);
    for(int j=0;j<dim;j++){ double r=randF(0,1); int I=(rng()%2)+1; newX[j]=x[j]+r*(SSA[j]-I*x[j]); if(newX[j]<LB)newX[j]=LB; if(newX[j]>UB)newX[j]=UB; }
    return newX;
}

vector<double> hide(const vector<double>& Xi,int t){ vector<double> newX=Xi; for(int j=0;j<DIM;j++){ double r=randF(0,1); newX[j]=Xi[j]+(1-2*r)*(UB-LB)/t; if(newX[j]<LB)newX[j]=LB; if(newX[j]>UB)newX[j]=UB;} return newX; }

int main(){
    auto pop = init_population(); vector<double> fitness(POP);
    for(int i=0;i<POP;i++) fitness[i]=rastrigin(pop[i]);
    double bestFit=1e18; vector<double> bestSol(DIM);
    for(int i=0;i<POP;i++) if(fitness[i]<bestFit){bestFit=fitness[i]; bestSol=pop[i];}
    for(int it=1; it<=MAX_IT; ++it){
        for(int i=0;i<POP;i++){
            int betterIdx=-1; vector<int> betterList;
            for(int j=0;j<POP;j++) if(fitness[j]<fitness[i]) betterList.push_back(j);
            if(!betterList.empty()) betterIdx=betterList[rng()%betterList.size()];
            vector<double> candidate; double rp=randF(0,1);
            if(rp<0.5 && betterIdx!=-1) candidate=escape(pop[i], pop[betterIdx]); else candidate=hide(pop[i], it);
            double f = rastrigin(candidate);
            if(f < fitness[i]){ pop[i]=candidate; fitness[i]=f; if(f<bestFit){ bestFit=f; bestSol=candidate; } }
        }
        cout<<"Iteration "<<it<<" | Best = "<<bestFit<<"\n";
    }
    cout<<"\nBest Solution Found:\n"; for(int d=0;d<DIM;d++) cout<<"x"<<d+1<<" = "<<bestSol[d]<<"\n";
    cout<<"\nBest Rastrigin Value = "<<bestFit<<"\n";
    return 0;
}

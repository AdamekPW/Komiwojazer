#include <iostream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <stack>
#include <limits>
#include <string>
#include <sstream>
#include <set>
#include <random>
#include <unordered_set>
#include <fstream>
#include <utility>
#include <numeric>
#include "SimpleAlgorithms.h"
#include "adam.h"

bool CmpFullGenome(const FullGenome &FG1, const FullGenome &FG2)
{
    return FG1.PathValue < FG2.PathValue;
}

void SortFullGenomes(vector<FullGenome> &FullGenomes)
{
    sort(FullGenomes.begin(), FullGenomes.end(), CmpFullGenome);
}
// do zrobienia: funkcja pomiaru czasu, optymalizacja wyboru najlepszych potomkow -
//  moze odbywac sie poprzez posortowanie vector<Genome> Genomes, a nastepnie resize()

vector<FullGenome> GenereteFullGenomes(vector<vector<double>> &Matrix, int HowMany)
{
    vector<FullGenome> FullGenomes;
    for (int i = 0; i < Matrix.size(); i++)
    {
        Genome G(Matrix.size(), -1);
        Przejdz3(Matrix, G, i);
        // Przejdz2(Matrix, G);
        double Cost = FindTotalCost(G, Matrix);
        FullGenomes.push_back({Cost, G});
    }
    SortFullGenomes(FullGenomes);
    return FullGenomes;
}

double FindBestGenome(vector<Genome> &Genomes, vector<vector<double>> &Matrix)
{
    double minCost = numeric_limits<int>::max();
    for (int i = 0; i < Genomes.size(); i++)
    {
        double k = FindTotalCost(Genomes[i], Matrix);
        if (k < minCost)
            minCost = k;
    }
    return minCost;
}

bool IsInVector(int k, vector<int> vec)
{
    for (int i = 0; i < vec.size(); i++)
        if (vec[i] == k)
            return true;
    return false;
}

pair<Genome, Genome> GenomesCrossover(Genome &P1, Genome &P2)
{
    int n = P1.size();
    int CuttingPoint = RandomNumber(2, n - 1);

    Genome K1, K2;
    K2.push_back(P2[0]);
    for (int i = 0; i < CuttingPoint; i++)
        K1.push_back(P1[i]);
    for (int i = CuttingPoint; i < n; i++)
        if (K2[0] != P1[i])
            K2.push_back(P1[i]);

    for (int i = 0; i < n; i++)
    {
        int Gen = P2[i];
        if (!IsInVector(Gen, K1))
        {
            K1.push_back(Gen);
        }
        if (!IsInVector(Gen, K2))
        {
            K2.push_back(Gen);
        }
    }

    return make_pair(K1, K2);
}

void MakeRandomMutation(Genome &G)
{
    int number = RandomNumber(0, 2);
    if (number == 0)
        Mutation1(G);
    if (number == 1)
        Mutation2(G);
}

double FullASGenetics(vector<vector<double>> Matrix, int ParentNumber, int CrossoversPerGeneration, int MutationsPerGeneration, int Iterations)
{
    vector<FullGenome> FullGenomes = GenereteFullGenomes(Matrix, ParentNumber);

    cout << "Full Na starcie: " << FullGenomes[0].PathValue << "\n";
    for (int generacja = 0; generacja < Iterations; generacja++)
    {
        for (int i = 0; i < CrossoversPerGeneration; i++)
        {
            int First = RandomNumber(0, FullGenomes.size());
            int Second = RandomNumber(0, FullGenomes.size());
            while (First == Second)
                Second = RandomNumber(0, FullGenomes.size());
            pair<Genome, Genome> Kids = GenomesCrossover(FullGenomes[First].G, FullGenomes[Second].G);
            FullGenomes.push_back({FindTotalCost(Kids.first, Matrix), Kids.first});
            FullGenomes.push_back({FindTotalCost(Kids.second, Matrix), Kids.second});
        }

        for (int i = 0; i < MutationsPerGeneration; i++)
        {
            int ToMutate = RandomNumber(0, FullGenomes.size());
            FullGenome FG = FullGenomes[ToMutate];
            FullGenomes.push_back(FG);
            MakeRandomMutation(FullGenomes[ToMutate].G);
            FullGenomes[ToMutate].PathValue = FindTotalCost(FullGenomes[ToMutate].G, Matrix);
        }

        SortFullGenomes(FullGenomes);
        FullGenomes.resize(ParentNumber);
    }

    return FullGenomes[0].PathValue;
}

#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include <stack>
#include <limits>
#include <string>
#include <sstream>
#include <set>
#include <random>
#include <unordered_set>
#include <fstream>
#include <ctime>

using namespace std;

struct Edge
{
    int V1, V2;
    double Distance;
};

typedef vector<int> Genome;

void PrintMatrix(vector<vector<double>> &Matrix);
void PrintMatrix(vector<vector<int>> &Matrix);

void PrintVector(vector<int> &vec);
void PrintVector(vector<double> &vec);

void PrintMST(vector<Edge> &Edges);

void SortEdgesByDistance(std::vector<Edge> &Edges);

vector<vector<double>> ConvertMatrix(vector<vector<double>> Matrix);

double FindTotalCost(vector<int> &Path, vector<vector<double>> &Matrix);

// zachlanny--------------
vector<vector<double>> BuildMST(vector<vector<double>> &Matrix);
void DFS(vector<vector<double>> &Matrix, stack<int> &Stos, vector<bool> &Visited, vector<int> &Wynik, double *Koszt);
double TSP(vector<vector<double>> &Matrix);

// heurystyka-------------------------------------
int RandomNumber(int LowerLimit, int UpperLimit);
int Przejdz(vector<vector<double>> &Matrix, vector<int> &genome, vector<bool> &Visited, double *Cost, int i, bool Randomize);
void Przejdz2(vector<vector<double>> &Matrix, vector<int> &Wynik);
void Przejdz3(vector<vector<double>> &Matrix, vector<int> &Wynik, int Start);
double NaiveNearestUnvisited(vector<vector<double>> &Matrix);
void wstawOdleglosc(vector<vector<double>> &matrix, int skad, int dokad, double wartosc);
double calcDistance(int x1, int y1, int x2, int y2);

double getDistance(vector<vector<double>> &matrix, int v_source, int v_dest);

double calcDistance(int x1, int y1, int x2, int y2);

double calcPathLength(vector<vector<double>> &matrix, vector<int> &path);

vector<vector<int>> GenereteGenomes(vector<vector<double>> &Matrix, int k);

void Swap(int *l1, int *l2);
void Mutation1(Genome &G);
void Mutation2(Genome &G);

vector<vector<int>> selectBestN(vector<vector<int>> &population, vector<vector<double>> &matrix, int n);

vector<vector<int>> selectBestN(vector<vector<int>> &vec, vector<double> &calculatedKeys, int n);
// dokladny ---------------------------------------

void BnB(vector<vector<double>> &Matrix, int Current, double CurrentCost,
         vector<bool> Visited, int VisitedCount, double *MinCost, int Start);

double Optimim(vector<vector<double>> &matrix);

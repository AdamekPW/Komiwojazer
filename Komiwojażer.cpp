#include <iostream>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <limits>
#include <string>
#include <sstream>
#include <fstream>
#include <utility>
#include "SimpleAlgorithms.h"
#include "adam.h"
#include "tomek.h"

using namespace std;

vector<vector<double>> loadDataFromFile(ifstream &Plik)
{
    int n;
    Plik >> n;
    int v;
    double x, y;
    vector<int> coor_x(n, 0);
    vector<int> coor_y(n, 0);
    for (int i = 0; i < n; i++)
    {
        Plik >> v >> x >> y;
        coor_x.at(v - 1) = x;
        coor_y.at(v - 1) = y;
    }
    vector<vector<double>> matrix(n, vector<double>(n, 0));

    for (int skad = 0; skad < n; skad++)
    {
        for (int dokad = 0; dokad < n; dokad++)
        {
            double odl = calcDistance(coor_x[skad], coor_y[skad], coor_x[dokad], coor_y[dokad]);
            wstawOdleglosc(matrix, skad, dokad, odl);
        }
    }

    for (int i = 0; i < matrix.size(); i++)
    {
        for (int j = 0; j < matrix[0].size(); j++)
        {
            matrix[j][i] = matrix[i][j];
        }
    }

    return matrix;
}

void RunTests()
{
    string Katalog = "Tests/";
    vector<string> Pliki = {"Berlin52.txt", "bier127.txt", "tsp250.txt", "tsp500.txt", "tsp1000.txt"};
    // vector<string> Pliki = {"Berlin52.txt"};
    for (auto x : Pliki)
    {
        cout << x << "\n";
        stringstream ss;
        ss << Katalog << x;
        string FileName = ss.str();
        ifstream plik(FileName);
        if (!plik.is_open())
        {
            cout << "Blad otwarcia\n";
            return;
        }
        // Berlin52
        vector<vector<double>> Matrix = loadDataFromFile(plik);
        int n = Matrix.size();

        pair<double, Genome> wynik = geneticAlgorithm(Matrix);
        cout << "Genetyczny: " << wynik.first << "\n";

        plik.close();
        cout << "------------------\n";
    }
}
void RunTest(string FileName)
{
    string Katalog = "Tests/";
    stringstream ss;
    ss << Katalog << FileName;
    ifstream plik(ss.str());
    if (!plik.is_open())
    {
        cout << "Blad otwarcia\n";
        return;
    }
    vector<vector<double>> Matrix = loadDataFromFile(plik);
    int n = Matrix.size();

    const int ileMinut = 3;
    auto start = chrono::high_resolution_clock::now();
    double durationSec;
    double best = numeric_limits<double>::max();

    while (durationSec < (ileMinut * 60))
    {
        pair<double, Genome> wynik = geneticAlgorithm(Matrix);
        auto end = chrono::high_resolution_clock::now();
        durationSec = chrono::duration_cast<chrono::duration<double>>(end - start).count();
        // cout << durationSec << " Wynik: " << wynik.first << "\n";
        if (wynik.first < best)
        {
            best = wynik.first;
            printf("Czas:\t%.2f\tWynik:\t%.2f\tBEST\n", durationSec, wynik.first);
        }
        else
        {
            printf("Czas:\t%.2f\tWynik:\t%.2f\n", durationSec, wynik.first);
        }
    }
    printf("All best\t%.2f\n", best);

    // double wynik = NaiveNearestUnvisited(Matrix);
    // cout << "Wynik: " << wynik << "\n";

    plik.close();
}

int main(int argc, char *argv[])
{
    srand(time(NULL));
    if (argc == 1)
    {
        RunTests();
    }
    else if (argc == 2)
    {
        cout << "Odpalanie na pliku " << argv[1] << "\n";
        RunTest(argv[1]);
    }
}
#include <limits>
#include "SimpleAlgorithms.h"

void PrintMatrix(vector<vector<double>> &Matrix)
{
    for (int i = 0; i < Matrix.size(); i++)
    {
        for (int j = 0; j < Matrix[i].size(); j++)
        {
            printf(" %.2f ", Matrix[i][j]);
        }
        cout << "\n";
    }
};
void PrintMatrix(vector<vector<int>> &Matrix)
{
    for (int i = 0; i < Matrix.size(); i++)
    {
        for (int j = 0; j < Matrix[i].size(); j++)
        {
            cout << Matrix[i][j] << " ";
        }
        cout << "\n";
    }
};
void PrintVector(vector<int> &vec)
{
    for (int i = 0; i < vec.size(); i++)
        cout << vec[i] << " ";
    cout << "\n";
};
void PrintVector(vector<double> &vec)
{
    for (int i = 0; i < vec.size(); i++)
        cout << vec[i] << " ";
    cout << "\n";
};

void PrintMST(vector<Edge> &Edges)
{
    for (int i = 0; i < Edges.size(); i++)
    {
        cout << Edges[i].V1 << "-" << Edges[i].V2 << " | " << Edges[i].Distance << "\n";
    }
}

void SortEdgesByDistance(std::vector<Edge> &Edges)
{
    std::sort(Edges.begin(), Edges.end(), [](const Edge &a, const Edge &b)
              { return a.Distance < b.Distance; });
}

vector<vector<double>> ConvertMatrix(vector<vector<double>> Matrix)
{
    int n = Matrix.size();
    vector<vector<double>> M1(n - 1, vector<double>(n - 1, 0));

    for (int i = 1; i < n; i++)
    {
        for (int j = 1; j < n; j++)
        {
            M1[i - 1][j - 1] = Matrix[i][j];
        }
    }
    for (int i = 0; i < M1.size(); i++)
    {
        for (int j = 0; j < M1[0].size(); j++)
        {
            M1[j][i] = M1[i][j];
        }
    }
    return M1;
}

double FindTotalCost(vector<int> &Path, vector<vector<double>> &Matrix)
{
    // Dziala na macierzy Tomasza (od V = 1)
    double Cost = 0;
    for (int i = 0; i < Path.size() - 1; i++)
    {
        int From = Path[i];
        int To = Path[i + 1];
        if (Matrix[From][To] != 0)
            Cost += Matrix[From][To];
        else
            Cost += Matrix[To][From];
    }

    int From = Path[Path.size() - 1];
    int To = Path[0];
    if (Matrix[From][To] != 0)
        Cost += Matrix[From][To];
    else
        Cost += Matrix[To][From];
    return Cost;
}

// zachlanny----------------------------------------------------------------------------
vector<vector<double>> BuildMST(vector<vector<double>> &Matrix)
{

    vector<vector<double>> FinalMatrix(Matrix.size(), vector<double>(Matrix.size(), 0));
    vector<Edge> Edges;
    for (int i = 0; i < Matrix.size(); i++)
    {
        for (int j = i; j < Matrix[i].size(); j++)
        {
            if (Matrix[i][j] != 0)
                Edges.push_back({i, j, Matrix[i][j]});
        }
    }

    // Sortowanie kraw?dzi po wagach
    SortEdgesByDistance(Edges);

    // tworzenie zbiorow, przy czym na poczatku kazdy zbior ma tylko jeden wierzcholek
    vector<int> ListOfRepresentatives;
    for (int i = 0; i < Matrix.size(); i++)
    {
        ListOfRepresentatives.push_back(i); // chwilowo kazdy wierzcholek jest swoim wlasnym reprezentantem
    }

    int i = 0;
    int HowManySetsLeft = Matrix.size();
    bool End = false;
    while (i < Edges.size() && !End)
    {

        int V1 = Edges[i].V1;
        int V2 = Edges[i].V2;
        if (ListOfRepresentatives[V1] != ListOfRepresentatives[V2])
        {
            // maja roznych reprezentantow, mozna dodac
            FinalMatrix[V1][V2] = Edges[i].Distance;
            FinalMatrix[V2][V1] = Edges[i].Distance;

            int NowyReprezentant = ListOfRepresentatives[V1];
            int StaryReprezentant = ListOfRepresentatives[V2];
            if (ListOfRepresentatives[V2] < ListOfRepresentatives[V1])
            {
                NowyReprezentant = ListOfRepresentatives[V2];
                StaryReprezentant = ListOfRepresentatives[V1];
            }
            // cout << "Krawedz: " << StaryReprezentant << "-" << NowyReprezentant << "\n";
            for (int j = 0; j < ListOfRepresentatives.size(); j++)
            {

                if (ListOfRepresentatives[j] == StaryReprezentant)
                {
                    ListOfRepresentatives[j] = NowyReprezentant;
                    // cout << "Podmieniono: " << j << " na " << NowyReprezentant << '\n';
                }
            }

            // PrintVector(ListOfRepresentatives);
            int j = 0;
            while (j < ListOfRepresentatives.size() - 1 && ListOfRepresentatives[j] == ListOfRepresentatives[j + 1])
                j++;
            if (j == ListOfRepresentatives.size() - 1)
                End = true;
        }
        i++;
    }
    // PrintVector(ListOfRepresentatives);
    // cout << "I = " << i << " | " << HowManySetsLeft << "\n";

    // cout << "---------\n";
    // PrintMatrix(FinalMatrix);

    // cout << "---------\n";
    return FinalMatrix;
}

void DFS(vector<vector<double>> &Matrix, stack<int> &Stos, vector<bool> &Visited, vector<int> &Wynik, double *Koszt)
{

    int i = Stos.top();
    // PrintVector(Matrix[i]);
    // cout << "Wierzcholek: " << i << "\n";
    int j = 0;
    while (j < Matrix[i].size() && (Matrix[i][j] <= 0.0f || (Visited[j] && Matrix[i][j] > 0)))
        j++;

    if (j == Matrix[i].size())
    {
        if (!Stos.empty())
        {

            Wynik.push_back(Stos.top());
            int V1 = Stos.top();
            Stos.pop();
            if (!Stos.empty())
            {
                int V2 = Stos.top();
                (*Koszt) += Matrix[V1][V2];
                DFS(Matrix, Stos, Visited, Wynik, Koszt);
            }
        }
        else
        {
            return;
        }
    }
    else
    {
        (*Koszt) += Matrix[i][j];
        Stos.push(j);

        Visited[j] = true;

        DFS(Matrix, Stos, Visited, Wynik, Koszt);
    }
}

double TSP(vector<vector<double>> &Matrix)
{
    vector<vector<double>> MatrixMST = BuildMST(Matrix);
    stack<int> Stos;
    vector<int> Wynik;
    vector<bool> Visited(MatrixMST.size(), false);
    Stos.push(0);
    Visited[0] = true;
    double Cost = 0;
    DFS(MatrixMST, Stos, Visited, Wynik, &Cost);
    // PrintVector(Wynik);

    return Cost;
}

int RandomNumber(int LowerLimit, int UpperLimit)
{
    int losowaLiczba = rand() % (UpperLimit - LowerLimit) + LowerLimit;
    return losowaLiczba;
}

int Przejdz(vector<vector<double>> &Matrix, vector<int> &genome, vector<bool> &Visited, double *Cost, int i, bool Randomize)
{
    // cout << i << "\n";
    genome.push_back(i);
    Visited[i] = true;
    int min = numeric_limits<int>::max();
    int Index = -1;
    for (int j = 0; j < Matrix.size(); j++)
    {
        if (Matrix[i][j] > 0 && Matrix[i][j] <= min && !Visited[j])
        {
            min = Matrix[i][j];
            Index = j;
        }
    }
    if (Randomize)
    {
        int SecondMin = numeric_limits<int>::max();
        int SecondIndex = -1;
        for (int j = 0; j < Matrix.size(); j++)
        {

            if (Matrix[i][j] > 0 && !Visited[j] && Matrix[i][j] <= SecondMin && j != Index)
            {
                SecondMin = -1;
                SecondIndex = j;
            }
        }
        if (SecondIndex != -1)
            Index = SecondIndex;
    }

    if (Index != -1)
    {
        (*Cost) += Matrix[i][Index];
        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<> dis(0, 1);
        Randomize = dis(gen);
        return Przejdz(Matrix, genome, Visited, Cost, Index, Randomize);
    }
    return i;
}

void Przejdz2(vector<vector<double>> &Matrix, vector<int> &Wynik)
{
    int n = Matrix.size();
    int Start = RandomNumber(0, n);
    // vector<int> Wynik(n, -1);
    vector<bool> Visited(n, false);
    Visited[Start] = true;
    int VisitedCount = 1;

    int i = Start;
    while (VisitedCount <= n)
    {
        Wynik[VisitedCount - 1] = i;
        Visited[i] = true;

        double min1Value = numeric_limits<int>::max();
        int min1Index = -1;
        double min2Value = numeric_limits<int>::max();
        int min2Index = -1;

        for (int j = 0; j < n; j++)
        {
            if (Matrix[i][j] > 0 && !Visited[j] && Matrix[i][j] <= min1Value)
            {
                min2Value = min1Value;
                min2Index = min1Index;
                min1Value = Matrix[i][j];
                min1Index = j;
            }
        }
        int Randomize = RandomNumber(0, 2);
        if (Randomize == 0 && min2Index != -1)
            i = min2Index;
        else
            i = min1Index;
        VisitedCount++;
    }
}

void Przejdz3(vector<vector<double>> &Matrix, vector<int> &Wynik, int Start)
{
    int n = Matrix.size();
    // int Start = RandomNumber(0, n);
    // vector<int> Wynik(n, -1);
    vector<bool> Visited(n, false);
    Visited[Start] = true;
    int VisitedCount = 1;

    int i = Start;
    while (VisitedCount <= n)
    {
        Wynik[VisitedCount - 1] = i;
        Visited[i] = true;

        double min1Value = numeric_limits<int>::max();
        int min1Index = -1;

        for (int j = 0; j < n; j++)
        {
            if (Matrix[i][j] > 0 && !Visited[j] && Matrix[i][j] <= min1Value)
            {
                min1Value = Matrix[i][j];
                min1Index = j;
            }
        }

        i = min1Index;
        VisitedCount++;
    }
}

double NaiveNearestUnvisited(vector<vector<double>> &Matrix)
{
    int size = Matrix.size();
    mt19937 gen = mt19937{random_device{}()};
    uniform_int_distribution<> randomVertex(0, size);

    vector<bool> visited(Matrix.size(), false);
    int initialVertex = randomVertex(gen);
    int v = initialVertex;
    int next;
    double distance;
    double minEdgeWeight = numeric_limits<double>::max();
    double cost = 0;

    while (!visited[v])
    {
        for (int i = 0; i < Matrix.size(); i++)
        {
            if (v <= i)
            {
                distance = Matrix[v][i];
            }
            else
            {
                distance = Matrix[i][v];
            }
            if (i == v || visited[i] || distance == 0)
            {
                continue;
            }
            if (distance < minEdgeWeight)
            {
                next = i;
                minEdgeWeight = distance;
            }
        }
        if (v == next)
        {
            break;
        }
        cost += minEdgeWeight;
        visited[v] = true;
        v = next;
        minEdgeWeight = numeric_limits<int>::max();
    }
    cost += getDistance(Matrix, v, initialVertex);
    return cost;
}
void wstawOdleglosc(vector<vector<double>> &matrix, int skad, int dokad, double wartosc)
{
    if (skad <= dokad)
        matrix[skad][dokad] = wartosc;
    else
        matrix[dokad][skad] = wartosc;
}

double getDistance(vector<vector<double>> &matrix, int v_source, int v_dest)
{
    if (v_source == v_dest)
    {
        cout << "Co to za chodzenie do siebie?!" << '\n';
        exit(1);
    }
    return (v_source < v_dest) ? matrix[v_source][v_dest] : matrix[v_dest][v_source];
}

double calcDistance(int x1, int y1, int x2, int y2)
{
    return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
}

double calcPathLength(vector<vector<double>> &matrix, vector<int> &path)
{
    double sum = 0;
    for (int i{0}; i < path.size() - 1; i++)
    {
        int v_source = path.at(i);
        int v_dest = path.at(i + 1);
        sum += getDistance(matrix, v_source, v_dest);
    }
    sum += getDistance(matrix, path.at(path.size() - 1), path.at(0));
    return sum;
}

vector<vector<int>> GenereteGenomes(vector<vector<double>> &Matrix, int k)
{
    vector<vector<int>> Genomes;
    for (int i = 0; i < Matrix.size(); i++)
    {
        vector<int> G(Matrix.size(), -1);
        Przejdz3(Matrix, G, i);
        Genomes.push_back(G);
    }
    vector<vector<int>> bestGenomes = selectBestN(Genomes, Matrix, k);
    return bestGenomes;
}

void Swap(int *l1, int *l2)
{
    int tmp = *l1;
    *l1 = *l2;
    *l2 = tmp;
}

void Mutation1(Genome &G)
{
    // Mutacja przez zamiane dwoch losowo wybranych miast
    int n = G.size();
    int FirstGen = RandomNumber(1, n);
    int SecondGen = RandomNumber(1, n);
    while (FirstGen == SecondGen)
        SecondGen = RandomNumber(1, n);

    Swap(&G[FirstGen], &G[SecondGen]);
}

void Mutation2(Genome &G)
{
    // Mutacja przez wyciecie ciaglu miast i wstawienie w losowe miejsce
    int n = G.size();
    int SizeOfCut = RandomNumber(1, (int)(n / 2));
    int StartingElement = RandomNumber(1, n - SizeOfCut);
    int EndingElement = StartingElement + SizeOfCut;

    Genome G1;
    G1.push_back(G[0]);

    int NewPosition = 0;
    int FirstPartSize = StartingElement;
    int SecondPartSize = n - EndingElement;
    int Proportion = FirstPartSize + SecondPartSize;
    int FirstProportion = (int)(100 * (float)FirstPartSize / (float)Proportion);
    int SecondProportion = (int)(100 * (float)SecondPartSize / (float)Proportion);

    int FromWhichPart = RandomNumber(0, 100);
    if (FromWhichPart <= FirstProportion)
    {
        if (StartingElement != 1)
            NewPosition = RandomNumber(1, StartingElement);
        else
            NewPosition = 1;

        for (int i = 1; i < NewPosition; i++)
            G1.push_back(G[i]);
        for (int i = StartingElement; i < EndingElement; i++)
            G1.push_back(G[i]);
        for (int i = NewPosition; i < StartingElement; i++)
            G1.push_back(G[i]);
        for (int i = EndingElement; i < n; i++)
            G1.push_back(G[i]);
    }
    else
    {
        NewPosition = RandomNumber(EndingElement, n + 1);

        for (int i = 1; i < StartingElement; i++)
            G1.push_back(G[i]);
        for (int i = EndingElement; i < NewPosition; i++)
            G1.push_back(G[i]);
        for (int i = StartingElement; i < EndingElement; i++)
            G1.push_back(G[i]);
        for (int i = NewPosition; i < n; i++)
            G1.push_back(G[i]);
    }

    // cout << "Cutting position: " << NewPosition << "\n";

    for (int i = 0; i < n; i++)
        G[i] = G1[i];
}

vector<vector<int>> selectBestN(vector<vector<int>> &population, vector<vector<double>> &matrix, int n)
{
    nth_element(population.begin(), population.begin() + n, population.end(),
                [&matrix](vector<int> first, vector<int> second)
                { return calcPathLength(matrix, first) < calcPathLength(matrix, second); });

    vector<vector<int>> populationTrimmed;
    copy(population.begin(), population.begin() + n, back_inserter(populationTrimmed));
    return populationTrimmed;
}

vector<vector<int>> selectBestN(vector<vector<int>> &vec, vector<double> &calculatedKeys, int n)
{
    int vecSize = vec.size();
    vector<pair<double, vector<int>>> vectorsToBeSorted;
    for (int i{0}; i < vecSize; i++)
    {
        pair<double, vector<int>> vecWithSortingKey = make_pair(calculatedKeys[i], vec[i]);
        vectorsToBeSorted.push_back(vecWithSortingKey);
    }

    nth_element(vectorsToBeSorted.begin(), vectorsToBeSorted.begin() + n, vectorsToBeSorted.end(),
                [](pair<double, vector<int>> &first, pair<double, vector<int>> &second)
                { return first.first < second.first; });

    vector<vector<int>> result;
    for (int i{0}; i < n; i++)
    {
        result.push_back(vectorsToBeSorted[i].second);
    }
    return result;
}

// dokladny-----------------------------------------------------------------------------

void BnB(vector<vector<double>> &Matrix, int Current, double CurrentCost,
         vector<bool> Visited, int VisitedCount, double *MinCost, int Start)
{

    if (VisitedCount == Matrix.size() && Current == Start && CurrentCost <= *MinCost)
    {
        *MinCost = CurrentCost;
        cout << "New Value = " << CurrentCost << "\n";
        return;
    }
    for (int j = 0; j < Matrix.size(); j++)
    {
        if (Matrix[Current][j] > 0 && !Visited[j] && CurrentCost + Matrix[Current][j] <= *MinCost)
        {
            Visited[j] = true;
            VisitedCount++;
            BnB(Matrix, j, CurrentCost + Matrix[Current][j], Visited, VisitedCount, MinCost, Start);
            Visited[j] = false;
            VisitedCount--;
        }
    }
}

double Optimim(vector<vector<double>> &matrix)
{
    double MinCost = NaiveNearestUnvisited(matrix);
    cout << "Optimum | Z heurystyki: " << MinCost << "\n";
    vector<vector<double>> Matrix = ConvertMatrix(matrix);
    vector<bool> Visited(Matrix.size(), false);
    // Visited[0] = true;
    // MinCost = 200000;
    int n = Matrix.size() - 1;
    BnB(Matrix, n, 0, Visited, 0, &MinCost, n);

    return MinCost;
}

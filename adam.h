#include <vector>

using namespace std;

typedef vector<int> Genome;

struct FullGenome
{
    double PathValue;
    Genome G;
};

bool CmpFullGenome(const FullGenome &FG1, const FullGenome &FG2);
void SortFullGenomes(vector<FullGenome> &FullGenomes);
vector<FullGenome> GenereteFullGenomes(vector<vector<double>> &Matrix, int HowMany);
double FindBestGenome(vector<Genome> &Genomes, vector<vector<double>> &Matrix);
bool IsInVector(int k, vector<int> vec);
pair<Genome, Genome> GenomesCrossover(Genome &P1, Genome &P2);
void MakeRandomMutation(Genome &G);
double FullASGenetics(vector<vector<double>> Matrix, int ParentNumber, int CrossoversPerGeneration, int MutationsPerGeneration, int Iterations);

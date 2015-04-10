#include <istream>
#include <string>
#include <sstream>
#include <vector>
#include <omp.h>
#include <fstream>
#include <iostream>
extern "C"
{
       #include <cblas.h>
}

using namespace std;


// load matrix from an ascii text file.
void load_matrix(std::istream* is,
        std::vector< std::vector<double> >* matrix,
        const std::string& delim = " \t")
{
    using namespace std;

    string      line;
    string      strnum;

    // clear first
    matrix->clear();

    // parse line by line
    while (getline(*is, line))
    {
        matrix->push_back(vector<double>());

        for (string::const_iterator i = line.begin(); i != line.end(); ++ i)
        {
            // If i is not a delim, then append it to strnum
            if (delim.find(*i) == string::npos)
            {
                strnum += *i;
                if (i + 1 != line.end()) // If it's the last char, do not continue
                    continue;
            }

            // if strnum is still empty, it means the previous char is also a
            // delim (several delims appear together). Ignore this char.
            if (strnum.empty())
                continue;

            // If we reach here, we got a number. Convert it to double.
            double       number;

            istringstream(strnum) >> number;
            matrix->back().push_back(number);

            strnum.clear();
        }
    }
}

void columnG_fast(int p, double * G_small, double * w, int Nsen, int Nsrc, double * G_col)
{
  int q = p % (Nsrc * Nsrc);
  int half = (p - 1) / (Nsrc * Nsrc); /* We have columns of two types: those which end with zeros and those which start with zeros. half determines the type*/
  int left = 0, right = 1;
  if(!q)
    q = Nsrc * Nsrc;
  int i = q % Nsrc;
  if(!i)
    i = Nsrc;
  i --;
  int j = (q - i) / Nsrc;
  double w_p = w[p-1];
  int k, l;
  for (l = 0; l < Nsen; ++l)
    for (k = 0; k < Nsen; ++k)
    {
      if(half == left)
      {
        G_col[k + Nsen * l] = G_small[Nsrc * k + i] * G_small[Nsrc * l + j] * w_p;
        G_col[Nsen * Nsen + k + Nsen * l] = 0.;
      }
      else if(half == right)
      {
        G_col[k + Nsen * l] = 0.;
        G_col[Nsen * Nsen + k + Nsen * l] = G_small[Nsrc * k + i] * G_small[Nsrc * l + j] * w_p;
      }
    }
}

void calc_violations(double * G_small, double * W, double * R, double * violations, int Ch, int Sr, int T)
{
    using namespace std;
    int Nsrc = Sr * Sr * 2;
    int Nsen = Ch * Ch * 2; /* Number of elements in a column */
    int s, i;
    double norm = 0.;
    double result[T];
    double * column = new double[Nsen];
     #pragma omp parallel for  /*ordered */
    for (s = 0; s < Nsrc; ++s)
    {
      columnG_fast(s+1, G_small, W, Ch, Sr, column);
      cblas_dgemv(CblasRowMajor, CblasTrans, Nsen, T, 1, R, T, column, 1, 0., result, 1);
      norm = cblas_dnrm2(T, result, 1);
      violations[s] = norm;
      if(!(s % 10000))
        cout << "s = "<< s << endl;
    }
}

// example


int main()
{  
    // initialize params
    using namespace std;
    // read the file
    std::ifstream Gin("G_small.txt");
    std::ifstream Win("w.txt");
    std::ifstream Rin("R.txt");

    // load matrices and sizes
    std::vector< std::vector<double> > G_v;
    std::vector< std::vector<double> > W_v;
    std::vector< std::vector<double> > R_v;
    load_matrix(&Gin, &G_v);
    load_matrix(&Rin, &R_v);
    load_matrix(&Win, &W_v);
    Gin.close();
    Win.close();
    Rin.close();
    int Src = G_v[0].size();
    int Ch = G_v.size();
    long int Nsrc = 2 * Src * Src;
    int Nch = 2 * Ch* Ch;
    int T = R_v[0].size();
    // ---------------------- //
    cout << "G Nraws = " << Ch << endl;
    cout << "G Ncolumns = " << Src << endl;
    cout << "T = " << T << endl;
    // Iitialize matrices //
    int i,j;
    double * G= new double[Ch * Src];
    for (i = 0; i < G_v.size(); i++)
        for (j = 0; j < G_v[0].size(); ++j)
            G[Src * i + j] = G_v[i][j];
    double * R = new double  [Nch * T];
    for (i = 0; i < R_v.size(); i++)
        for (j = 0; j < R_v[0].size(); ++j)
            R[T * i + j] = R_v[i][j];
    double * W = new double[Nsrc];
    for (i = 0; i < W_v.size(); i++)
        W[i] = W_v[i][0];
    double * V = new double[Nsrc];
    // ------------------------ //
    calc_violations(G, W, R, V, Ch, Src, T); 
     // cout << "V:\n";

    ofstream Vout("V.txt");
     for (int i = 0; i < Nsrc; ++i)
         Vout << V[i] <<endl;
    Vout.close();
    cout << "G Nraws = " << Ch << endl;
    cout << "G Ncolumns = " << Src << endl;
    cout << "T = " << T << endl;
    return 0;
}

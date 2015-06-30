#include <istream>
#include <string>
#include <sstream>
#include <vector>
#include <omp.h>
#include <fstream>
#include <iostream>
#include <math.h>
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



inline void Get_i_j_from_s(int s, int N, int & i, int & j)
{
  i = s % N;
  if(!i)
    i = N;
  i --;
  j = (s - i) / N;
}

void calc_violations(double * G_small, double * R, double * violations, int Ch, int Sr, int T)
{
    cout << "   Parallel section";
    #pragma omp parallel num_threads(8)
    {
      int Nsrc_pairs = Sr * Sr;
      int Nsites = Sr / 2;
      int Nsite_pairs = Nsites * (Nsites + 1) / 2; /* For each time instant this is number 
of source-space cross-spectrum elements above the diagonal plus on the diagonal itself. */
      int Nsen_pairs = Ch * Ch; /* Number of elements in a column */
      int s, t, k, l, m; // iterators
      int i, j; // indices
      double norm_sq = 0.;
      double result[T];
      double prod;
      double * temp = new double[Ch];
      double ** colRmat = new double * [T];

      s = 0;
      int IND[Nsite_pairs][2];
      for (i = 0; i < Nsites; ++i)
        for (j = i; j < Nsites; ++j)
        {
            IND[s][0] = i;
            IND[s][1] = j;
            s++;
        }

      for (t = 0; t < T; ++t)   
      {
        colRmat[t] = new double[Ch*Ch];
        s = 0;
        for (k = 0; k < Ch; ++k)
          for (l = k; l < Ch; ++l)
          {
            colRmat[t][Ch * k + l] = R[T * s + t];
            s++;
          }
      }

      double ** G = new double * [Sr];
      for (s = 0; s < Sr; ++s)
      { 
        G[s] = new double[Ch];
        for (m = 0; m < Ch; ++m)
            G[s][m] = G_small[Sr * m + s];
      }
      #pragma omp for 
      for (s = 0; s < Nsite_pairs; ++s)
      {     
        if(!(s % 800000))
          {
            printf(".");
            fflush(stdout);
          }
        // ---------------------- //
        // Get_i_j_from_s(s+1, Nsites, i, j);
        i = IND[s][0];
        j = IND[s][1];
        // cout << "i = " << i << " j = " << j << endl;
        
        for (int t = 0; t < T; ++t)
        {
          temp = G[2 * i];
          cblas_dtrmv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit, Ch, colRmat[t], Ch, temp, 1);
          prod = cblas_ddot(Ch, G[2 * j], 1, temp, 1);
          norm_sq += prod*prod;
          temp = G[2 * i + 1];
          cblas_dtrmv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit, Ch, colRmat[t], Ch, temp, 1);
          prod = cblas_ddot(Ch, G[2 * j], 1, temp, 1);
          norm_sq += prod*prod;
          temp = G[2 * i];
          cblas_dtrmv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit, Ch, colRmat[t], Ch, temp, 1);
          prod = cblas_ddot(Ch, G[2 * j + 1], 1, temp, 1);
          norm_sq += prod*prod;
          temp = G[2 * i + 1];
          cblas_dtrmv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit, Ch, colRmat[t], Ch, temp, 1);
          prod = cblas_ddot(Ch, G[2 * j + 1], 1, temp, 1);
          norm_sq += prod*prod;
        }
        violations[i * Nsites + j] = sqrt(norm_sq);
        norm_sq = 0;
        // ---------------------- //
      }
    }
}

// example


int main()
{  
    // initialize params
    using namespace std;
    // read the file
    std::ifstream Gin("../aux/G_small.txt");
    std::ifstream Rin("../aux/R.txt");

    // load matrices and sizes
    cout << "   Loading data..." << endl;
    std::vector< std::vector<double> > G_v;
    std::vector< std::vector<double> > R_v;
    load_matrix(&Gin, &G_v);
    load_matrix(&Rin, &R_v);
    Gin.close();
    Rin.close();
    // cout << "   Done." << endl;
    int Src = G_v[0].size();
    int Ch = G_v.size();
    long int Nsrc_pairs = Src * Src;
    int Nsites = Src / 2;
    int Nsite_pairs = Nsites * Nsites;
    int Nch = Ch* Ch;
    int T = R_v[0].size();
    // ---------------------- //
    // cout << "   G Nraws = " << Ch << endl;
    // cout << "   G Ncolumns = " << Src << endl;
    // cout << "   T = " << T << endl;
    // Iitialize matrices //
    int i,j;
    double * G = new double[Ch * Src];
    for (i = 0; i < G_v.size(); i++)
        for (j = 0; j < G_v[0].size(); ++j)
            G[Src * i + j] = G_v[i][j];
    double * R = new double  [Nch * T];
    for (i = 0; i < R_v.size(); i++)
        for (j = 0; j < R_v[0].size(); ++j)
            R[T * i + j] = R_v[i][j];
    double * V = new double[Nsite_pairs];
    // ------------------------ //
    calc_violations(G, R, V, Ch, Src, T); 
     // cout << "V:\n";
    cout << endl << "   Writing data...";
    ofstream Vout("../aux/V.txt");
     for (int i = 0; i < Nsite_pairs; ++i)
         Vout << V[i] <<endl;
    Vout.close();
    cout << "\n";
    return 0;
}

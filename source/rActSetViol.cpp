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

void columnG_fast(int iDi, int iDj, double * dmG_small, double w, int iCh, int iSr, double* dvG_col);
void calcGpair(int iSi, int iSj, double* dmG_small, double w, int iCh, int iSr, double* G_s);
void load_matrix(std::istream* is, std::vector< std::vector<double> >* matrix, const std::string& delim = " \t");

// load matrix from an ascii text file.


void CpGcolToGs(int iColG, int ixs, double* Gs, double* Gcol)
{
  int cS;
  for (cS = 0; cS < iColG; ++cS)
  {
    Gs[ixs] = Gcol[cS];
    ixs++;
  }
}


void calcGpair(int iSi, int iSj, double* dmG_small, double w, int iCh, int iSr, double* G_s)
{
  int iColG = iCh * iCh;
  int ixs = 0;
  int cS;
  double* dvG_col = new double[iColG];
  columnG_fast(iSi * 2, iSj * 2, dmG_small, w, iCh, iSr, dvG_col);
  CpGcolToGs(iColG, ixs, G_s, dvG_col);
  columnG_fast(iSi * 2 + 1, iSj * 2, dmG_small, w, iCh, iSr, dvG_col);
  CpGcolToGs(iColG, ixs, G_s, dvG_col);
  columnG_fast(iSi * 2 , iSj * 2 + 1, dmG_small, w, iCh, iSr, dvG_col);
  CpGcolToGs(iColG, ixs, G_s, dvG_col);
  columnG_fast(iSi * 2 + 1, iSj * 2 + 1, dmG_small, w, iCh, iSr, dvG_col);
  CpGcolToGs(iColG, ixs, G_s, dvG_col);
}

void columnG_fast(int iDi, int iDj, double * dmG_small, double w, int iCh, int iSr, double* dvG_col)
{
  int ixk, ixl;
  for (ixk = 0; ixk < iCh; ++ixk)
    for (ixl = 0; ixl < iCh; ++ixl)
      dvG_col[iCh * ixk +  ixl] = dmG_small[iSr * ixk + iDi] * dmG_small[iSr * ixl + iDj] * w;
}


inline void Get_i_j_from_s(int s, int N, int & i, int & j)
{
  i = s % N;
  if(!i)
    i = N;
  i --;
  j = (s - i) / N;
}

void calc_violations(double * dmG_small, double * W, double * R, double * violations, int iCh, int iSr, int T)
{
    cout << "   Parallel section";
    #pragma omp parallel num_threads(8)
    {
      int Nsrc_pairs = iSr * iSr;
      int Nsites = iSr / 2;
      int Nsite_pairs = Nsites * (Nsites + 1) / 2; /* For each time instant this is number 
of source-space cross-spectrum elements above the diagonal plus on the diagonal itself. */
      int colG = iCh * iCh; /* Number of elements in a column */
      int s, t, k, l, m; // iterators
      int i, j; // indices
      double norm_sq = 0.;
      double result[T];
      double prod;
      double * dmTemp = new double[4 * T];
     
      s = 0;
      int IND[Nsite_pairs][2];
      for (i = 0; i < Nsites; ++i)
        for (j = i; j < Nsites; ++j)
        {
            IND[s][0] = i;
            IND[s][1] = j;
            s++;
        }

      double* G_s = new double[4 * colG];
      
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
        calcGpair(i, j, dmG_small, W[s], iCh, iSr, G_s);
/*        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,\
4, T, colG, 1., G_s, colG, R, T, 0., dmTemp, T);
*/        
        // violations[i * Nsites + j] = cblas_dnrm2(4 * T, dmTemp, 1);
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
    std::ifstream Win("../aux/w.txt");
    std::ifstream Rin("../aux/R.txt");

    // load matrices and sizes
    cout << "   Loading data..." << endl;
    std::vector< std::vector<double> > G_v;
    std::vector< std::vector<double> > W_v;
    std::vector< std::vector<double> > R_v;
    load_matrix(&Gin, &G_v);
    load_matrix(&Rin, &R_v);
    load_matrix(&Win, &W_v);
    Gin.close();
    Win.close();
    Rin.close();
    // cout << "   Done." << endl;
    int Src = G_v[0].size();
    int iCh = G_v.size();
    long int Nsrc_pairs = Src * Src;
    int Nsites = Src / 2;
    int Nsite_pairs = Nsites * Nsites;
    int Nch = iCh* iCh;
    int T = R_v[0].size();
    // ---------------------- //
    // cout << "   G Nraws = " << iCh << endl;
    // cout << "   G Ncolumns = " << Src << endl;
    // cout << "   T = " << T << endl;
    // Iitialize matrices //
    int i,j;
    double * G = new double[iCh * Src];
    for (i = 0; i < G_v.size(); i++)
        for (j = 0; j < G_v[0].size(); ++j)
            G[Src * i + j] = G_v[i][j];
    double * R = new double  [Nch * T];
    for (i = 0; i < R_v.size(); i++)
        for (j = 0; j < R_v[0].size(); ++j)
            R[T * i + j] = R_v[i][j];
    double * W = new double[Nsrc_pairs];
    for (i = 0; i < W_v.size(); i++)
        W[i] = W_v[i][0];
    double * V = new double[Nsite_pairs];
    // ------------------------ //
    calc_violations(G, W, R, V, iCh, Src, T); 
     // cout << "V:\n";
    cout << endl << "   Writing data..." << endl;
    ofstream Vout("../aux/V.txt");
     for (int i = 0; i < Nsite_pairs; ++i)
         Vout << V[i] <<endl;
    Vout.close();
    cout << "\n";
    return 0;
}


void load_matrix(std::istream* is,
        std::vector< std::vector<double> >* matrix,
        const std::string& delim)
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


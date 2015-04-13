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

inline void columnG_fast(int p, double * G_small, double * w, int Nsen, int Nsrc, double * G_col)
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

inline void Get_i_j_from_s(int s, int Nsrc, int & i, int & j, int & half)
{
  int q = s % (Nsrc * Nsrc);
  half = (s - 1) / (Nsrc * Nsrc); /* We have columns of two types: those which end with zeros and those which start with zeros. half determines the type*/
  if(!q)
    q = Nsrc * Nsrc;
  i = q % Nsrc;
  if(!i)
    i = Nsrc;
  i --;
  j = (q - i) / Nsrc;
}

void calc_violations(double * G_small, double * W, double * R, double * violations, int Ch, int Sr, int T)
{
    using namespace std;
    int Nsrc = Sr * Sr * 2;
    int Nsen = Ch * Ch * 2; /* Number of elements in a column */
    int s, t, k, l, m; // iterators
    int i, j; // indices
    int half, left = 0, right = 1;
    double norm = 0., norm_ = 0.;
    double result[T];
    double prod;
    double * column = new double[Nsen];
    double * temp = new double[Ch];
    double * colRmat = new double[Ch*Ch];
    double * Gj = new double[Ch];
    double * Gi = new double[Ch];
    cout << "Entering parallel section..." << endl;
    // #pragma omp parallel num_threads(8)
    // {
      // cout << "nthreads = " << omp_get_max_threads();
      // #pragma omp parallel for num_threads(8)
    for (s = 0; s < Nsrc; ++s)
    {
      columnG_fast(s+1, G_small, W, Ch, Sr, column);
      cblas_dgemv(CblasRowMajor, CblasTrans, Nsen, T, 1, R, T, column, 1, 0., result, 1);
      norm = cblas_dnrm2(T, result, 1);
      violations[s] = norm;
      // if(!(s % 10000))
        cout << "s = "<< s << endl;


      // ---------------------- //
      Get_i_j_from_s(s+1, Sr, i, j, half);
      cout << "i = " << i << " j = " << j << " half = " << half << endl;
      for (m = 0; m < Ch; ++m)
      {
        Gj[m] = G_small[Sr * m + j];
        Gi[m] = G_small[Sr * m + i];
        cout << "Gj: " << Gj[m] << " Gi: " << Gi[m] << endl; 
      }
      for (int t = 0; t < T; ++t)
      {
        if(half == left)
           for (k = 0; k < Ch; ++k)
             for (l = 0; l < Ch; ++l)
               colRmat[Ch * k + l] = R[T * (Ch * k + l) + t];
        else if(half == right)
           for (k = 0; k < Ch; ++k)
             for (l = 0; l < Ch; ++l)
               colRmat[Ch * k + l] = R[T * (Ch * k + l + Ch * Ch) + t];
        else 
          cout << "ERROR! calc_violations: half was calculated with error\n";
        cblas_dgemv(CblasRowMajor, CblasNoTrans, Ch, Ch, 1, colRmat, Ch, Gi, 1, 0., temp, 1);
        prod = cblas_ddot(Ch, Gj, 1, temp, 1);
        prod *= W[s];
        norm_ += prod*prod;
      }
      // norm_ = sqrt(norm_);
      cout << "norm = " << norm * norm << " norm_ = " << norm_ << endl;
      norm_ = 0;
      // ---------------------- //
    }
    // }
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

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
    cout << "Entering parallel section..." << endl;
    #pragma omp parallel num_threads(8)
    {
      int Nsrc = Sr * Sr * 2;
      int Nsen = Ch * Ch * 2; /* Number of elements in a column */
      int s, t, k, l, m; // iterators
      int i, j; // indices
      int half, left = 0, right = 1;
      double norm_sq = 0.;
      double result[T];
      double prod;
      double * temp = new double[Ch];
      double ** colRmat_l = new double * [T];
      for (t = 0; t < T; ++t)   
      {
        colRmat_l[t] = new double[Ch*Ch];
        for (k = 0; k < Ch; ++k)
          for (l = 0; l < Ch; ++l)
            colRmat_l[t][Ch * k + l] = R[T * (Ch * k + l) + t];
      }   
      double ** colRmat_r = new double * [T];
      for (t = 0; t < T; ++t)   
      {
        colRmat_r[t] = new double[Ch*Ch];
        for (k = 0; k < Ch; ++k)
          for (l = 0; l < Ch; ++l)
            colRmat_r[t][Ch * k + l] = R[T * (Ch * k + l + Ch * Ch) + t];
      }
      double ** colRmat;
      double ** G = new double * [Sr];
      for (s = 0; s < Sr; ++s)
      { 
        G[s] = new double[Ch];
        for (m = 0; m < Ch; ++m)
            G[s][m] = G_small[Sr * m + s];
      }
        // cout << "nthreads = " << omp_get_max_threads();
      #pragma omp  for 
      for (s = 0; s < Nsrc; ++s)
      {     
        if(!(s % 800000))
          {
            printf(".");
            fflush(stdout);
          }
        // ---------------------- //
        Get_i_j_from_s(s+1, Sr, i, j, half);
        // cout << "i = " << i << " j = " << j << " half = " << half << endl;
        if(half == left)
          colRmat = colRmat_l;
        else if(half == right)
          colRmat = colRmat_r;
        else
          cout << "ERROR! calc_violations: half was calculated with error\n";
        
        for (int t = 0; t < T; ++t)
        {
          cblas_dgemv(CblasRowMajor, CblasNoTrans, Ch, Ch, W[s], colRmat[t], Ch, G[i], 1, 0., temp, 1);
          prod = cblas_ddot(Ch, G[j], 1, temp, 1);
          norm_sq += prod*prod;
        }
        violations[s] = sqrt(norm_sq);
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

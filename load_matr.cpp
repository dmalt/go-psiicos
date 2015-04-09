#include <istream>
#include <string>
#include <sstream>
#include <vector>
#include <omp.h>
#include <Eigen/Dense>

void calc_violations(const Eigen::MatrixXd & G_small, const Eigen::MatrixXd & W, const Eigen::MatrixXd & R, Eigen::MatrixXd & violations, int Nsen, int Nsrc, int T);

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
    for (k = 0; k < Nsen; ++k)
        for (l = 0; l < Nsen; ++l)
        {
            if(half == left)
            {
                G_col[k + Nsen * l] = G_small[k + Nsen * i] * G_small[l + Nsen * j] * w_p;
                G_col[Nsen * Nsen + k + Nsen * l] = 0.;
            }
            else if(half == right)
            {
                G_col[k + Nsen * l] = 0.;
                G_col[Nsen * Nsen + k + Nsen * l] = G_small[k + Nsen * i] * G_small[l + Nsen * j] * w_p;
            }
        }
}

void calc_violations(const Eigen::MatrixXd & G_small, const Eigen::MatrixXd & W, const Eigen::MatrixXd & R, Eigen::MatrixXd & violations, int Nsen, int Nsrc, int T)
{
    int Sr = Nsrc * Nsrc * 2;
    int Sn = Nsen * Nsen * 2; /* Number of elements in a column */
    int s, i, nthreads;
    double norm = 0.;
    /*mexPrintf("T = %d\n", T);*/
/*  #pragma omp parallel firstprivate(column, result, norm) shared(s, Sr, Sn, T, R, lambda, violations) private(nthreads) num_threads(8)
    {*/ 
        /*nthreads = omp_get_num_threads();
        printf("Number of threads = %d\n", nthreads);*/
        #pragma omp parallel for private(s) num_threads(8)/*ordered */
        for (s = 0; s < 10; ++s)
        {
            
            // cblas_dgemv(101, 111, T, Sn, 1, R, Sn, column, 1, 0., result, 1);
            /*for(i = 0; i < T; i++)
                mexPrintf("%f,", result[i]);*/
            /*mexPrintf("\n");*/
            // norm = cblas_dnrm2(T, result, 1);
            violations(s) = norm;
            printf("s = %d\n", s);
        }
    /*}*//* mxDestroyArray(column);
*/}

// example
#include <fstream>
#include <iostream>

int main()
{
    using namespace std;
    using namespace Eigen;
    // initialize params

    // read the file
    std::ifstream Gs("G_small.txt");
    std::ifstream Ws("w.txt");
    std::ifstream Rs("R.txt");

    // load matrices and sizes
    std::vector< std::vector<double> > G_v;
    std::vector< std::vector<double> > W_v;
    std::vector< std::vector<double> > R_v;
    load_matrix(&Gs, &G_v);
    load_matrix(&Rs, &R_v);
    load_matrix(&Ws, &W_v);
    int Src = G_v[0].size();
    int Ch = G_v.size();
    long int Nsrc = 2 * Src * Src;
    int Nch = 2 * Ch* Ch;
    int T = R_v[0].size();
    // ---------------------- //

    // Iitialize Eigen matrices //
    MatrixXd G(Ch, Src);
    for (int i = 0; i < Ch; i++)
        G.row(i) = VectorXd::Map(&G_v[i][0], Src);
    MatrixXd R(Nch, T);
    for (int i = 0; i < Nch; i++)
        R.row(i) = VectorXd::Map(&R_v[i][0], T);
    VectorXd W(Nsrc);
    for (int i = 0; i < Nsrc; i++)
        W(i) = W_v[i][0];
    VectorXd V(Nsrc);
    // ------------------------ //

    cout << "G Nraws = " << Ch << endl;
    cout << "G Ncolumns = " << Src << endl;
    cout << "T = " << T << endl;
    return 0;
}

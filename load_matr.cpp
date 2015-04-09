#include <istream>
#include <string>
#include <sstream>
#include <vector>
#include <Eigen/Dense>

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
    // ------------------------ //

    // cout << G << endl;
    cout << W << endl;
    // print out the matrix

    // cout << "The matrix is:" << endl;
    // for (std::vector< std::vector<double> >::const_iterator it = G_v.begin(); it != G_v.end(); ++ it)
    // {
    //     for (std::vector<double>::const_iterator itit = it->begin(); itit != it->end(); ++ itit)
    //         cout << *itit << '\t';

    //     cout << endl;
    // }
    cout << "W_v.size = " << W_v.size() << endl;
    cout << "G Nraws = " << Ch << endl;
    cout << "G Ncolumns = " << Src << endl;
    cout << "T = " << T << endl;
    return 0;
}

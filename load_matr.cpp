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
    std::ifstream G_smalls("G_small.txt");
    std::ifstream ws("w.txt");
    std::ifstream Rs("R.txt");

    // load the matrix
    std::vector< std::vector<double> > G_small;
    std::vector< std::vector<double> > w;
    std::vector< std::vector<double> > R;
    load_matrix(&G_smalls, &G_small);
    load_matrix(&Rs, &R);
    load_matrix(&ws, &w);
    int Src = G_small[0].size(), Ch = G_small.size();
    long int Nsrc = 2 * Src * Src;
    int Nch = 2 * Ch* Ch;
    int T = R[0].size();
    // print out the matrix

    // cout << "The matrix is:" << endl;
    // for (std::vector< std::vector<double> >::const_iterator it = G_small.begin(); it != G_small.end(); ++ it)
    // {
    //     for (std::vector<double>::const_iterator itit = it->begin(); itit != it->end(); ++ itit)
    //         cout << *itit << '\t';

    //     cout << endl;
    // }
    cout << "w.size = " << w.size() << endl;
    cout << "G_small.size x = " << Ch << endl;
    cout << "G_small.size y = " << Src << endl;
    cout << "T = " << T << endl;
    return 0;
}

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>

using namespace std;

int main(int argc, char *argv[])
{
    cout << "Hello World!" << endl;

    Eigen::VectorXd v(3);
    Eigen::VectorXd w(3);


    v[2] = 40;
    w[2] = 0.5;
    std::cout << v << std::endl;
    std::cout << v.norm();
    return 0;
}

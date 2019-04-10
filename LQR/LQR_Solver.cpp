#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <math.h>

using namespace Eigen;

MatrixXd Kronecker(MatrixXd A, MatrixXd B)
{
    MatrixXd C(A.rows() * B.rows(), A.cols() * B.cols());
    for (int a = 0; a < A.rows(); a++)
        for (int b = 0; b < A.cols(); b++)
            for (int c = 0; c < B.rows(); c++)
                for (int d = 0; d < B.cols(); d++)
                    C(a*B.rows()+c, b*B.cols()+d) = A(a,b) * B(c,d);
    return C;
}

MatrixXd PlaceEigenvalues(MatrixXd A, VectorXd B)
{
    int rows = A.rows();
    int cols = A.cols();

    MatrixXd C(rows, cols);
    MatrixXd Ap = MatrixXd::Identity(rows,cols);
    //Controllability Matrix berechnen
    for (int b = 0; b < cols; b++)
    {
        for (int p = 0; p < b; p++) Ap *= A;
        for (int i = 0; i < rows; i++)
        {
            C(i,b) = (Ap * B)(i);
        }
    }
    MatrixXd Ci(rows, cols);
    //Matrix inversieren
    Ci = C.inverse();

    //Vorfaktor der Gleichung erzeugen
    RowVectorXd pre(cols);
    for (int i = 0; i < cols-1; i++)
    {
       pre(i) = 0;
    }
    pre(cols-1) = 1;

    MatrixXd I = MatrixXd::Identity(rows, cols);
    //gewünschte charakteristische Gleichung erzeugen
    MatrixXd Delta(rows, cols);
    Delta = (A + I);
    for (int i = 0; i < rows-1; i++)
    {
        Delta = Delta * (A + 1 * I);
    }
    //Ackermann Formel anwenden
    VectorXd kT(rows);
    kT = pre * Ci * Delta;
    RowVectorXd k(cols);
    k = kT.transpose();
    return k;
}

MatrixXd LyapunovSolver(MatrixXd A, MatrixXd Q)
{
    MatrixXd I = MatrixXd::Identity(A.rows(), A.cols());
    VectorXd Qv(Q.rows()*Q.cols());
    for (int i = 0; i < Q.rows(); i++)
        for (int b = 0; b < Q.cols(); b++)
            Qv(i*Q.cols() + b) = Q(i,b);

    MatrixXd IA = Kronecker(I,A.transpose());
    MatrixXd AI = Kronecker(A.transpose(),I);

    MatrixXd pre = (IA+AI).inverse();
    MatrixXd Xv = pre * -Qv;

    Map<MatrixXd>X(Xv.data(), A.rows(), A.cols());
    return X;
}

MatrixXd IterateRiccati(MatrixXd A ,MatrixXd B ,MatrixXd Q ,MatrixXd R ,MatrixXd K, double stop_condition)
{
    std::cout << ":Start Riccati Assignments:" << std::endl;
    MatrixXd Ak = A;
    MatrixXd Qk = Q;
    MatrixXd Pk(A.rows(),A.cols());
    MatrixXd Pk0(A.rows(),A.cols());
    MatrixXd Kk = K;
    MatrixXd Pks = 1000 * MatrixXd::Identity(A.rows(), A.cols());
    double d;
    int i = 0;
    do
    {
    std::cout << ":Start Riccati Iteration:" << std::endl;
    Pk0 = Pk;
    Ak = A - B * Kk;
    Qk = Q + Kk.transpose() * R * Kk;
    std::cout << ":Start Lyapunov:" << std::endl;
    Pk = LyapunovSolver(Ak,Qk);
    std::cout << ":Lyapunov solved:" << std::endl;
//  std::cout << R.inverse() << std::endl;
//  std::cout << Pk << std::endl;
    Kk = R.inverse() * B.transpose() * Pk;
    std::cout << ":K solved:" << std::endl;
    if (Pks.sum() > Pk.sum()) Pks = Pk;
    if(Pk0 != Pk) d = (Pk-Pk0).sum()/(Pk0).sum();
    std::cout << "Derivative of cost function : " << sqrt(std::abs(d)) << std::endl;
    i++;
    }while ( std::abs(d) >= stop_condition && i <= 50);
    Kk = R.inverse() * B.transpose() * Pks;
    std::cout << Kk << std::endl;
    return Kk;
}

int main()
{
    double Rm = 2;
    double L = 0.5;
    double Kphi = 0.02;
    double J = 0.02;
    double b = 0.01;
    Matrix<double, 3,3, RowMajor> A;
    VectorXd B(3);
    VectorXd Br(3);
    Matrix<double, 3,3, RowMajor> Q;
    Matrix<double, 1,1, RowMajor> R;
    RowVectorXd C(3);

    /*
    A << 0, 1, 0.1, 10;
    B << 0, 0.1;
    Q << 1, 0, 0, 1;
    R << 1;

    A << 0,1,9.81,0;
    B << 0, -1;
    Q << 10000,0,0,0;
    R << 1;
*//*
    A << 0, 1, 0, 0, 0, Kphi/J,  0, -Kphi/L, -Rm/L;
    B << 0, 0, 1/L;
    Q << 0, 0, 0,
         0, 100, 0,
         0, 0, 10;
    R << 1;
*/
    A << -0, Kphi/J, 0,
        -Kphi/L, -Rm/L, 0,
         1, 0, 0;
    B << 0, 1/L, 0;
    Br << 0, 0, -1;
    C << 1,0,0;
    Q << 1,0,0,
         0,1,0,
         0,0,1;
    R << 1;


    MatrixXd k = PlaceEigenvalues(A,B);
    std::cout << k << std::endl;
    RowVectorXd nk = IterateRiccati(A,B,Q,R,k,0.0001);
    //std::cout << "K : " << K << std::endl;
    //std::cout << "ik : " << ik << std::endl;
    std::cout << 1 << std::endl;
    //MatrixXd stable = A + B * nk;
    //EigenSolver<MatrixXd>eig(stable, false);
    std::cout << 2 << std::endl;
    //std::cout << nk << std::endl;
    //std::cout << eig.eigenvalues() << std::endl;

    VectorXd x(A.rows());
    VectorXd xt(A.rows());
    double xi = 0;
    VectorXd xp(A.rows());
    std::cout << 3 << std::endl;
     x << 0,0,0;
    xt << 1500,0,0;
    std::cout << 4 << std::endl;
    for(int i = 0; i < 100; i++)
    {
    std::cout << x(0) << "      " << x(1) << "      " << x(2) << "        " << nk*x << "    " << i << std::endl;
    xp = A * x + B * (nk * (x-xt)) + Br * xt(0);
    x = x + xp * 0.01;
    }
    int a;
    std::cin >> a;
    return 1;
}





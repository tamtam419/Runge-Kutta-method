#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
using namespace std;

const double pi = 3.14159246;
const double G = 6.67384e-11; //[m^3 kg^-1 s^-2]
const double c = 3.0e8; // [m/s]
const double Ps = 100; // [kg m^-1 s^-2]
int column = 0;

double dm_dr (double r, double m, double rho, double P){
    return 4*pi*rho*r*r;
}

double dp_dr (double r, double m, double rho, double P){
    return - (G*(m + 4 * pi * r * r * r * P / (c * c))*(rho + P / (c * c)) )/ (r * r * (1 - 2 * G * m / ( r * c * c)));
}

double find_rho1_rho2 (double P, vector<double> p, vector<double> rho) {
    double p1 = p[0], p2 = p[1], rho1 = rho[0], rho2 = rho[1];

    for (int i = 0; i < column; i++ ) {
        // cout << p[i] << " " << rho[i] << endl;
        if (p[i] <= P && P< p[i+1]){
            p1 = p[i];
            p2 = p[i+1];
            rho1 = rho[i];
            rho2 = rho[i+1];
            break;
        }
    }
    // // printf("p1 = {}, p2 = {}, rho1 = {}, rho2 = {}", p1, p2, rho1, rho2);
    // cout << "P = "<< P << "p1 = " << p1 << ", p2 = " << p2 << " , rho1 = " << rho1 << " , rho2 = " << rho2 << " " << endl;
    return ((rho2 - rho1) * (P - p1) + (p2 - p1) * rho1) / (p2 - p1);
}

int main(){

    vector<double> rho(110,0), p(110,0);
    // double rho[110], p[110];

    ifstream fin("./BPS.dat");
    ofstream fout("q14_output.dat");

    while (fin >> rho[column] >> p[column]){
        rho[column] = rho[column]*1e3; // kg / m^3
        p[column] = p[column]*1e-1; // kg m^-1 s^-2
        column++;
    }

    // 刻み幅
    double dr = 1.0;

    for (int i = 94; i <= column; i++){

        //　初期値
        double r = dr;
        double P = p[i];
        double RHO = rho[i];
        double m = 4.0 * pi * r * r * r * RHO / 3.0;
        // assert(i == 0 && K >=0.12 && K<=0.13);

        double k1[2],k2[2],k3[2],k4[2];

        // fout << "R  P  RHO  m  M_b  m/M_b dp_dr dm_dr 2*G*m/rc^2" << endl;
        while (P >= Ps){ // Ps = 100


            if (abs(dp_dr(r,m,RHO,P))*1.0e2 < P) {
                dr = 1.0;
            }
            else if (P > abs(dp_dr(r,m,RHO,P))*1.0e1 && abs(dp_dr(r,m,RHO,P))*1.0e2 >= P) {
                dr = 1.0;
            }
            else if (P > abs(dp_dr(r,m,RHO,P))*5.0 && abs(dp_dr(r,m,RHO,P))*1.0e1 >= P) {
                dr = 1.0e-1;
            }    
            else if (P > abs(dp_dr(r,m,RHO,P)) && abs(dp_dr(r,m,RHO,P))*5.0 >= P) {
                dr = 1.0e-2;
            }
            else if (abs(dp_dr(r,m,RHO,P)) >= P) {
                dr = 1.0e-3;
            }

            k1[0] = dr * dm_dr(r,m,RHO,P);
            k1[1] = dr * dp_dr(r,m,RHO,P);

            k2[0] = dr * dm_dr(r + dr/2.0, m + k1[0]/2.0, RHO, P + k1[1]/2.0);
            k2[1] = dr * dp_dr(r + dr/2.0, m + k1[0]/2.0, RHO, P + k1[1]/2.0);

            k3[0] = dr * dm_dr(r + dr/2.0, m + k2[0]/2.0, RHO, P + k2[1]/2.0);
            k3[1] = dr * dp_dr(r + dr/2.0, m + k2[0]/2.0, RHO, P + k2[1]/2.0);

            k4[0] = dr * dm_dr(r + dr, m + k3[0], RHO, P + k3[1]);
            k4[1] = dr * dp_dr(r + dr, m + k3[0], RHO, P + k3[1]);
            
            m += (k1[0] + 2.0*k2[0] + 2.0*k3[0] + k4[0])/6.0;
            P += (k1[1] + 2.0*k2[1] + 2.0*k3[1] + k4[1])/6.0;

            double M_b = (c * c * r) / (2.0 * G);

            RHO = find_rho1_rho2(P,p,rho);
            r += dr;
            // break;
            fout << m << " " << 4 * pi * r * r * r * P / (c * c) << " " <<  m + 4 * pi * r * r * r * P / (c * c) << " " << RHO + P / (c * c) << " " << 1 - 2 * G * m / ( r * c * c) << endl;
            // printf("%.20ef\n", m + 4 * pi * r * r * r * P / (c * c));
            // fout << r << " " << P << " " << RHO << " " << m << " " << M_b << " " << m / M_b << " " << dp_dr(r,m,RHO,P) << " " << dm_dr(r,m,RHO,P)<< " " << 2*G*m/(r*c*c) << endl;
            // cout << r << " " << P << " " << RHO << " " << m << " " << M_b << " " << m / M_b << " " << dp_dr(r,m,RHO,P) << " " << dm_dr(r,m,RHO,P)<<  endl;
        }
        // fout << rho[i] << " " << m << " " << r << endl;
        break;
        // fout << r << " " << m << endl;
        // fout << endl;

    }
}
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>

double alpha(int g_i, int g_a, double p, double T, double I){
    double k_B = 8.6173e-5;
    double n_0 = 7.34e21*p/T;
    double tmp = 21.1*g_i/g_a*pow(T/1000.,2.5)/p*exp(-I/T/k_B);
    double a = pow(tmp/(1.+tmp),0.5);
    double diff = 1.;
    while (abs(diff)>0.0001) {
        double n_i = n_0*a/(1.+a);
        double n_a = n_0*(1.-a)/(1.+a);
        double av_sq_z = 2.*n_i/(2.*n_i+n_a);
        double Gamma = 3.42e-4*pow(n_i,0.5)/pow(T,1.5);
        tmp = 21.1*g_i/g_a*pow(T/1000.,2.5)/p*(1.-Gamma/6.*av_sq_z)*exp(-I/k_B/T+Gamma);
        double new_a = pow(tmp/(1.+tmp),0.5);
        diff = new_a-a;
        n_0 = 7.34e21*p/T/(1.-Gamma/6.*av_sq_z);
        a = new_a;
    }
    return a;
}

double Gamma_D(int g_i, int g_a, double p, double T, double I){
    double k_B = 8.6173e-5;
    double n_0 = 7.34e21*p/T;
    double tmp = 21.1*g_i/g_a*pow(T/1000.,2.5)/p*exp(-I/T/k_B);
    double a = pow(tmp/(1.+tmp),0.5);
    double diff = 1.;
    double Gamma;
    while (abs(diff)>0.0001) {
        double n_i = n_0*a/(1.+a);
        double n_a = n_0*(1.-a)/(1.+a);
        double av_sq_z = 2.*n_i/(2.*n_i+n_a);
        Gamma = 3.42e-4*pow(n_i,0.5)/pow(T,1.5);
        tmp = 21.1*g_i/g_a*pow(T/1000.,2.5)/p*(1.-Gamma/6.*av_sq_z)*exp(-I/k_B/T+Gamma);
        double new_a = pow(tmp/(1.+tmp),0.5);
        diff = new_a-a;
        n_0 = 7.34e21*p/T/(1.-Gamma/6.*av_sq_z);
        a = new_a;
    }
    return Gamma;
}

int main() {
    int g_ion = 2;
    int g_atom = 1;
    double ion_pot = 24.5876;
    double max_T = 5000;
    int N = int(max_T/100.);
    std::vector<double> T(N);
    std::vector<double> Gamma(N);
    for (int i = 0; i<N; i++) {
        T[i]=max_T/N*(i+1);
        Gamma[i]=Gamma_D(g_ion,g_atom,1.,T[i],ion_pot);
    }
    std::ofstream fout;
	fout.open("Plotting.txt");
	for (int i = 0; i < N; i++) {
		fout << T[i] << "," << Gamma[i] <<std::endl;
	}
	fout.close();
}
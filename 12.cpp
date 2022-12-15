#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>

// Функция, вычисляющая степень диссоциации
double alpha(int g_i, int g_a, double p, double T, double I){
    double k_B = 8.6173e-5; // постоянная Больцмана [эВ / К]
    double n_0 = 7.34e21 * p / T; // общая концентрация частиц
    double tmp = 21.1 * g_i / g_a * pow(T / 1000., 2.5) / p * exp(-I / T / k_B); // используем формулу Саха (начальное приближение - идеальная плазма)
    double a = pow(tmp / (1. + tmp), 0.5);
    double diff = 1.;
    while (abs(diff) > 0.0001) {
        double n_i = n_0 * a / (1. + a); // концентрация ионов
        double n_a = n_0 * (1. - a) / (1. + a); // концентрация атомов
        double av_sq_z = 2. * n_i / (2. * n_i + n_a); // средний квадрат заряда
        double Gamma = 3.42e-4 * pow(n_i, 0.5) / pow(T, 1.5); // параметр неидеальности Дебая
        tmp = 21.1 * g_i / g_a * pow(T / 1000., 2.5) / p * (1. - Gamma / 6. * av_sq_z) * exp(-I / k_B / T + Gamma); // учёт неидеальности
        double new_a = pow(tmp / (1. + tmp), 0.5);
        diff = new_a - a; // ошибка в определении степени диссоциации
        n_0 = 7.34e21 * p / T / (1. - Gamma / 6. * av_sq_z);
        a = new_a;
    }
    return a;
}

// Функция, вычисляющая электропроводность
double Sigma(double x) {
    double k_B = 1.381e-16; // постоянная Больцмана [эрг / К]
    double e = 4.803e-10; // заряд электрона [ед. СГСЭ]
    double m = 9.109e-34; // масса электрона [г]
    double p_sum = 1.; // давление смеси [атм]
    double T = 2300.; // температура [К]
    double N_sum = 7.34e21 * p_sum / T; // общая концентрация частиц [см^-3]
    double N_Ar = N_sum / (1. + x); // концентрации и парциальные давления атомов аргона и калия
    double N_K = N_sum - N_Ar;
    double p_Ar = p_sum / (1. + x);
    double p_K = p_sum - p_Ar;
    double alpha_K = alpha(1., 2., p_K, 2300., 4.3407); // степень диссоциации калия
    double n_e = N_K * alpha_K / (1. + alpha_K); // концентрация электронов (их рождают атомы калия)
    double sigma_Ar = 1.e-18; // транспортные сечения столкновения с атомами калия и аргона [см^2]
    double sigma_K = 1.e-15;
    double v = sqrt(3*k_B*T/m); // тепловая скорость электронов, [см/с]
    double sigma = n_e * pow(e,2.) / (m * v * N_Ar * (sigma_Ar + x * sigma_K));
    return sigma;
}

int main() {
    int N = 2000;
    double x_max = 0.;
    double s_max = 0.;
    std::vector<double> vec_x (N);
    std::vector<double> vec_s (N);
    for (int e = 0; e < N; e++) {
        vec_x[e] = 0.00001*e;
        vec_s[e] = Sigma(vec_x[e]);
        if (vec_s[e] > s_max){
            s_max = vec_s[e];
            x_max = vec_x[e];
        }
    }
    std::ofstream fout;
	fout.open("Plotting.txt");
	for (int e = 0; e < N; e++) {
		fout << vec_x[e] << "," << vec_s[e] <<std::endl;
	}
	fout.close();
    std::cout<<"x_max = "<<x_max<<std::endl;
    std::cout<<"s_max = "<<s_max<<std::endl;
}
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>

// Функция для вычисления безразмерной энергии уровня по безразмерной координате точки поворота
double epsilon(double bound){
    return pow(cosh(bound), 2.);
}

// Функция для оценки номера энергетического уровня по безразмерной координате точки поворота
double approx_level_number(double coord, double coef, int M) {
    std::vector<double> vec_x(M);
    double h = 2. * coord / M; // шаг расчётной сетки
    for (int i = 0; i < M; i++){
        vec_x[i] = -coord + i * h; // генерация сетки
    }
    double sum = 0.;
    for (int i = 0; i < M; i++){
        sum += pow((epsilon(coord)-epsilon(vec_x[i])), 0.5) * h; // интегрирование
    }
    return sum * coef - 0.5; // вычисляем приближённое значение n, используя правило квантования Бора-Зоммерфельда
}

// Вычисление точного значения координаты точки поворота по номеру уровня
double get_bound(int num, double coef, int M){
    double par = 1.0; // начальное приближение для любого num
    while (approx_level_number(par, coef, M) < num){
        par *= 2; // увеличиваем отрезок интегрирования, пока интеграл меньше нужного
    }
    while (abs(approx_level_number(par, coef, M) - num) > 0.001){
        par -= par / 2. * (approx_level_number(par, coef, M) - num) / abs(approx_level_number(par, coef, M) - num); // метод дихотомии
    }
    return par;
}

int main() {
    setlocale(LC_ALL, "Russian");
    int M = 1000; // число узлов сетки
    double const h = 1.05457e-34; // приведённая постоянная Планка
    double const PI = 3.14159265;
    double mass = 1.6749e-27; // масса (нейтрона)
    double U0 = 8.01e-20; // параметры потенциальной ямы
    double a = 1.e-9;
    double coef = a*pow(2*mass*U0,0.5)/h/PI; // безразмерный коэффициент
    std::vector<double> x_coord(51); // рассчитаем координаты (безразмерные) точек поворота для первых 50 уровней
    std::ofstream fout; // для записи в файл
	fout.open("Graph.txt");
    for (int i = 0; i < 51; i++) {
        x_coord[i] = get_bound(i, coef, M);
        fout << i <<"," <<x_coord[i] << std::endl;
    }
    fout.close();
    return 0;
}
#include <iostream>
#include <math.h>
#include <map>
#include <fstream>

const double R = 8.3144626; // ������������� ������� ����������
const double k_B = 1.380649e-23; // ���������� ���������
const double aem = 1.660539069e-27; // ������� ������� �����

// ������� sgn(x)
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

// ������������� ���������� ���������������� ����������
double Omega11(double tau){
    return 1.2067*pow(tau, -0.1946);
}

double Omega12(double tau){
    return 1.0632*pow(tau, -0.1736);
}

double Omega13(double tau){
    return 0.9908*pow(tau, -0.1643);
}

double Omega22(double tau){
    return 1.3104*pow(tau, -0.1848); 
}

double A_int(double tau){
    return Omega22(tau)/Omega11(tau);
}

double B_int(double tau){
    return (5.*Omega12(tau) - 4.*Omega13(tau))/Omega11(tau);
}

// ������� ������������� ��� ���������� �������� � ����������������
double viscosity_matrix_A[5][5] = {{3.06290, 8.81837, 6.57562, 11.2472, 3.40998},
                                   {8.81837, 8.45280, 10.4450, 15.2944, 18.2681},
                                   {6.57562, 10.4450, 6.98910, 18.5326, 16.9684},
                                   {11.2472, 15.2944, 18.5326, 6.96290, 11.2303},
                                   {3.40998, 18.2681, 16.9684, 11.2303, 7.56830}};

double viscosity_matrix_T[5][5] = {{-21.33, 63.630, 51.870, 121.27, 45.890},
                                   {63.630, 16.470, 47.660, 94.860, 125.00},
                                   {51.870, 47.660, 65.700, 148.93, 151.58},
                                   {121.27, 94.860, 148.93, 71.070, 116.71},
                                   {45.890, 125.00, 151.58, 116.71, 112.31}};

double viscosity_matrix_exp[5][5] = {{0.724300, 0.614098, 0.602712, 0.508158, 0.658754},
                                     {0.614098, 0.642584, 0.620956, 0.565158, 0.522034},
                                     {0.602712, 0.620956, 0.639770, 0.541593, 0.542416},
                                     {0.508158, 0.565158, 0.541593, 0.667000, 0.631571},
                                     {0.658754, 0.522034, 0.542416, 0.631571, 0.655473}};

double conductivity_correction[5][5] = {{1.000, 0.895, 0.940, 1.020, 1.060},
                                        {0.895, 1.000, 0.845, 0.900, 0.940},
                                        {0.940, 0.845, 1.000, 0.850, 0.935},
                                        {1.020, 0.900, 0.850, 1.000, 0.871},
                                        {1.060, 0.940, 0.935, 0.871, 1.000}};

// �����, ����������� �������� ���
class Gas{
public:
    Gas() = default;
    Gas(double M, double T_cr, double p_cr, double rho_cr, double T_boil, double sigma, double eps, bool is_He, int period);
    double coef_2(double T);
    double coef_3(double T);
    double viscosity_CE(double T);
    double conductivity_CE(double T);
    double c_p_gas(double T, double p);
    double c_v_gas(double T, double p);
    double rho_M(double T, double p);
    double M, sigma, eps, m0, V_cr, T_cr;
    int period;
    bool is_He;
private:
    double p_cr, rho_cr, T_boil;
};

Gas::Gas(double M, double T_cr, double p_cr, double rho_cr, double T_boil, double sigma, double eps, bool is_He, int period){
    this->M = M; // �������� ����� 
    this->T_cr = T_cr; // ����������� �����������
    this->p_cr = p_cr; // ����������� ��������
    this->rho_cr = rho_cr; // ����������� ���������
    this->T_boil = T_boil; // ����������� �������
    this->sigma = sigma; // ����������� ������� ������������ (��������� ��������-������)
    this->eps = eps; // �������� ������� ���������� (��������� ��������-������)
    this->m0 = this->M*1.e3*aem; // ����� �����
    this->V_cr = R*this->T_cr/this->p_cr; // �������� ����������� �����
    this->is_He = is_He; // ����� �� ��� (true ��� �����, false ��� ���������)
    this->period = period; // �� ������� ������ ������ ������� � ������� ����������
}

double Gas::coef_2(double T){ // ������ ���������� ����������� ������� ����
    if (this->is_He == true){
        return 1.e-6*(8.4 - 0.0018*T + 115./pow(T,0.5) - 835./T); 
    }
    else{
        double theta = T/this->T_cr; // ���������� �����������
        return (-102.6 + (102.732 - 0.001*theta - 0.44/pow(theta, 1.22))*tanh(4.5*pow(theta, 0.5)))*this->V_cr;
    }
}

double Gas::coef_3(double T){ // ������ ���������� ����������� ������� ����
    if (this->is_He == true){
        return 0.;
    }
    else{
        double theta = T/this->T_cr; // ���������� �����������
        return (0.0757 + (-0.0862 - 3.6e-5*theta + 0.0237/pow(theta, 0.059))*tanh(0.84*theta))*pow(this->V_cr, 2.);
    }
}

double Gas::viscosity_CE(double T){ // ������������ �������� �� ������ ������� -- �������
    // double tau = T/this->eps; // ���������� ����������� ��� ���������� ���������� ���������������� ����������
    // return 5./16.*pow(M_PI*this->m0*k_B*T, 0.5)/(M_PI*pow(this->sigma, 2.)*Omega22(tau));
    return 1.e-7*viscosity_matrix_A[this->period][this->period]*pow((T - viscosity_matrix_T[this->period][this->period]),viscosity_matrix_exp[this->period][this->period]);
}

double Gas::conductivity_CE(double T){ // ���������������� �� ������ ������� -- �������
    // double tau = T/this->eps; // ���������� ����������� ��� ���������� ���������� ���������������� ����������
    // return 25./32.*pow(M_PI*this->m0*k_B*T, 0.5)/(M_PI*pow(this->sigma, 2.)*Omega22(tau))*3./2.*k_B/this->m0;
    return 15./4.*k_B/this->m0*this->viscosity_CE(T);
}

double Gas::c_p_gas(double T, double p){ // �������� �������� ����������� ��� ���������� ��������
    double dB, d2B, dC, d2C; // ����������� ������� � �������� ����������� ������������� ������������ ������������
    double B = this->coef_2(T);
    double C = this->coef_3(T);
    double r = this->rho_M(T, p);
    double theta = T/this->T_cr; // ���������� �����������
    if (this->is_He == true){
        dB = (-0.0018 - 115./2./pow(T, 1.5) + 835./pow(T, 2.))*1.e-6;
        d2B = (115.*3./4./pow(T, 2.5) - 835.*2./pow(T, 3.))*1.e-6;
        dC = 0.;
        d2C = 0.;
    }
    else{
        double a = 102.732; 
        double b = -0.001; 
        double c = -0.44;
        double d = 4.5;
        double alp = 1.22;
        double bet = 0.5;
        dB = this->V_cr*((b-alp*c/pow(theta, alp+1.))*tanh(d*pow(theta, bet))+(a+b*theta+c/pow(theta, alp))*d*bet*pow(theta, bet-1.)/pow(cosh(d*pow(theta, bet)), 2.))/this->T_cr;
        d2B = this->V_cr*(alp*(alp+1.)*c/pow(theta, alp+2.)*tanh(d*pow(theta, bet))+2.*(b-alp*c/pow(theta, alp+1.))*d*bet*pow(theta, bet-1.)/pow(cosh(d*pow(theta, bet)), 2.)+(a+b*theta+c/pow(theta, alp))*d*bet*((bet-1.)*pow(theta, bet-2.)*pow(cosh(d*pow(theta, bet)), 2.)-pow(theta, bet-1.)*2.*sinh(d*pow(theta, bet))*cosh(d*pow(theta, bet))*d*bet*pow(theta, bet-1.))/pow(cosh(d*pow(theta, bet)), 4.))/pow(this->T_cr, 2.);
        a = -0.0862; 
        b = -3.6e-5; 
        c = 0.0237;
        d = 0.84;
        alp = 0.059;
        bet = 1.;
        dC = pow(this->V_cr, 2.)*((b-alp*c/pow(theta, alp+1.))*tanh(d*pow(theta, bet))+(a+b*theta+c/pow(theta, alp))*d*bet*pow(theta, bet-1.)/pow(cosh(d*pow(theta, bet)), 2.))/this->T_cr;
        d2C = pow(this->V_cr, 2.)*(alp*(alp+1.)*c/pow(theta, alp+2.)*tanh(d*pow(theta, bet))+2.*(b-alp*c/pow(theta, alp+1.))*d*bet*pow(theta, bet-1.)/pow(cosh(d*pow(theta, bet)), 2.)+(a+b*theta+c/pow(theta, alp))*d*bet*((bet-1.)*pow(theta, bet-2.)*pow(cosh(d*pow(theta, bet)), 2.)-pow(theta, bet-1.)*2.*sinh(d*pow(theta, bet))*cosh(d*pow(theta, bet))*d*bet*pow(theta, bet-1.))/pow(cosh(d*pow(theta, bet)), 4.))/pow(this->T_cr, 2.);
    }
    return 5./2.*R+r*R*((B-T*dB-pow(T,2.)*d2B)+r*(C-pow(T,2.)/2.*d2C))+R*T*((B-T*dB)+r*(2.*C-T*dC))*(-1.)*((r+B*pow(r,2.)+C*pow(r,3.))/T+dB*pow(r,2.)+dC*pow(r,3.))/(1.+2.*B*r+3.*C*pow(r,2.));
}

double Gas::c_v_gas(double T, double p){ // �������� �������� ����������� ��� ���������� ������
    double dB, d2B, dC, d2C;
    double B = this->coef_2(T);
    double C = this->coef_3(T);
    double theta = T/this->T_cr;
    double r = this->rho_M(T, p);
    if (this->is_He == true){
        dB = (-0.0018-115./2./pow(T, 1.5)+835./pow(T, 2.))*1.e-6;
        d2B = (115.*3./4./pow(T, 2.5)-835.*2./pow(T, 3.))*1.e-6;
        dC = 0.;
        d2C = 0.;
    }
    else{
        double a = 102.732; 
        double b = -0.001; 
        double c = -0.44;
        double d = 4.5;
        double alp = 1.22;
        double bet = 0.5;
        dB = this->V_cr*((b-alp*c/pow(theta, alp+1.))*tanh(d*pow(theta, bet))+(a+b*theta+c/pow(theta, alp))*d*bet*pow(theta, bet-1.)/pow(cosh(d*pow(theta, bet)), 2.))/this->T_cr;
        d2B = this->V_cr*(alp*(alp+1.)*c/pow(theta, alp+2.)*tanh(d*pow(theta, bet))+2.*(b-alp*c/pow(theta, alp+1.))*d*bet*pow(theta, bet-1.)/pow(cosh(d*pow(theta, bet)), 2.)+(a+b*theta+c/pow(theta, alp))*d*bet*((bet-1.)*pow(theta, bet-2.)*pow(cosh(d*pow(theta, bet)), 2.)-pow(theta, bet-1.)*2.*sinh(d*pow(theta, bet))*cosh(d*pow(theta, bet))*d*bet*pow(theta, bet-1.))/pow(cosh(d*pow(theta, bet)), 4.))/pow(this->T_cr, 2.);
        a = -0.0862; 
        b = -3.6e-5; 
        c = 0.0237;
        d = 0.84;
        alp = 0.059;
        bet = 1.;
        dC = pow(this->V_cr, 2.)*((b-alp*c/pow(theta, alp+1.))*tanh(d*pow(theta, bet))+(a+b*theta+c/pow(theta, alp))*d*bet*pow(theta, bet-1.)/pow(cosh(d*pow(theta, bet)), 2.))/this->T_cr;
        d2C = pow(this->V_cr, 2.)*(alp*(alp+1.)*c/pow(theta, alp+2.)*tanh(d*pow(theta, bet))+2.*(b-alp*c/pow(theta, alp+1.))*d*bet*pow(theta, bet-1.)/pow(cosh(d*pow(theta, bet)), 2.)+(a+b*theta+c/pow(theta, alp))*d*bet*((bet-1.)*pow(theta, bet-2.)*pow(cosh(d*pow(theta, bet)), 2.)-pow(theta, bet-1.)*2.*sinh(d*pow(theta, bet))*cosh(d*pow(theta, bet))*d*bet*pow(theta, bet-1.))/pow(cosh(d*pow(theta, bet)), 4.))/pow(this->T_cr, 2.);
    }
    return 3./2.*R-r*R*T*(2.*dB+T*d2B+r*(dC+T/2.*d2C));
}

double Gas::rho_M(double T, double p){ // ���������� �������� ��������� ������� ���� �� ����������� � ��������
    double a_coef = this->coef_3(T)*R*T/p;
    double b_coef = this->coef_2(T)*R*T/p;
    double c_coef = R*T/p;
    double rho0 = c_coef; // ��������� ����������� -- ����������� ���
    double er = 1.;
    while (abs(er)>1.e-7){ // ���������� ��������� ������ �������� ������� �������
        er = -(a_coef*pow(rho0, 3.)+b_coef*pow(rho0, 2.)+c_coef*rho0-1.)/(3.*a_coef*pow(rho0, 2.)+2.*b_coef*rho0+c_coef);
        rho0 += er;
    }
    return rho0;
}

// ������� ������� ������ �������� �����
Gas He(4.003e-3, 5.20, 2.275e5, 69.64, 13.78, 2.576e-10, 10.22, true, 0);
Gas Ne(20.179e-3, 44.50, 26.78e5, 481.9, 27.09, 2.789e-10, 35.7, false, 1);
Gas Ar(39.948e-3, 150.7, 48.63e5, 535.6, 87.29, 3.418e-10, 124.0, false, 2);
Gas Kr(83.800e-3, 209.5, 55.1e5, 908.4, 119.78, 3.61e-10, 190., false, 3);
Gas Xe(131.29e-3, 289.7, 58.4e5, 1110., 160.05, 4.055e-10, 229., false, 4);

// ������� "(������ � ������� ����������) -- ������ ���� Gas"
std::map<int, Gas> map_of_gases{
    {1,He}, {2, Ne}, {3, Ar}, {4, Kr}, {5, Xe}
};

// �����, ����������� �������� ����� �������� �����
class Mixture{
public:
    Mixture() = default;
    double mix_coef_2(double T);
    double mix_coef_3(double T);
    double molar_density(double T, double p);
    double density(double T, double p);
    double mix_add_mu(double T, double p);
    double mix_add_lambda(double T, double p);
    double mix_mu_CE(double T);
    double mix_lambda_CE(double T);
    double mix_R(double T, double p);
    double viscosity(double T, double p);
    double kinematic_viscosity(double T, double p);
    double conductivity(double T, double p);
    double c_p(double T, double p);
    double c_v(double T, double p);
    double kappa(double T, double p);
    double compress(double T, double p);
    double Prandtl (double T, double p);
    double friction (double T, double p);
    double ideal_friction (double T, double p);
    double delta_friction (double T, double p);
    double delta_viscosity(double T,double p);
    double mix_M;

    Gas gas1, gas2;
    int period1, period2;
    double x1, x2, mix_sigma, mix_eps;
};

// Mixture::Mixture(int gas_name_1, int gas_name_2, double x1){
//     num1 = gas_name_1;
//     num2 = gas_name_2;
//     this->gas1 = map_of_gases[gas_name_1]; // ���������� �����
//     this->gas2 = map_of_gases[gas_name_2];
//     this->x1 = x1; // ������� ���� �����������
//     this->x2 = 1.-x1;
//     this->mix_sigma = 0.5*(gas1.sigma+gas2.sigma); // ����������� ������� ������������ (��������� �������� -- ������)
//     this->mix_eps = pow(gas1.eps*gas2.eps, 0.5); // ����������� ������� �������������� ������ (��������� ��������-������)
//     this->mix_M = x1*gas1.M+(1.-x1)*gas2.M; // �������� ����� �����
// };

double Mixture::mix_coef_2(double T){ // ������ ���������� ����������� ��� �����
    double B1 = this->gas1.coef_2(T);
    double B2 = this->gas2.coef_2(T);
    double beta = this->gas1.V_cr/this->gas2.V_cr; // ��������� �������� ����������� �������
    double T12 = 4.*beta/pow(1.+beta, 2.)*pow(this->gas1.T_cr*this->gas2.T_cr, 0.5);
    double V_crit_mix = 0.5*(this->gas1.V_cr+this->gas2.V_cr);
    double theta_mix = T/T12;
    double B12 = (-102.6 + (102.732 - 0.001*theta_mix - 0.44/pow(theta_mix, 1.22))*tanh(4.5*pow(theta_mix, 0.5)))*V_crit_mix;
    return pow(this->x1, 2.)*B1+2.*this->x1*this->x2*B12+pow(this->x2, 2.)*B2;
}

double Mixture::mix_coef_3(double T){ // ������ ���������� ����������� ��� ����� �����
    double C1 = this->gas1.coef_3(T);
    double C2 = this->gas2.coef_3(T);
    double C112 = pow(abs(C1), 2./3.)*pow(abs(C2), 1./3.)*sgn(C2);
    double C221 = pow(abs(C2), 2./3.)*pow(abs(C1), 1./3.)*sgn(C1);
    return pow(this->x1, 3.)*C1+3.*pow(this->x1, 2.)*this->x2*C112+3.*this->x1*pow(this->x2, 2.)*C221+pow(this->x2, 3.)*C2;
}

double Mixture::c_p(double T, double p){ // �������� �������� ����������� ����� ��� ���������� ��������
    double B = this->mix_coef_2(T);
    double C = this->mix_coef_3(T);
    double r = this->molar_density(T, p);
    // ����������� ���������� ������������� �� ����������� ��������� ��������
    double dT = 0.01;
    double dB = (this->mix_coef_2(T+dT)-this->mix_coef_2(T-dT))/(2.*dT);
    double dC = (this->mix_coef_3(T+dT)-this->mix_coef_3(T-dT))/(2.*dT);
    double d2B = (this->mix_coef_2(T+dT)-2.*this->mix_coef_2(T)+this->mix_coef_2(T-dT))/(pow(dT,2.));
    double d2C = (this->mix_coef_3(T+dT)-2.*this->mix_coef_3(T)+this->mix_coef_3(T-dT))/(pow(dT,2.));
    return 5./2.*R+r*R*((B-T*dB-pow(T,2.)*d2B)+r*(C-pow(T,2.)/2.*d2C))+R*T*((B-T*dB)+r*(2.*C-T*dC))*(-1.)*((r+B*pow(r,2.)+C*pow(r,3.))/T+dB*pow(r,2.)+dC*pow(r,3.))/(1.+2.*B*r+3.*C*pow(r,2.));
}

double Mixture::c_v(double T, double p) { // �������� �������� ����������� ����� ��� ���������� ������
    double B = this->mix_coef_2(T);
    double C = this->mix_coef_3(T);
    double r = this->molar_density(T, p);
    // ����������� ���������� ������������� �� ����������� ��������� ��������
    double dT = 0.01;
    double dB = (this->mix_coef_2(T+dT)-this->mix_coef_2(T-dT))/(2.*dT);
    double dC = (this->mix_coef_3(T+dT)-this->mix_coef_3(T-dT))/(2.*dT);
    double d2B = (this->mix_coef_2(T+dT)-2.*this->mix_coef_2(T)+this->mix_coef_2(T-dT))/(pow(dT,2.));
    double d2C = (this->mix_coef_3(T+dT)-2.*this->mix_coef_3(T)+this->mix_coef_3(T-dT))/(pow(dT,2.));
    return 3./2.*R-r*R*T*(2.*dB+T*d2B+r*(dC+T/2.*d2C));
}

double Mixture::molar_density(double T, double p){ // ������ �������� ��������� �����
    double a_coef = this->mix_coef_3(T)*R*T/p;
    double b_coef = this->mix_coef_2(T)*R*T/p;
    double c_coef = R*T/p;
    double rho0 = c_coef; // ��������� ����������� -- ����������� ���
    double er = 1.;
    while (abs(er)>1.e-7){ // ���������� ��������� ������ �������� ������� �������
        er = -(a_coef*pow(rho0, 3.)+b_coef*pow(rho0, 2.)+c_coef*rho0-1.)/(3.*a_coef*pow(rho0, 2.)+2*b_coef*rho0+c_coef);
        rho0 += er;
    }
    return rho0;
}

double Mixture::density(double T, double p){ // ��������� �������� �����
    return this->molar_density(T, p)*this->mix_M;
}

double Mixture::mix_add_mu(double T, double p){ // ������������� ������� ��������� �������� � ��������
    double beta = this->gas1.V_cr/this->gas2.V_cr; // ��������� �������� ����������� ������� �����

    double V12 = 0.5*(this->gas1.V_cr+this->gas2.V_cr); // ������������������ ����� � ����������� ����������� �� ������ ������������� ��� �������
    double T12 = 4.*beta/pow(1.+beta, 2.)*pow(this->gas1.T_cr*this->gas2.T_cr, 0.5);

    double V_tmp = pow(this->x1,2.)*this->gas1.V_cr+pow(this->x2,2.)*this->gas2.V_cr+2.*this->x1*this->x2*V12; // ����������������� ����� � ����������� ������������ �� �������� ��� ��� �������
    double T_tmp = 1./V_tmp*(pow(this->x1,2.)*this->gas1.V_cr*this->gas1.T_cr+pow(this->x2, 2.)*this->gas2.V_cr*this->gas2.T_cr+2.*this->x1*this->x2*V12*T12);
    double a_mu = 0.291*V_tmp*this->molar_density(T,p);
    double psi_mu = 0.221*a_mu+1.062*pow(a_mu, 2.)-0.509*pow(a_mu, 3.)+0.225*pow(a_mu, 4.);
    double mu_star = 0.204e-7*pow(this->mix_M*T_tmp, 0.5)/pow(0.291*V_tmp, 2./3.);
    return (1. - this->x1*gas1.is_He - this->x2*gas2.is_He)*(1.-1./2.3)*mu_star*psi_mu; // changed 0n 06.05.2025
}

double Mixture::mix_add_lambda(double T, double p){ // ������������� ������� ��������� �������� � ����������������
    double beta = this->gas1.V_cr/this->gas2.V_cr; // ��������� �������� ����������� ������� �����

    double V12 = 0.5*(this->gas1.V_cr+this->gas2.V_cr); // ������������������ ����� � ����������� ����������� �� ������ ������������� ��� �������
    double T12 = 4.*beta/pow(1.+beta, 2.)*pow(this->gas1.T_cr*this->gas2.T_cr, 0.5);

    double V_tmp = pow(this->x1,2.)*this->gas1.V_cr+pow(this->x2,2.)*this->gas2.V_cr+2.*this->x1*this->x2*V12; // ����������������� ����� � ����������� ������������ �� �������� ��� ��� �������
    double T_tmp = 1./V_tmp*(pow(this->x1,2.)*this->gas1.V_cr*this->gas1.T_cr+pow(this->x2, 2.)*this->gas2.V_cr*this->gas2.T_cr+2.*this->x1*this->x2*V12*T12);
    double a_lam = 0.291*V_tmp*this->molar_density(T,p);
    double psi_lam = 0.645*a_lam+0.331*pow(a_lam, 2.)+0.0368*pow(a_lam, 3.)-0.0128*pow(a_lam, 4.);
    double lam_star = 0.304e-4*pow(T_tmp, 0.277)/pow(this->mix_M, 0.465)/pow(0.291*V_tmp, 0.415);
    return (1.-1./2.94)*lam_star*psi_lam;
}

double Mixture::mix_mu_CE(double T){ // ������������ �������� ����� � ������ ������� -- �������
    double tau = T/this->mix_eps;
    double m1 = this->gas1.m0;
    double m2 = this->gas2.m0;
    double mu1 = this->gas1.viscosity_CE(T);
    double mu2 = this->gas2.viscosity_CE(T);
    double mu12 = 1.e-7*viscosity_matrix_A[period1][period2]*pow((T-viscosity_matrix_T[period1][period2]),viscosity_matrix_exp[period1][period2]);
    double phi12 = mu1/mu12*(2.*m1*m2)/(pow(m1+m2,2.))*(5./3./A_int(tau)+m2/m1);
    double phi21 = mu2/mu12*(2.*m1*m2)/(pow(m1+m2,2.))*(5./3./A_int(tau)+m1/m2);
    return mu1/(1.+phi12*x2/x1)+mu2/(1.+phi21*x1/x2);
}

double Mixture::mix_lambda_CE(double T){ // ���������������� ����� � ������ ������� -- �������
    double tau = T/this->mix_eps;
    double m1 = this->gas1.m0;
    double m2 = this->gas2.m0;
    double m12 = 2.*m1*m2/(m1 + m2);
    double lam1 = this->gas1.conductivity_CE(T);
    double lam2 = this->gas2.conductivity_CE(T);
    if (x1*x2 < 1.e-6){
        return x1*lam1 + x2*lam2;
    }
    else{
        double mu12 = 1.e-7*viscosity_matrix_A[period1][period2]*pow((T-viscosity_matrix_T[period1][period2]),viscosity_matrix_exp[period1][period2]);
        double lam12 = 15./4.*k_B/m12*mu12*conductivity_correction[period1][period2];
        double L11 = pow(x1,2.)/lam1 + x1*x2/2./lam12*(15./2.*pow(m1,2.)+25./4.*pow(m2,2.)-3.*pow(m2,2.)*B_int(tau)+4.*m1*m2*A_int(tau))/(pow(m1+m2,2.)*A_int(tau));
        double L22 = pow(x2,2.)/lam2 + x1*x2/2./lam12*(15./2.*pow(m2,2.)+25./4.*pow(m1,2.)-3.*pow(m1,2.)*B_int(tau)+4.*m1*m2*A_int(tau))/(pow(m1+m2,2.)*A_int(tau));
        double L12 = -x1*x2/2./lam12*m1*m2/pow(m1+m2,2.)/A_int(tau)*(55./4.-3.*B_int(tau)-4.*A_int(tau));
        return (pow(x1,2.)/L11-2.*x1*x2*L12/L11/L22+pow(x2,2.)/L22)/(1.-pow(L12,2.)/(L11*L22));
    }
}

double Mixture::mix_R(double T, double p){ // ������� ���������� �����
    return p/T/this->density(T,p);
}

double Mixture::viscosity(double T, double p){ // ������������ �������� �����
    return this->mix_mu_CE(T)+this->mix_add_mu(T,p);
}

double Mixture::conductivity(double T, double p){ // ���������������� �����
    return this->mix_lambda_CE(T)+this->mix_add_lambda(T,p);
}

double Mixture::kappa(double T, double p){ // ���������� �������� �����
    return this->c_p(T,p)/this->c_v(T,p);
}

double Mixture::compress(double T, double p) { // ����������� ����������� �����
    return 1.+this->mix_coef_2(T)*this->molar_density(T,p)+this->mix_coef_3(T)*pow(this->molar_density(T,p),2.);
}

double Mixture::kinematic_viscosity(double T, double p) {
    return this->viscosity(T,p)/this->density(T,p);
}

double Mixture::Prandtl(double T, double p) {
    return this->viscosity(T,p)*this->c_p(T,p)/this->conductivity(T,p)/this->mix_M;
}

double Mixture::friction(double T, double p) {
    return pow(this->kinematic_viscosity(T,p),0.2)*this->density(T,p); // changed on 06.05.25
}

double Mixture::ideal_friction(double T, double p) {
    double r = p*mix_M/(R*T);
    return pow((this->mix_mu_CE(T)/r),0.2)*r; // added on 06.05.25
}

double Mixture::delta_friction(double T, double p){
    return (friction(T,p) - ideal_friction(T,p)) / ideal_friction(T,p) * 100.;
}

double Mixture::delta_viscosity(double T, double p){
    return (viscosity(T,p) - mix_mu_CE(T)) / mix_mu_CE(T) * 100.;
}



class properties_by_percentage : public Mixture{
public:
    properties_by_percentage(Gas gas1, Gas gas2, double x){
        this->gas1 = gas1; // ���������� �����
        this->gas2 = gas2;
        this->period1 = gas1.period;
        this->period2 = gas2.period;
        this->x1 = x; // ������� ���� �����������
        this->x2 = 1.-x1;
        this->mix_sigma = 0.5*(gas1.sigma+gas2.sigma); // ����������� ������� ������������ (��������� �������� -- ������)
        this->mix_eps = pow(gas1.eps*gas2.eps, 0.5); // ����������� ������� �������������� ������ (��������� ��������-������)
        this->mix_M = x1*gas1.M + (1.-x1)*gas2.M; // �������� ����� �����
    }
};

class properties_by_molar_mass : public Mixture{
public:
    properties_by_molar_mass(Gas gas1, Gas gas2, double x){
        this->gas1 = gas1; // ���������� �����
        this->gas2 = gas2;
        this->period1 = gas1.period;
        this->period2 = gas2.period;
        this->x1 = (x-gas2.M)/(gas1.M-gas2.M); // ������� ���� �����������
        this->x2 = 1.-x1;
        this->mix_sigma = 0.5*(gas1.sigma+gas2.sigma); // ����������� ������� ������������ (��������� �������� -- ������)
        this->mix_eps = pow(gas1.eps*gas2.eps, 0.5); // ����������� ������� �������������� ������ (��������� ��������-������)
        this->mix_M = x; // �������� ����� �����
    }
};


double data_to_write(int name_gas_1, int name_gas_2, double molar, double T, double p){
    Gas G1 = map_of_gases[name_gas_1];
    Gas G2 = map_of_gases[name_gas_2];
    return properties_by_molar_mass(G1, G2, 0.001*molar).delta_viscosity(T,p);
}



int main() {

    std::ofstream fout;

    fout.open("DeltaViscoXe.txt");

    for (int i = 0; i<100; i++){
        double t = 400. + (1200. - 400.)/99*i;
        fout<<t<<","<<properties_by_percentage(Xe, He, 1.).delta_viscosity(t,1.e6)<<","<<properties_by_percentage(Xe, He, 1.).delta_viscosity(t,2.e6)<<","<<properties_by_percentage(Xe, He, 1.).delta_viscosity(t,3.e6)<<","<<properties_by_percentage(Xe, He, 1.).delta_viscosity(t,4.e6)<<","<<properties_by_percentage(Xe, He, 1.).delta_viscosity(t,5.e6)<<","<<properties_by_percentage(Xe, He, 1.).delta_viscosity(t,6.e6)<<","<<properties_by_percentage(Xe, He, 1.).delta_viscosity(t,7.e6)<<std::endl;
    }

    // int num_for_plot = 100; // ����� ����� ��� ���������� �������

    // std::ofstream fout; // ��� ������ � ����

    // double T_test, p_test;

    // std::string property_path = "Visco_CE"; // ����� �� ��������� ����

    // std::string parameter_path; // ����� � ����������� (����������� � ��������)

    // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // parameter_path = "/400_2/"; 

    // T_test = 400.; // �����������, �
    // p_test = 2.e6; // ��������, ��

   	// fout.open(property_path + parameter_path + "HeNe.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(1,2,4.003+(20.183-4.003)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "HeAr.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(1,3,4.003+(39.948-4.003)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "HeKr.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(1,4,4.003+(83.80-4.003)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "HeXe.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(1,5,4.003+(131.3-4.003)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "NeAr.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(2,3,20.183+(39.948-20.183)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "NeKr.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(2,4,20.183+(83.80-20.183)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "NeXe.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(2,5,20.183+(131.3-20.183)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "ArKr.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(3,4,39.948+(83.80-39.948)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "ArXe.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(3,5,39.948+(131.3-39.948)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "KrXe.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(4,5,83.8+(131.3-83.8)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // parameter_path = "/400_7/"; 

    // T_test = 400.; 
    // p_test = 7.e6; 

   	// fout.open(property_path + parameter_path + "HeNe.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(1,2,4.003+(20.183-4.003)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "HeAr.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(1,3,4.003+(39.948-4.003)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "HeKr.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(1,4,4.003+(83.80-4.003)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "HeXe.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(1,5,4.003+(131.3-4.003)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "NeAr.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(2,3,20.183+(39.948-20.183)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "NeKr.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(2,4,20.183+(83.80-20.183)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "NeXe.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(2,5,20.183+(131.3-20.183)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "ArKr.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(3,4,39.948+(83.80-39.948)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "ArXe.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(3,5,39.948+(131.3-39.948)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "KrXe.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(4,5,83.8+(131.3-83.8)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // parameter_path = "/1200_2/"; 

    // T_test = 1200.; 
    // p_test = 2.e6; 

   	// fout.open(property_path + parameter_path + "HeNe.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(1,2,4.003+(20.183-4.003)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "HeAr.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(1,3,4.003+(39.948-4.003)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "HeKr.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(1,4,4.003+(83.80-4.003)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "HeXe.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(1,5,4.003+(131.3-4.003)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "NeAr.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(2,3,20.183+(39.948-20.183)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "NeKr.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(2,4,20.183+(83.80-20.183)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "NeXe.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(2,5,20.183+(131.3-20.183)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "ArKr.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(3,4,39.948+(83.80-39.948)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "ArXe.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(3,5,39.948+(131.3-39.948)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "KrXe.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(4,5,83.8+(131.3-83.8)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // parameter_path = "/1200_7/"; 

    // T_test = 1200.; 
    // p_test = 7.e6; 

   	// fout.open(property_path + parameter_path + "HeNe.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(1,2,4.003+(20.183-4.003)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "HeAr.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(1,3,4.003+(39.948-4.003)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "HeKr.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(1,4,4.003+(83.80-4.003)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "HeXe.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(1,5,4.003+(131.3-4.003)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "NeAr.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(2,3,20.183+(39.948-20.183)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "NeKr.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(2,4,20.183+(83.80-20.183)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "NeXe.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(2,5,20.183+(131.3-20.183)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "ArKr.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(3,4,39.948+(83.80-39.948)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "ArXe.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(3,5,39.948+(131.3-39.948)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;        
    // }
    // fout.close();

    // fout.open(property_path + parameter_path + "KrXe.txt");
    // for (int i = 0; i<num_for_plot; i++){
    //     fout<<data_to_write(4,5,83.8+(131.3-83.8)/(num_for_plot-1.)*i, T_test, p_test)<<std::endl;  
    // }
    // fout.close();
    
    return 0;
}
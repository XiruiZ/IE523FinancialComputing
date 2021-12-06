//
//  main.cpp
//  project8
//
//  Created by 仲崽 on 2021/11/16.
//  Copyright © 2021年 xx. All rights reserved.
//

#include <random>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <algorithm>


using namespace std;

mt19937 mt_rand(1234);
normal_distribution<double> dis_standard_normal(0.0, 1.0);

double normalCDF(double x) // Phi(-∞, x) aka N(x)
{
    return std::erfc(-x / sqrt(2)) / 2;
}

double p_u, p_d, R, u;
double risk_free_rate, initial_stock_price, expiration_time, volatility, strike_price;
int no_of_divisions;



double option_price_call_black_scholes(const double S,     // spot (underlying) price
                                       const double K,     // strike (exercise) price,
                                       const double r,     // interest rate
                                       const double sigma,   // volatility
                                       const double time)    // time to maturity
{
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+r*time)/(sigma*time_sqrt)+0.5*sigma*time_sqrt;
    double d2 = d1-(sigma*time_sqrt);
    double C = S*normalCDF(d1) - K*exp(-r*time)*normalCDF(d2);
    return C;
};

double option_price_put_black_scholes(const double S,      // spot price
                                      const double K,      // Strike (exercise) price,
                                      const double r,      // interest rate
                                      const double sigma,  // volatility
                                      const double time)
{
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+r*time)/(sigma*time_sqrt) + 0.5*sigma*time_sqrt;
    double d2 = d1-(sigma*time_sqrt);
    return K*exp(-r*time)*normalCDF(-d2) - S*normalCDF(-d1);
};


double max(double a, double b) {
    return (b < a) ? a:b;
}


double european_call_option_dyn_prog(int k, int j, double** V_european_call_option){
    
    if (k == no_of_divisions){ // if it is a terminal node
        int i=j-no_of_divisions;
        V_european_call_option[k][j] = max(0.0 , initial_stock_price * pow(u,i) - strike_price);
    }
    
    else{
        //cout << "test" << V_european_call_option.ncols() << endl;
        if (V_european_call_option[(k+1)][(j+1)]==-1)
            V_european_call_option[(k+1)][(j+1)] = european_call_option_dyn_prog(k+1,j+1,V_european_call_option);
        if (V_european_call_option[(k+1)][j]==-1)
            V_european_call_option[(k+1)][j] = european_call_option_dyn_prog(k+1,j,V_european_call_option);
        if (V_european_call_option[(k+1)][(j-1)]==-1)
            V_european_call_option[(k+1)][(j-1)] = european_call_option_dyn_prog(k+1,j-1,V_european_call_option);
        
        V_european_call_option[k][j] = (p_u*V_european_call_option[(k+1)][(j+1)] + (1.0-p_u-p_d)*V_european_call_option[(k+1)][j] + p_d*V_european_call_option[(k+1)][(j-1)]) / R;
    }
    
    return V_european_call_option[k][j];
}

double get_european_call(){
    double ** V_european_call_option = new double*[(no_of_divisions+1)];
    for (int i=0; i<=no_of_divisions; i++)
        V_european_call_option[i] = new double[(2*no_of_divisions+1)];
    for (int i=0; i<=no_of_divisions; i++)
        for (int j=0; j<=2*no_of_divisions; j++)
            V_european_call_option[i][j]=-1;
    
    double call = european_call_option_dyn_prog(0, no_of_divisions, V_european_call_option);
    
    for (int i=0; i<=no_of_divisions; i++)
        delete [] V_european_call_option[i];
    delete [] V_european_call_option;
    
    return call;
}


double european_put_option_dyn_prog(int k, int j, double** V_european_put_option){

    if (k == no_of_divisions){ // if it is a terminal node
        int i=j-no_of_divisions;
        V_european_put_option[k][j] = max(0.0 , strike_price - initial_stock_price * pow(u,i));
    }
    
    else{
        //cout << "test" << V_european_call_option.ncols() << endl;
        if (V_european_put_option[(k+1)][(j+1)]==-1)
            V_european_put_option[(k+1)][(j+1)] = european_put_option_dyn_prog(k+1,j+1,V_european_put_option);
        if (V_european_put_option[(k+1)][j]==-1)
            V_european_put_option[(k+1)][j] = european_put_option_dyn_prog(k+1,j,V_european_put_option);
        if (V_european_put_option[(k+1)][(j-1)]==-1)
            V_european_put_option[(k+1)][(j-1)] = european_put_option_dyn_prog(k+1,j-1,V_european_put_option);
        
        V_european_put_option[k][j] = (p_u*V_european_put_option[(k+1)][(j+1)] + (1.0-p_u-p_d)*V_european_put_option[(k+1)][j] + p_d*V_european_put_option[(k+1)][(j-1)]) / R;
    }
    
    return V_european_put_option[k][j];
}

double get_european_put(){
    double **V_european_put_option = new double*[(no_of_divisions+1)];
    for (int i=0; i<=no_of_divisions; i++)
        V_european_put_option[i] = new double[(2*no_of_divisions+1)];
    for (int i=0; i<=no_of_divisions; i++)
        for (int j=0; j<=2*no_of_divisions; j++)
            V_european_put_option[i][j]=-1;
    
    double put = european_put_option_dyn_prog(0, no_of_divisions, V_european_put_option);
    
    for (int i=0; i<=no_of_divisions; i++)
        delete [] V_european_put_option[i];
    delete [] V_european_put_option;
    
    return put;
}



double american_call_option_dyn_prog(int k, int j, double** V_american_call_option){
    
    int i=j-no_of_divisions;
    
    if (k == no_of_divisions){ // if it is a terminal node
        V_american_call_option[k][j] = max(0.0 , initial_stock_price * pow(u,i) - strike_price);
    }
    
    else{
        //cout << "test" << V_european_call_option.ncols() << endl;
        if (V_american_call_option[(k+1)][(j+1)]==-1)
            V_american_call_option[(k+1)][(j+1)] = american_call_option_dyn_prog(k+1,j+1,V_american_call_option);
        if (V_american_call_option[(k+1)][j]==-1)
            V_american_call_option[(k+1)][j] = american_call_option_dyn_prog(k+1,j,V_american_call_option);
        if (V_american_call_option[(k+1)][(j-1)]==-1)
            V_american_call_option[(k+1)][(j-1)] = american_call_option_dyn_prog(k+1,j-1,V_american_call_option);
        
        double discounted_one_step_forward_value = (p_u*V_american_call_option[(k+1)][(j+1)] + (1.0-p_u-p_d)*V_american_call_option[(k+1)][j] + p_d*V_american_call_option[(k+1)][(j-1)]) / R;
        double value_if_option_is_exercised_now = max(0, initial_stock_price*pow(u,i)-strike_price);
        
        V_american_call_option[k][j] = max(discounted_one_step_forward_value, value_if_option_is_exercised_now);
    }
    
    return V_american_call_option[k][j];
}

double get_american_call(){
    double ** V_american_call_option = new double*[(no_of_divisions+1)];
    for (int i=0; i<=no_of_divisions; i++)
        V_american_call_option[i] = new double[(2*no_of_divisions+1)];
    for (int i=0; i<=no_of_divisions; i++)
        for (int j=0; j<=2*no_of_divisions; j++)
            V_american_call_option[i][j]=-1;
    
    double call = american_call_option_dyn_prog(0, no_of_divisions, V_american_call_option);
    
    for (int i=0; i<=no_of_divisions; i++)
        delete [] V_american_call_option[i];
    delete [] V_american_call_option;
    
    return call;
}


double american_put_option_dyn_prog(int k, int j, double** V_american_put_option){
    
    int i=j-no_of_divisions;
    
    if (k == no_of_divisions){ // if it is a terminal node
        V_american_put_option[k][j] = max(0.0 , strike_price - initial_stock_price * pow(u,i));
    }
    
    else{
        //cout << "test" << V_european_call_option.ncols() << endl;
        if (V_american_put_option[(k+1)][(j+1)]==-1)
            V_american_put_option[(k+1)][(j+1)] = american_put_option_dyn_prog(k+1,j+1,V_american_put_option);
        if (V_american_put_option[(k+1)][j]==-1)
            V_american_put_option[(k+1)][j] = american_put_option_dyn_prog(k+1,j,V_american_put_option);
        if (V_american_put_option[(k+1)][(j-1)]==-1)
            V_american_put_option[(k+1)][(j-1)] = american_put_option_dyn_prog(k+1,j-1,V_american_put_option);
        
        double discounted_one_step_forward_value = (p_u*V_american_put_option[(k+1)][(j+1)] + (1.0-p_u-p_d)*V_american_put_option[(k+1)][j] + p_d*V_american_put_option[(k+1)][(j-1)]) / R;
        double value_if_option_is_exercised_now = max(0, strike_price-initial_stock_price*pow(u,i));
        
        V_american_put_option[k][j] = max(discounted_one_step_forward_value, value_if_option_is_exercised_now);
    }
    
    return V_american_put_option[k][j];
}

double get_american_put(){
    double **V_american_put_option = new double*[(no_of_divisions+1)];
    for (int i=0; i<=no_of_divisions; i++)
        V_american_put_option[i] = new double[(2*no_of_divisions+1)];
    for (int i=0; i<=no_of_divisions; i++)
        for (int j=0; j<=2*no_of_divisions; j++)
            V_american_put_option[i][j]=-1;
    
    double put = american_put_option_dyn_prog(0, no_of_divisions, V_american_put_option);
    
    for (int i=0; i<=no_of_divisions; i++)
        delete [] V_american_put_option[i];
    delete [] V_american_put_option;
    
    return put;
}



int main(int argc, const char * argv[]) {
//    cout << "Please use cin to input the variables. Thank you! ^v^ \n";
//    cout << "Expiration:";
//    cin >> expiration_time;
//    cout << "Number of Stages:";
//    cin >> no_of_divisions;
//    cout << "Risk-free rate:";
//    cin >> risk_free_rate;
//    cout << "Volatility:";
//    cin >> volatility;
//    cout << "Initial Price:";
//    cin >> initial_stock_price;
//    cout << "Strike Price:";
//    cin >> strike_price;
    
    sscanf (argv[1], "%lf", &expiration_time);
    sscanf (argv[2], "%d", &no_of_divisions);
    sscanf (argv[3], "%lf", &risk_free_rate);
    sscanf (argv[4], "%lf", &volatility);
    sscanf (argv[5], "%lf", &initial_stock_price);
    sscanf (argv[6], "%lf", &strike_price);
    
    //volatility = volatility/100.0;
    
    R = exp(risk_free_rate * expiration_time / ((float) no_of_divisions));
    u = exp(volatility * sqrt(2 * (expiration_time / ((float) no_of_divisions))));
    p_u = pow((sqrt(R) - 1.0/sqrt(u)) / (sqrt(u) - 1.0/sqrt(u)) , 2);
    p_d = pow((sqrt(u) - sqrt(R)) / (sqrt(u) - 1.0/sqrt(u)) , 2);
    
    
    double trinomial_European_call, trinomial_European_put, BS_call, BS_put;
    double trinomial_American_call, trinomial_American_put;
    
    BS_call = option_price_call_black_scholes(initial_stock_price, strike_price, risk_free_rate, volatility, expiration_time);
    BS_put = option_price_put_black_scholes(initial_stock_price, strike_price, risk_free_rate, volatility, expiration_time);
    
    trinomial_European_call = get_european_call();
    trinomial_European_put = get_european_put();
    
    trinomial_American_call = get_american_call();
    trinomial_American_put = get_american_put();
    
    cout << "---------------------------" << endl;
    cout << "Expiration Time (Years) = " << expiration_time << "\n" << "Number of Divisions = " << no_of_divisions << "\n" << "Risk Free Interest Rate = " << risk_free_rate << "\n" << "Volatility (% age of stock value) = " << (int)(volatility*100) << "\n" << "Initial Stock Price = " << initial_stock_price << "\n" << "Strike Price = " << strike_price << endl;
    cout << "---------------------------" << endl;
    cout << "Up Factor = " << u << endl;
    cout << "Uptick Probablity = " << p_u << endl;
    cout << "Notick Probablity = " << 1.0-p_u-p_d << endl;
    cout << "Downtick Probablity = " << p_d << endl;
    cout << "---------------------------" << endl;
    cout << "Trinomial Price of an European Call Option = " << trinomial_European_call << endl;
    cout << "Call Price according to Black-Scholes = " << BS_call << endl;
    cout << "---------------------------" << endl;
    cout << "Trinomial Price of an European Put Option = " << trinomial_European_put << endl;
    cout << "Put Price according to Black-Scholes = " << BS_put << endl;
    cout << "---------------------------" << endl;
    cout << "Verifying Put-Call Parity: S+P-C = K * exp(-r*T)" << endl;
    cout << initial_stock_price << " + " << trinomial_European_put << " - " << trinomial_European_call << " = " << strike_price << " * exp(-" << risk_free_rate << "*" << expiration_time << ")" << endl;
    cout << initial_stock_price+trinomial_European_put-trinomial_European_call << " = " << strike_price * exp(-risk_free_rate*expiration_time) << endl;
    cout << "---------------------------" << endl;
    cout << "Trinomial Price of an American Call Option = " << trinomial_American_call << endl;
    cout << "Trinomial Price of an American Put Option = " << trinomial_American_put << endl;
    cout << "---------------------------" << endl;
    cout << "Let Us Verify the Put-Call Parity: S+P-C = K*exp(-r*T)" << endl;
    cout << initial_stock_price << " + " << trinomial_American_put << " - " << trinomial_American_call << " ?= " << strike_price << " * exp(-" << risk_free_rate << "*" << expiration_time << ")" << endl;
    if ((initial_stock_price+trinomial_American_put-trinomial_American_call) == (strike_price * exp(-risk_free_rate*expiration_time))){
         cout << initial_stock_price+trinomial_American_put-trinomial_American_call << " = " << strike_price * exp(-risk_free_rate*expiration_time) << endl;
        cout << "The Put-Call Parity Holds!" << endl;
    }
    else{
        cout << initial_stock_price+trinomial_American_put-trinomial_American_call << " != " << strike_price * exp(-risk_free_rate*expiration_time) << endl;
        cout << "Seems That The Put-Call Parity Does Not Hold!" << endl;
    }
    
    return 0;
}

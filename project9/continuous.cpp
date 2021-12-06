//
//  main.cpp
//  project9
//
//  Created by 仲崽 on 2021/11/30.
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

//mt19937 mt_rand(1234);
normal_distribution<double> dis_standard_normal(0.0, 1.0);
default_random_engine generator;


double normalCDF(double x) // Phi(-∞, x) aka N(x)
{
    return std::erfc(-x / sqrt(2)) / 2;
}

double get_uniform()
{
    uniform_real_distribution<double> distribution(0.0,1.0);
    double number = distribution(generator);
    return (number);
}

double max(double a, double b) {
    return (b < a) ? a:b;
}

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



double T, r, sigma, S0, K, B;
int n, m;


int main(int argc, char * argv[]) {
    
    sscanf (argv[1], "%lf", &T); // Duration
    sscanf (argv[2], "%lf", &r); // Risk-free interest rate
    sscanf (argv[3], "%lf", &sigma); // Volatility
    sscanf (argv[4], "%lf", &S0); // Initial stock price
    sscanf (argv[5], "%lf", &K); // Strike Price
    sscanf (argv[6], "%d", &n); // the number of trials/repetitions of the Monte Carlo simulations
    sscanf (argv[7], "%d", &m); // the number of sample-points in each sample price-path
    sscanf (argv[8], "%lf", &B); // Barrier Price
    
    double delta_T = T/((double) m);
    double delta_R = (r - 0.5*pow(sigma,2))*delta_T;
    double delta_SD = sigma*sqrt(delta_T);
    
    double call_option_price=0.0;
    double put_option_price=0.0;
    double adjusted_call_option_price=0.0;
    double adjusted_put_option_price=0.0;
    
    double theoretical_call_price, theoretical_put_price;
    double lambda = (r + pow(sigma,2) / 2.0) / pow(sigma,2);

    
    
    // by sharing random variables we create 4 paths
    
    for (int i = 0; i < n; i++){
        
        double current_stock_price[4] = {S0, S0, S0, S0};
        double adjust_stock_price[4] = {S0, S0, S0, S0};
        double not_knocked[4]={1,1,1,1};
        
        for (int j = 0; j < m; j++)
        {

            // create the unit normal variates using the Box-Muller Transform
            double x = get_uniform();
            double y = get_uniform();
            double a =  sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
            double b =  sqrt(-2.0*log(x)) * sin(6.283185307999998*y);
            
            adjust_stock_price[0] = adjust_stock_price[0]*exp(delta_R + delta_SD*a);
            current_stock_price[0] = current_stock_price[0]*exp(delta_R + delta_SD*a);
            if (current_stock_price[0] <= B){
                current_stock_price[0] = 0;
                not_knocked[0]=0;
            }
            

            adjust_stock_price[1] = adjust_stock_price[1]*exp(delta_R - delta_SD*a);
            current_stock_price[1] = current_stock_price[1]*exp(delta_R - delta_SD*a);
            if (current_stock_price[1] <= B){
                current_stock_price[1] = 0;
                not_knocked[1]=0;
            }
            
            adjust_stock_price[2] = adjust_stock_price[2]*exp(delta_R + delta_SD*b);
            current_stock_price[2] = current_stock_price[2]*exp(delta_R + delta_SD*b);
            if (current_stock_price[2] <= B){
                current_stock_price[2] = 0;
                not_knocked[2]=0;
            }
            
            adjust_stock_price[3] = adjust_stock_price[3]*exp(delta_R - delta_SD*b);
            current_stock_price[3] = current_stock_price[3]*exp(delta_R - delta_SD*b);
            if (current_stock_price[3] <= B){
                current_stock_price[3] = 0;
                not_knocked[3]=0;
            }
        }
        
        double sum_call=0.0;
        double sum_put=0.0;
        double adj_sum_call=0.0;
        double adj_sum_put=0.0;
        
        for (int j = 0; j < 4; j++)
            sum_call += max(0.0, current_stock_price[j] - K) * not_knocked[j];
        call_option_price += sum_call / 4.0;
        
        for (int j = 0; j < 4; j++)
            sum_put += max(0.0, K - current_stock_price[j]) * not_knocked[j];
        put_option_price += sum_put / 4.0;

        double p_c[4];
        for (int j = 0; j < 4; j++){
            if (S0 <= B || adjust_stock_price[j] <= B)
                p_c[j] = 1.0;
            else
                p_c[j]=exp(-2*log(S0/B) * log(adjust_stock_price[j]/B) / (pow(sigma,2) * T));
            
        }
        
        for (int j = 0; j < 4; j++)
            adj_sum_call += max(0.0, adjust_stock_price[j] - K) * (1-p_c[j]);
        adjusted_call_option_price += adj_sum_call / 4.0;
        
        for (int j = 0; j < 4; j++)
            adj_sum_put += max(0.0, K - adjust_stock_price[j]) * (1-p_c[j]);
        adjusted_put_option_price += adj_sum_put / 4.0;
        
    }
    
    call_option_price = exp(-r*T) * (call_option_price / ((double) n));
    adjusted_call_option_price = exp(-r*T) * (adjusted_call_option_price / ((double) n));
    put_option_price = exp(-r*T) * (put_option_price / ((double) n));
    adjusted_put_option_price = exp(-r*T) * (adjusted_put_option_price / ((double) n));
    
    //Theoretical Prices
    if (B <= K) {
        double y = log(pow(B,2)/S0/K) / sigma / sqrt(T) + lambda * sigma * sqrt(T);
        double c_di= S0 * pow((B/S0), 2*lambda) * normalCDF(y) - K * exp(-1*r*T) * pow((B/S0), 2*lambda-2) * normalCDF(y - sigma*sqrt(T));
        double c = option_price_call_black_scholes(S0, K, r, sigma, T);
        theoretical_call_price = c - c_di;
    }
    else {
        double x1 = log(S0/B) / sigma / sqrt(T) + lambda * sigma * sqrt(T);
        double y1 = log(B/S0) / sigma / sqrt(T) + lambda * sigma * sqrt(T);
        double c_do = S0 * normalCDF(x1) -  K * exp(-1*r*T) * normalCDF(x1-sigma*sqrt(T)) - S0 * pow((B/S0), 2*lambda) * normalCDF(y1) + K * exp(-1*r*T) * pow((B/S0), 2*lambda-2) * normalCDF(y1-sigma*sqrt(T));
        theoretical_call_price = c_do;
    }

    if (B <= K) {
        double y = log(pow(B,2)/S0/K) / sigma / sqrt(T) + lambda * sigma * sqrt(T);
        double x1 = log(S0/B) / sigma / sqrt(T) + lambda * sigma * sqrt(T);
        double y1 = log(B/S0) / sigma / sqrt(T) + lambda * sigma * sqrt(T);
        double p_di = - S0 * normalCDF(-x1) + K * exp(-1*r*T)*normalCDF(-x1+sigma*sqrt(T)) + S0 * pow((B/S0), 2*lambda) * (normalCDF(y) - normalCDF(y1)) - K * exp(-1*r*T) * pow((B/S0), 2*lambda-2) * (normalCDF(y-sigma*sqrt(T))-normalCDF(y1-sigma*sqrt(T)));
        double p = option_price_put_black_scholes(S0, K, r, sigma, T);
        theoretical_put_price = p - p_di;
    }
    else {
        theoretical_put_price = 0;
    }
    
    
    
    cout << "-----------------------------" << endl;
    cout << "European Down-and-Out Continuous Barrier Options Pricing via Monte Carlo Simulation" << endl;
    cout << "Expiration Time (Years) = " << T << endl;
    cout << "Risk Free Interest Rate = " << r << endl;
    cout << "Volatility (%age of stock value) = " << sigma*100 << endl;
    cout << "Initial Stock Price = " << S0 << endl;
    cout << "Strike Price = " << K << endl;
    cout << "Barrier Price = " << B << endl;
    cout << "Number of Trials = " << n << endl;
    cout << "Number of Divisions = " << m << endl;
    cout << "-----------------------------" << endl;
    cout << "-----------------------------" << endl;
    cout << "The average Call Price by explicit simulation = " << call_option_price<< endl;
    cout << "The call price using the (1-p)-adjustment term = " << adjusted_call_option_price<< endl;
    cout << "Theoretical Call Price = " << theoretical_call_price << endl;
    cout << "-----------------------------" << endl;
    cout << "The average Put Price by explicit simulation = " << put_option_price<< endl;
    cout << "The put price using the (1-p)-adjustment term = " <<adjusted_put_option_price << endl;
    cout << "Theoretical Put Price = " << theoretical_put_price<< endl;
    cout << "-----------------------------" << endl;
    
    

    return 0;
}

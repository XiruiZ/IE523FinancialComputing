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

//mt19937 mt_rand(10);
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




double T, r, sigma, S0, K, B;
int n, m;


int main(int argc, char * argv[]) {
    
    sscanf (argv[1], "%lf", &T); // Duration
    sscanf (argv[2], "%lf", &r); // Risk-free interest rate
    sscanf (argv[3], "%lf", &sigma); // Volatility
    sscanf (argv[4], "%lf", &S0); // Initial stock price
    sscanf (argv[5], "%lf", &K); // Strike Price
    sscanf (argv[6], "%d", &n); // the number of trials/repetitions of the Monte Carlo simulations
    sscanf (argv[7], "%d", &m); // the number of (equally-spaced) discrete barriers from 0 to T
    sscanf (argv[8], "%lf", &B); // Barrier Price
    
    double delta_T = T/((double) m);
    double delta_R = (r - 0.5*pow(sigma,2))*delta_T;
    double delta_SD = sigma*sqrt(delta_T);
    
    double call_option_price=0.0;
    double put_option_price=0.0;
    double adjusted_call_option_price=0.0;
    double adjusted_put_option_price=0.0;
    
    
    for (int i = 0; i < n; i++){
        
        double current_stock_price[4] = {S0, S0, S0, S0};
        double adjust_stock_price[4] = {S0, S0, S0, S0};
        double not_knocked[4]={1,1,1,1};
        double p_d[4]={1,1,1,1};
        
        for (int j = 0; j < m; j++)
        {
            
            // create the unit normal variates using the Box-Muller Transform
            double x = get_uniform();
            double y = get_uniform();
            double a =  sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
            double b =  sqrt(-2.0*log(x)) * sin(6.283185307999998*y);
            
            current_stock_price[0] = current_stock_price[0]*exp(delta_R + delta_SD*a);
            adjust_stock_price[0] = adjust_stock_price[0]*exp(delta_R + delta_SD*a);
            if (current_stock_price[0] <= B){
                current_stock_price[0] = 0;
                not_knocked[0]=0;
            }
            
            current_stock_price[1] = current_stock_price[1]*exp(delta_R - delta_SD*a);
            adjust_stock_price[1] = adjust_stock_price[1]*exp(delta_R - delta_SD*a);
            if (current_stock_price[1] <= B){
                current_stock_price[1] = 0;
                not_knocked[1]=0;
            }
            
            current_stock_price[2] = current_stock_price[2]*exp(delta_R + delta_SD*b);
            adjust_stock_price[2] = adjust_stock_price[2]*exp(delta_R + delta_SD*b);
            if (current_stock_price[2] <= B){
                current_stock_price[2] = 0;
                not_knocked[2]=0;
            }
            
            current_stock_price[3] = current_stock_price[3]*exp(delta_R - delta_SD*b);
            adjust_stock_price[3] = adjust_stock_price[3]*exp(delta_R - delta_SD*b);
            if (current_stock_price[3] <= B){
                current_stock_price[3] = 0;
                not_knocked[3]=0;
            }
            
        }
        
        for (int j = 1; j < m; j++){
            double mu[4];
            for (int k = 0; k < 4; k++) {
                mu[k] = S0 + ((double) j)/((double) m) * (adjust_stock_price[k]-S0);
            }
            double sigma = sqrt(((double) j) * ((double) m - j)
                                /((double) pow(m,2)));
            for (int k = 0; k < 4; k++) {
                p_d[k] = p_d[k] * (1 - normalCDF((B - mu[k])/sigma));
            }
        }
        
        double sum_call=0.0;
        double sum_put=0.0;
        double adj_sum_call=0.0;
        double adj_sum_put=0.0;
        
        for (int j = 0; j < 4; j++)
            sum_call += not_knocked[j] * max(0.0, current_stock_price[j] - K);
        call_option_price += sum_call / 4.0;
        
        for (int j = 0; j < 4; j++)
            sum_put += not_knocked[j] * max(0.0, K - current_stock_price[j]);
        put_option_price += sum_put / 4.0;
        /*
        double p_c[4];
        for (int j = 0; j < 4; j++){
            if (S0 <= B || adjust_stock_price[j] <= B)
                p_c[j] = 1.0;
            else
                p_c[j]=exp(-2*log(S0/B)* log(adjust_stock_price[j]/B)/(pow(sigma,2) * T));
            
        }
        */
        for (int j = 0; j < 4; j++){
            if (S0 <= B || adjust_stock_price[j] <= B)
                p_d[j] = 0;
            adj_sum_call += max(0.0, adjust_stock_price[j] - K) * p_d[j];
        }
        adjusted_call_option_price += adj_sum_call / 4.0;
        
        for (int j = 0; j < 4; j++){
            if (S0 <= B || adjust_stock_price[j] <= B)
                p_d[j] = 0;
            adj_sum_put += max(0.0, K - adjust_stock_price[j]) * p_d[j];
        }
        adjusted_put_option_price += adj_sum_put / 4.0;
        
    }
    
    call_option_price = exp(-r*T) * (call_option_price / ((double) n));
    adjusted_call_option_price = exp(-r*T) * (adjusted_call_option_price / ((double) n));
    put_option_price = exp(-r*T) * (put_option_price / ((double) n));
    adjusted_put_option_price = exp(-r*T) * (adjusted_put_option_price / ((double) n));
    

    
    
    cout << "-----------------------------" << endl;
    cout << "Expiration Time (Years) = " << T << endl;
    cout << "Risk Free Interest Rate = " << r << endl;
    cout << "Volatility (%age of stock value) = " << sigma*100 << endl;
    cout << "Initial Stock Price = " << S0 << endl;
    cout << "Strike Price = " << K << endl;
    cout << "Barrier Price = " << B << endl;
    cout << "Number of Trials = " << n << endl;
    cout << "Number of Discrete Barriers = " << m << endl;
    cout << "-----------------------------" << endl;
    cout << "The average Call Price via explicit simulation of price paths = " << call_option_price<< endl;
    cout << "The average Call Price with Brownian-Bridge correction on the final price = " << adjusted_call_option_price<< endl;
    cout << "The average Put Price via explicit simulation = " << put_option_price<< endl;
    cout << "The average Put Price with Brownian-Bridge correction on the final price = " << adjusted_put_option_price<< endl;
    cout << "-----------------------------" << endl;
    
    
    
    return 0;
}


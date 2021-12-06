// Written by Prof. Sreenivas for IE523: Financial Computing

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <numeric>

#include "lp_lib.h"

using namespace std;

const double ERROR_ = 1e-10;
int number_of_cash_flows;
vector <double> price_list;
vector <int> maturity_list;
vector <double> yield_to_maturity;
vector <double> duration;
vector <double> convexity;
double debt_obligation_amount;
double time_when_debt_is_due;
vector <double> percentage_of_cash_flow_to_meet_debt_obligation;

string input;

double func(vector <double> cash_flow, double price, int maturity, double rate)
{
    double temp=0.0;
    for (int i=0; i<maturity; i++){
        // t = i+1
        temp += cash_flow[i] * pow(1+rate, maturity-(i+1));
    }
    double f_r = price * pow(1+rate, maturity) - temp;
    return f_r;
}

double derivative_function(vector <double> cash_flow, double price, int maturity, double rate)
{
    double temp=0.0;
    for (int i=0; i<maturity-1; i++){
        // t = i+1
        temp += cash_flow[i] * (maturity-(i+1)) * pow(1+rate, maturity-(i+1)-1);
    }
    double f_r_prime = maturity * price * pow(1+rate, maturity-1) - temp;
    return f_r_prime;
}

double Newton_Raphson(vector <double> cash_flow, double price, int maturity, double rate)
{
    double r = rate;
    double f_r = func(cash_flow, price, maturity, r);
    double f_r_prime = derivative_function(cash_flow, price, maturity, r);
    while (abs(f_r) > ERROR_){
        r -= f_r/f_r_prime;
        f_r = func(cash_flow, price, maturity, r);
        f_r_prime = derivative_function(cash_flow, price, maturity, r);
    }
    return r;
}

double get_duration(vector <double> cash_flow, double price, int maturity, double rate)
{
    double N=0.0;
    for (int i=0; i<maturity; i++){
        // t = i+1
        N += (i+1) * cash_flow[i] / pow(1+rate, i+1);
    }
    N = N / price;
    return N;
}

double get_convexity(vector <double> cash_flow, double price, int maturity, double rate)
{
    double C=0.0;
    for (int i=0; i<maturity; i++){
        // t = i+1
        C += (i+1) * (i+1+1) * cash_flow[i] / pow(1+rate, i+1+2);
    }
    C = C / price;
    return C;
}

double present_value_of_debt()
{
    double V_0;
    // using the average-value-of-the-YTMs
    double r = accumulate(yield_to_maturity.begin(), yield_to_maturity.end(), 0.0) / yield_to_maturity.size();
    V_0 = debt_obligation_amount / pow(1+r, time_when_debt_is_due);
    return V_0;
}

//void print_data(char *filename)
void print_data(string filename)
{
	cout << "Input File: " << filename << endl;
	cout << "We owe " << debt_obligation_amount << " in " << time_when_debt_is_due << " years" << endl;
	cout << "Number of Cash Flows: " << number_of_cash_flows << endl;
	for (int i = 0; i < number_of_cash_flows; i++) 
	{
		cout << "---------------------------" << endl;
		cout << "Cash Flow #" << i+1 << endl;
		cout << "Price = " << price_list[i] << endl;
		cout << "Maturity = " << maturity_list[i] << endl;
		cout << "Yield to Maturity = " << yield_to_maturity[i] << endl;
		cout << "Duration = " << duration[i] << endl;
		cout << "Convexity = " << convexity[i] << endl;
		cout << "Percentage of Face Value that would meet the obligation = " << 
		percentage_of_cash_flow_to_meet_debt_obligation[i] << endl;
	}
	cout << "---------------------------" << endl;
}

void get_data(char* argv[])
{
    
    cout << "Input file: ";
    cin >> input;
    ifstream input_filename(input);
    
    if (input_filename.is_open()) {
        input_filename >> number_of_cash_flows;
        
        for (int i=0; i<number_of_cash_flows; i++){
            vector <double> cash_flow;
            double price, ytm, d, c;
            int maturity;
            
            input_filename >> price;
            input_filename >> maturity;
            
            string str;
            for (int j=0; j<maturity; j++)
            {
                double cf;
                input_filename >> cf;
                cash_flow.push_back(cf);
            }
            
            // guess from 0
            ytm = Newton_Raphson(cash_flow, price, maturity, 0);
            d = get_duration(cash_flow, price, maturity, ytm);
            c = get_convexity(cash_flow, price, maturity, ytm);
            
            price_list.push_back(price);
            maturity_list.push_back(maturity);
            yield_to_maturity.push_back(ytm);
            duration.push_back(d);
            convexity.push_back(c);
            
        }
        
        input_filename >> debt_obligation_amount;
        input_filename >> time_when_debt_is_due;
        
        double V_0=present_value_of_debt();
        for (int i=0; i<number_of_cash_flows; i++){
            double p = V_0 / price_list[i];
            percentage_of_cash_flow_to_meet_debt_obligation.push_back(p);
            // ?
        }
            }
    else {
        cout << "Input file missing" << endl;
        exit(0);
    }
}

void get_optimal_portfolio()
{
	lprec *lp;
    
	// Maximize portfolio convexity
    double* solution;
    solution = new double[number_of_cash_flows];
    lp = make_lp(0,number_of_cash_flows); // there are number_of_cash_flows variables
    set_verbose(lp, 3); // keeps the message reporting of lp_solve to a minimum
    
    {
        // sum(lambda_i * D_i) = N
        double* row;
        row = new double[number_of_cash_flows+1];
        for (int i=1; i<=number_of_cash_flows; i++){
            row[i] = duration[i-1];
        }
        add_constraint(lp, row, EQ, time_when_debt_is_due);
    }
    {
        // sum(lambda_i) = 1
        double* row;
        row = new double[number_of_cash_flows + 1];
        for (int i=1; i<=number_of_cash_flows; i++){
            row[i] = 1;
        }
        add_constraint(lp, row, EQ, 1);
    }
    {
        // find lambda that maximize sum(lambda_i * C_i),
        // which minimizes sum(lambda_i * -C_i)
        double* row;
        row = new double[number_of_cash_flows + 1];
        for (int i=1; i<=number_of_cash_flows; i++){
            row[i] = -convexity[i-1];
        }
        set_obj_fn(lp, row);
    }
    
    cout << "***************************" << endl;
    cout << "Model name:" << endl;
    print_lp(lp);
    
    int ret;
    ret = solve(lp);
    
    cout << "Returned Value from Solve is " << ret << endl;
    
    // Check if there is a solution
    if (ret == 0){
        cout << "Largest Convexity we can get is: " << -1*get_objective(lp) << endl;
    
        get_variables(lp, solution);
        cout << "Optimal portfolio:" << endl;
        for (int i = 0; i < number_of_cash_flows; i++)
            cout << "%Cash Flow:" << i+1 << "  " << solution[i] << endl;
    
        cout << "***************************" << endl;
        cout << "To immunize against small change in 'r' for each $1 of PV, you should buy:" << endl;
        for (int i = 0; i < number_of_cash_flows; i++){
            if (solution[i]!=0){
                cout << "$" << solution[i] << "  of  Cash Flow#" << i+1 << endl;
            }
        }
        cout << "If you need to immunize for a large PV-value, just buy an appropriate proportion." << endl;
    
        double PV[] = {500, 750, 1000, present_value_of_debt()};
        for (int i=0; i<4; i++){
            cout << "***************************" << endl;
            cout << "For example, if you want to immunize for $" << PV[i] << " of PV, buy:" << endl;
            for (int i = 0; i < number_of_cash_flows; i++){
                if (solution[i]!=0){
                    cout << "$" << solution[i]*PV[i] << "  of  Cash Flow#" << i+1 << endl;
                }
            }
        }
        cout << "***************************" << endl;
    }
    else{
        cout << "There is no portfolio that meets the duration constraint of " << time_when_debt_is_due << " years." << endl;
    }
    
    delete_lp(lp);
    //delete[] solution;
    //delete[] row;

}
	
int main (int argc, char* argv[])
{
    /*
	if (argc == 1) {
		cout << "Input filename missing" << endl;
	}
	else
     */
	{
        cout << "Please cin the input file name. Thank you! :)" << endl;
		get_data(argv);
		
		//print_data(argv[1]);
        print_data(input);
		
		get_optimal_portfolio();
	}
	return (0);
}


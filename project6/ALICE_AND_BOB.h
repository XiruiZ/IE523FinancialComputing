/*
 *  alice_and_bob.h
 *  Loosing as little as possible
 *
 *  Created by Ramavarapu Sreenivas on 9/2/12.
 *  Copyright 2012 University of Illinois at Urbana-Champaign. All rights reserved.
 *
 */
#ifndef ALICE_AND_BOB
#define ALICE_AND_BOB

#include <cmath>
#include <random>
#include <map>
//#include "matplotlibcpp.h"

using namespace std;
//namespace plt = matplotlibcpp;

mt19937 mt_rand(1234);
#define GAUSSIAN_DENSITY(x) 0.398942280375387*exp(-x*x/2.0)
default_random_engine generator;


class Alice_Bob_Game{
    
private:
	double alice_probability, bob_probability;
    // use map to save results: n_i_map(vector{n,i}) = nCi
    map <vector<int>, int> n_i_map = {};
    // matrix may be quicker?
    // but how to distribute memory?
	
	// private member function: uniform RV generator
	double get_uniform()
	{
        uniform_real_distribution<double> distribution(0.0,1.0);
        double number = distribution(generator);
        return (number);
	}
	
	// private member function: nCi (i.e. n-take-i) 
	int take(int n, int i)
	{
        if (i>n || n<=0)
            return 0;
        
        vector<int> temp = {n, i};
        auto it = n_i_map.find(temp); // check if nCi has already been calculated
        
        if ( it == n_i_map.end()  ) {
            // not found
            int nCi;
            if (i==0 || i==n)
                nCi = 1;
            else if (i==1 || i==(n-1))
                nCi = n;
            else // now 2<=i<=n-2
                nCi = take(n-1, i-1) + take(n-1, i);
            n_i_map.insert({temp, nCi});

            return nCi;
        }
        else { // found
            return  it->second; // return n_i_map(temp)
        }
	}
	
	// this routine implements the probability that Alice has more 
	// heads than Bob after n-many coin tosses
	double theoretical_value(double q, double p, int n)
	{
		// implement equation 1.1 of Addona-Wagon-Wilf paper
        double f=0;
        for (int r=0; r<n; r++){
            double tem=0;
            for (int s=r+1; s<=n; s++){
                tem = tem + take(n, s) * pow(q, s) * pow(1.0-q, n-s);
            }
            f = f + take(n, r) * pow(p, r) * pow(1.0-p, n-r) * tem;
        }
        return f;
	}

public: 
	// public function:
    
	void set_probability(double alice_p, double bob_p)
	{
		alice_probability = alice_p;
		bob_probability = bob_p;
	}
	
	// probability of Alice winning the game.
	double simulated_value(int number_of_coin_tosses_in_each_game, int no_of_trials)
	{
		int no_of_wins_for_alice = 0;
		for (int i = 0; i < no_of_trials; i++) 
		{
			int number_of_heads_for_alice = 0;
			int number_of_heads_for_bob = 0;
			for (int j = 0; j < number_of_coin_tosses_in_each_game; j++) 
			{
				if (get_uniform() < alice_probability) 
					number_of_heads_for_alice++;
				if (get_uniform() < bob_probability)
					number_of_heads_for_bob++;
			}
			if (number_of_heads_for_alice > number_of_heads_for_bob)
				no_of_wins_for_alice++;
		}
		return (((double) no_of_wins_for_alice)/((double) no_of_trials));
	}
		
	int search_result()
	{
		// implememt a discrete-search procedure for the optimal n-value. 
		// start with n = 1 and find the discrete-value of n that has 
		// the largest probability for Alice winning.  Why would this work?
		// See Theorem 2.2 of the paper for the reason!
        double p1, p2;
        int n = 1;
        p1 = theoretical_value(alice_probability, bob_probability, 1);
        p2 = theoretical_value(alice_probability, bob_probability, 2);
        
        while(p2>p1){
            n += 1;
            p1 = p2;
            p2 = theoretical_value(alice_probability, bob_probability, n+1);
        }
        
        return n;
	}
    
    
    void output(int n, int no_of_trials){
        
        ofstream myFile("output.csv");
        myFile << "N, P_experiment, P_theory" << "\n";
        
        //vector<double> P_experiment;
        //vector<double> P_theory;
        //vector<int> N;
        for (int i=1; i<=n; i++){
            double p_experiment = simulated_value(i, no_of_trials);
            //P_experiment.push_back(p_experiment);
            double p_theory = theoretical_value(alice_probability, bob_probability, i);
            //P_theory.push_back(p_theory);
            //N.push_back(i);
            myFile << i << "," << p_experiment << "," << p_theory << "\n";
        }
        
        myFile.close();
        
        /*
        plt::plot(N, P_experiment);
        plt::plot(N, P_theory,"r-");
        plt::title("Alice's Winning Rate VS n");
        plt::legend();
        plt::show();
    
     */
    }

};
#endif










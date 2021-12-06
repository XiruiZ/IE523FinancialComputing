// IE523: Financial Computation
// "How to lose as little as possible" by Addona, Wagon and Wilf
// Written by Prof. Sreenivas
// 

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include "ALICE_AND_BOB.h"
using namespace std;
	
int main (int argc, const char* argv[])
{
	Alice_Bob_Game x;
    cout << "Please cin p and q. Thank you! :)" << endl;
	//sscanf (argv[1], "%lf", &alice_success_prob);
	//sscanf (argv[2], "%lf", &bob_success_prob);
	
    double alice_success_prob, bob_success_prob;
    cout << "Probability of success for Alice = ";
    cin >> alice_success_prob;
    cout << "Probability of success for Bob = ";
	cin >> bob_success_prob;
	
	x.set_probability(alice_success_prob, bob_success_prob);
	
	int optimal = x.search_result();
    if (optimal > 0){
		cout << "The optimal number of coin tosses in each game is " << optimal << endl;
        int no_of_trials = pow(10, 6);
        x.output(optimal+5, no_of_trials);
    }
    
	else {
		cout << "The optimal number of coin tosses in each game exceeds 100... Quitting" << endl;
	}
}
	
	
	

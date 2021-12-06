//
//  main.cpp
//  project7
//
//  Created by 仲崽 on 2021/10/26.
//  Copyright © 2021年 xx. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <random>

#include "/Users/zxr/workspace/newmat11/newmatap.h"
#include "/Users/zxr/workspace/newmat11/newmat.h"
#include "/Users/zxr/workspace/newmat11/newmatio.h"

using namespace std;

//unsigned seed = (unsigned) std::chrono::system_clock::now().time_since_epoch().count();
mt19937 mt_rand(1234);
default_random_engine generator;


double get_uniform()
{
    uniform_real_distribution<double> distribution(0.0,1.0);
    double number = distribution(generator);
    return (number);
}

// This routine fills a n x n matrix, with random entries
Matrix create_random_matrix(int n)

{
    Matrix A(n, n);
    for (int i = 1; i <= n; i++){
        for (int j = 1; j <= n; j++){
            A(i,j) = 5.0 * get_uniform();
            double p_or_n = get_uniform();
            if (p_or_n < 0.5)
                A(i,j) = -1 * A(i,j);
        }
    }
    return A;
}

void print_matrix(Matrix A, int n)
{
    int i, j;
    
    //cout << endl;
    cout << "---------------";
    
    for (i = 1; i <= n; i++) {
        cout << endl;
        cout << "|  " ;
        for (j = 1; j <= n; j++)
            cout << A(i,j) << "  ";
        cout << "|";
    }
    cout << endl;
    cout << "---------------";
    cout << endl;
}


Matrix repeated_squaring(Matrix A, int exponent, int n){
    if (0==exponent){
        IdentityMatrix I(n);
        return I;
    }
    if (1==exponent%2){ // k is odd
        return (A * repeated_squaring(A*A, exponent/2, n));
    }
    else{
        return (repeated_squaring(A*A, exponent/2, n));
    }
}

Matrix direct_multiplication(Matrix A, int exponent, int n){
    if (0==exponent){
        IdentityMatrix I(n);
        return I;
    }
    if (1==exponent){
        return A;
    }
    return (A * direct_multiplication(A, exponent-1, n));
}

void output_exp(Matrix A, int n){
    string filename="";
    filename = "output_exp_" + to_string(n) + ".csv";
    ofstream myFile(filename);
    myFile << "exponent, t_repeated_squaring, t_Brute_Force" << "\n";
    
    for (int k=0; k<=100; k++){
        float t_repeated_squaring, t_Brute_Force;
        
        float time_before = clock();
        Matrix B = repeated_squaring(A, k, n);
        float time_after = clock();
        float diff = time_after - time_before;
        t_repeated_squaring = diff/CLOCKS_PER_SEC;
        
        time_before = clock();
        Matrix C = direct_multiplication(A, k, n);
        time_after = clock();
        diff = time_after - time_before;
        t_Brute_Force = diff/CLOCKS_PER_SEC;
        
        myFile << k << "," << t_repeated_squaring << "," << t_Brute_Force << "\n";
        
        //if (k%10==0)
            //cout << k << endl;
    }
    
    myFile.close();
}

void output_size(int exponent){
    ofstream myFile("output_size.csv");
    myFile << "Matrix Size, t_repeated_squaring, t_Brute_Force" << "\n";
    
    float t_repeated_squaring, t_Brute_Force;
    float time_before, time_after, diff;
    for (int i=1; i<=1000; i+=10){
        Matrix A(i, i);
        A = create_random_matrix(i);
        
        time_before = clock();
        Matrix B = repeated_squaring(A, exponent, i);
        time_after = clock();
        diff = time_after - time_before;
        t_repeated_squaring = diff/CLOCKS_PER_SEC;
        
        time_before = clock();
        Matrix C = direct_multiplication(A, exponent, i);
        time_after = clock();
        diff = time_after - time_before;
        t_Brute_Force = diff/CLOCKS_PER_SEC;
        
        myFile << i << "," << t_repeated_squaring << "," << t_Brute_Force << "\n";
        
        //if (i%100 == 1)
            //cout << i << endl;
    }
    
    myFile.close();
}



int main(int argc, const char * argv[]) {
    
    int n, k;
    cout << "Please use cin to input the size and exponent. Thank you! ^_^ \n";
    cout << "The size of matrix is:";
    cin >> n;
    cout << "The exponent is:";
    cin >> k;
    
    Matrix A(n, n);
    A = create_random_matrix(n);
    
    cout << "------------------------------" << endl;
    cout << "The number of rows/columns in the square matrix is: " << n << endl;
    //cout << endl;
    cout << "The exponent is: " << k << endl;
    //cout << "------------------------------" << endl;
    cout << endl;
    cout << "Repeated Squaring Result:" << endl;
    float time_before = clock();
    Matrix B = repeated_squaring(A, k, n);
    float time_after = clock();
    float diff = time_after - time_before;
    print_matrix(B , n);
    //cout << endl;
    cout << "It took " << diff/CLOCKS_PER_SEC << " seconds to complete" << endl;
    cout << endl;
    //cout << "------------------------------" << endl;
    cout << endl;
    cout << "Direct Multiplication Result:" << endl;
    time_before = clock();
    Matrix C = direct_multiplication(A, k, n);
    time_after = clock();
    diff = time_after - time_before;
    print_matrix(C , n);
    //cout << endl;
    cout << "It took " << diff/CLOCKS_PER_SEC << " seconds to complete" << endl;
    cout << endl;
    cout << "------------------------------" << endl;
    //output_size(10);
    //output_exp(A, n);
    
    return 0;
}


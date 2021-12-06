//
//  CardGame.h
//  project4
//
//  Created by 仲崽 on 2021/10/1.
//  Copyright © 2021年 xx. All rights reserved.
//

#ifndef CardGame_h
#define CardGame_h

#include <vector>
#include <fstream>


using std::vector;
using namespace std;

class CardGame{
    
    int n;
    double *M;
    
public:
    CardGame(int n) {
        this->n = n;
        M = new double[(n+1)*(n+1)];
        for (int i=0; i<(n+1)*(n+1); i++)
            //Initialize the matrix as -1.
            M[i]=-1;
    }
    
    ~CardGame() {
        delete[] M;
    }
    
    
    double getValue(int r, int b){
        
        //for (int j=0; j<=n; j++)
        //    M[0*(n+1)+j] = j;
        //for (int i=0; i<=n; i++)
        //    M[i*(n+1)+0] = 0;
        //for (int i=1; i<=n; i++){
        //    for (int j=1; j<=n; j++){
        //        double r=(double) i/(i+j);
        //        double b=(double) j/(i+j);
        //        M[i*(n+1)+j] = max(r*M[(i-1)*(n+1)+j]+b*M[i*(n+1)+(j-1)], (double) j-i);
        //    }
        //}
        
        // Set up the first row, where red-cards-left=0;
        // Obviously, the player will say 'Stop' when there are no red cards left.
        // And the value at this time should be equal to black-cards-left.
        if (0==r){
            M[0*(n+1)+b] = b;
        }
        // Set up the first col, where black-cards-left=0;
        // Of course, the player will not say 'Stop' until all cards are dealt now.
        // And the value at this time should be 0.
        else if (0==b){
            M[r*(n+1)+0] = 0;
        }
        // Then, fill the whole matrix using equation 1 in the instructions.
        else {
            double pr=(double) r/(r+b);
            double pb=(double) b/(r+b);
            // If the value of the previous step has not been calculated, do it.
            if (-1 == M[(r-1)*(n+1)+b])
                M[(r-1)*(n+1)+b] = getValue(r-1,b);
            if (-1 == M[r*(n+1)+(b-1)])
                M[r*(n+1)+(b-1)] = getValue(r,b-1);
            M[r*(n+1)+b] = max(pr*M[(r-1)*(n+1)+b]+pb*M[r*(n+1)+(b-1)], (double) b-r);
        }
        
        //cout << endl;
        //for (int i=0; i<=n; i++){
            //for (int j=0; j<=n; j++){
                //cout << M[i*(n+1)+j] << "   ";
            //}
            //cout << endl;
        //}
        
        return M[r*(n+1)+b];
    }
    
};


#endif /* CardGame_h */

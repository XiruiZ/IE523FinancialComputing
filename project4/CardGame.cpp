//
//  main.cpp
//  project4
//
//  Created by 仲崽 on 2021/9/30.
//  Copyright © 2021年 xx. All rights reserved.
//

#include <iostream>
#include "CardGame.h"



int main(int argc, const char * argv[]) {
    // insert code here...
    cout << "The size of deck is:";
    int n;
    cin >> n;
    
    CardGame card(n/2);
    
    cout << "Total Number of Cards = " << n << endl;
    cout << "Value of Game = " << card.getValue(n/2, n/2) << endl;
    cout << endl;
    return 0;
}

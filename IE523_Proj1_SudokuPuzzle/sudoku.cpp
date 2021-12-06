// Soduku Solver using Brute-Force Search implemted using 
// recursion.
// Written for IE523: Financial Computation by Prof. Sreenivas
// and GE523: Discrete Event Systems
//
#include <iostream>
#include "sudoku.h"

//int main()
//{
  //  Sudoku x;
  //x.set_puzzle(1);
  // for (int i=0; i<9; i++){
  //     x.set_puzzle(i,0,i+1);
  // }
  //  int n=0;
  //  for (int i=0; i<3; i++){
  //      for (int j=0; j<3; j++){
  //          x.set_puzzle(i,j,n+1);
  //          n++;
  //      }
//   }
    
//    x.print_puzzle();
//    cout << x.block_valid(3,1) << endl;
//    return 0;
//}


int main (int argc, char * const argv[])
{
	Sudoku x;
	x.read_puzzle(argc, argv);
	x.print_puzzle();
	//x.Solve(0,0);
    x.alternate_Solve(0, 0);
    cout << endl;
    cout << "There are " << x.get_solu_num() << " solution(s) in total." << endl;

    return 0;
}

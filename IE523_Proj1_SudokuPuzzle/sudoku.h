/*
 *  sudoku.h
 *  Sudoku
 *  Created by Prof. Ramavarapu Sreenivas 
 *  Inspired by: http://web.eecs.utk.edu/courses/spring2012/cs140/Notes/Sudoku/index.html
*/

#ifndef sudoku
#define sudoku

#include <vector>
#include <fstream>
#include <algorithm>

using std::vector;
using namespace std;



class Sudoku 
{ 
	// Private
	int puzzle[9][9];
    int solu_num=0;

//public:
	// Private member function that checks if the named row is valid
	bool row_valid(int row)
	{
        for (int i = 0; i < 9; i++) {
            if (0 == puzzle[row][i]){ // No worries about 0
                continue;
            }
            for (int j = i+1; j < 9; j++){
                if (puzzle[row][i] == puzzle[row][j]){
                    return false; // Return false is there is duplication
                }
            }
        }

        return true; // Return true is there is no duplication and the row is vaild
	}
	
	// Private member function that checks if the named column is valid
	bool col_valid(int col)
	{
        for (int i = 0; i < 8; i++) {
            if (0 == puzzle[i][col]){ // No need to check for 0
                continue;
            }
            for (int j = i+1; j < 9; j++){
                if (puzzle[i][col] == puzzle[j][col]){
                    return false; // Return false is there is duplication
                }
            }
        }
        
        return true; // Return true is there is no duplication and the column is vaild
	}
	
	// Private member function that checks if the named 3x3 block is valid
	bool block_valid(int row, int col)
	{
        // First define which block it belongs to
        int rA, rB, cA, cB;
        
        if (0 == row%3){ // That is, row==0, 3, 6
            rA = row;
            rB = row+3;
        }
        else if (1 == row%3){ // That is, row==1, 4, 7
            rA = row-1;
            rB = row+2;
        }
        else{ // That is, row==2, 5, 8
            rA = row-2;
            rB = row+1;
        }
        // Do the silimar thing to col
        if (0 == col%3){ // That is, col==0, 3, 6
            cA = col;
            cB = col+3;
        }
        else if (1 == col%3){ // That is, col==1, 4, 7
            cA = col-1;
            cB = col+2;
        }
        else{ // That is, col==2, 5, 8
            cA = col-2;
            cB = col+1;
        }
        
        // Then check if there is duplication
        int ele = puzzle[row][col];
        for (int i=rA; i<rB; i++){
            for (int j=cA; j<cB; j++){
                if (i==row && j==col){ // That is the element itself
                    continue;
                }
                if (puzzle[i][j]==ele){
                    return false; // If there is duplication, return False
                }
            }
        }
        
        return true; // If there is NO duplication, return True
	}
	
public:
    
	// Public member function that reads the incomplete puzzle
	// we are not doing any checks on the input puzzle -- that is,
	// we are assuming they are indeed valid
    
    // For checking
    //void set_puzzle(int n){
    //    for (int i=0; i<9; i++){
    //        for (int j=0; j<9; j++){
    //            puzzle[i][j] = n;
    //        }
    //    }
    //}
    //void set_puzzle(int row, int col, int n){
    //    puzzle[row][col] = n;
    //}
    
    
    
	void read_puzzle(int argc, char * const argv[])
	{
        double value_read_from_file;
        vector <double> P;
        
        cout << "Input Puzzle is:";
        string file = "";
        cin >> file;
        
        ifstream input_file(file);
        if (input_file.is_open()){
            while(input_file >> value_read_from_file){
                P.push_back(value_read_from_file);
            }
        }
        else{
            cout << "Input file does not exist" << endl;
        }
        
        for (int i = 0; i < P.size(); i++){
            int r=i/9;
            int c=i%9;
            puzzle[r][c] = P[i];
        }
	}
	
	// Public member function that prints the puzzle when called
	void print_puzzle()
	{
		std::cout << std::endl << "Board Position" << std::endl;
		for (int i = 0; i < 9; i++)
		{
			for (int j = 0; j < 9; j++)
			{
				// check if we have a legitimate integer between 1 and 9
				if ((puzzle[i][j] >= 1) && (puzzle[i][j] <= 9))
				{
					// printing initial value of the puzzle with some formatting
					std::cout << puzzle[i][j] << " ";
				}
				else {
					// printing initial value of the puzzle with some formatting
					std::cout << "X ";
				}
			}
			std::cout << std::endl;
		}
	}
	
    
    
	// Public member function that (recursively) implements the brute-force 
	// search for possible solutions to the incomplete Sudoku puzzle
	bool Solve(int row, int col)
	{
		// this part of the code identifies the row and col number of the 
		// first incomplete (i.e. 0) entry in the puzzle.  If the puzzle has
		// no zeros, the variable row will be 9 => the puzzle is done, as 
		// each entry is row-, col- and block-valid...
		
		// use the pseudo code of figure 3 of the description
        
        // Find the next 0 grid
        
        int r=row;
        int c=col;
        
        while (true){
            
            if (0==(puzzle[r][c])){ // Found it
                break; // Get out of the while loop
            }
            // If not found yet, checl the next grid
            if (c<8){
                c++;
            }
            else{
                r++;
                c=0;
            }
            
            if (r>8){ // No more grids are 0, solved!
                print_puzzle();
                return true;
            }
            
        }
        
        // Try from 1-9
        for (int n=1; n<10; n++){
            puzzle[r][c] = n;
            if (row_valid(r) && col_valid(c) && block_valid(r,c)){
                if (Solve(r,c)){
                    return true;
                }
            }
        }
        
        /* If we got here then all assignments made to puzzle[i][j] are invalid. So, we reset its value and return false*/
        puzzle[r][c] = 0;
        return false;
	}
    
    
    // Find all solutions
    bool alternate_Solve(int row, int col){
        int r=row;
        int c=col;
        //int sn = solu_num;
        
        while (true){
            
            if (0==(puzzle[r][c])){ // Found it
                break; // Get out of the while loop
            }
            // If not found yet, checl the next grid
            if (c<8){
                c++;
            }
            else{
                r++;
                c=0;
            }
            
            if (r>8){ // No more grids are 0, solved!
                print_puzzle(); // print the current solution
                solu_num ++;
                return true;
            }
            
        }
        
        // Try from 1-9
        bool found_solu = false;
        for (int n=1; n<10; n++){
            puzzle[r][c] = n;
            if (row_valid(r) && col_valid(c) && block_valid(r,c)){
                if (alternate_Solve(r,c)){
                    found_solu = true;
                    // That is, there is at least 1 solution if fond_solu==true
                }
            }
        }
        
        puzzle[r][c] = 0; // Reset puzzle[i][j] to 0 for further solutions searching
        
        return found_solu; // Return if at least 1 solution was found
        
    }
    
    int get_solu_num(){ // Show solution number
        return solu_num;
    }
    
};



#endif

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"tblib.h"



//-----------------------------------------------------------------------------------------------
// threeBody.c takes optimization requirements and then determines the optimal solution to a
// three body problem and writes that solution to a file in a directory titled 'Results'. This
// directory must exist for threeBody.c to run correctly.
//
// Inputs:
//	argv[1] - The objective to be optimized. This is 1 for objective 1 (minimum velocity)
//		  and 2 for objective 2 (minimum flight time)
// 	argv[2] - Moon-s/c clearance distance in m. Specifies how close the s/c is allowed to 
//		  to get to the Moon/
//	argv[3] - Search accuracy. Defines on how fine the search grid should ultimately become
//
// Output:
//	Prints the optimal solution to the specified problem to the command line	
//
//
// Note: This program is subject to find only local solution, and has no way of knowing whether
//	 or not the solution is a global solution. 
//----------------------------------------------------------------------------------------------
int main (int argc, char *argv[]){

double pi = 4 * atan(1.0);
int obj;
double cle, acc;
obj = atoi(argv[1]);
cle = atof(argv[2]);
acc = atof(argv[3]);


double f = pow(10, 7); // Initial relative Error
int iter = 0; // Initializing iterator
double guess[2] = { 50, pi };

while(fabs(f) > acc){ // while the search grid is less fine than acc specifies.
	printf("\rBeginning Grid Search iteration number %d", iter);
	fflush(stdout);
	double minVal = guessOpt(guess, iter, obj, cle, acc, &f); // Find minimum value in current grid	
	iter++; // Iterate so that guessOpt knows to create a smaller window and finer grid. 
}	

printf("Optimal Guess: %.2f m/s, %.2f deg\nFinal Accuracy: %.3f\n\n", guess[0], guess[1]*180/pi, fabs(f));

return 0;
}





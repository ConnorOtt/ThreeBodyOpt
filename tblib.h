// Header File for Assignment 6 - Three Body Problem


#include<stdio.h>
#include<stdlib.h>
#include<math.h>


//-----------------------------------------------------------------------------------------------
// yPrimeFunc Calculates rates of change at a given state and populates yPrime[] with them
//
// Inputs:
//	double stateArr[] - Array containing the current state of the s/c and Moon.
//	double yPrimep[]  - Array which will be populated with rates of change and accessed in
//			    rungeKutta. 
// Output:
//	double *yPrime 	  - Pointer to yPrime for access outside function. 
//
//----------------------------------------------------------------------------------------------
double* yPrimeFunc(double stateArr[],  double yPrime[]) {
	
double Xs, Ys, Xm, Ym, Vxs, Vys, Vxm, Vym;
Xs = stateArr[0];
Ys = stateArr[1];
Xm = stateArr[2];
Ym = stateArr[3];
Vxs = stateArr[4];
Vys = stateArr[5];
Vxm = stateArr[6];
Vym = stateArr[7];

double G=6.67408*pow(10,-11);
double mM=7.34767309e22;   //(kg) moon mass
double mE=5.97219e24;      //(kg) earth mass
double mS=28833;           //(kg) spacecraft mass

// distance between s/c and earth
double dSE = sqrt(pow(Xs,2) + pow(Ys,2));
// distance between moon and earth
double dME = sqrt(pow(Xm,2) + pow(Ym,2));

// distance between s/c and moon
double dSM = sqrt(pow((Xs - Xm),2) + pow((Ys - Ym),2));

// Compute accelarations based on current locations
double AseX = -G*mE*(Xs)/(pow(dSE,3)); // On S/C by Eath
double AseY = -G*mE*(Ys)/pow(dSE,3); // on S/C by Earth
	
double AmeX = -G*mE*(Xm)/pow(dME,3); // on Moon by Earth
double AmeY = -G*mE*(Ym)/pow(dME,3); // on Moon Eath

double AsmX = G*mM*(Xm - Xs)/pow(dSM,3); // on S/C by Moon
double AsmY = G*mM*(Ym - Ys)/pow(dSM,3); // on S/C by Moon
double AmsX = -AsmX*mS/mM; // on Moon by S/C
double AmsY = -AsmY*mS/mM; // on Moon by S/C

yPrime[0] = Vxs;	
yPrime[1] = Vys;
yPrime[2] = Vxm;
yPrime[3] = Vym;
                           
yPrime[4] = AseX + AsmX;
yPrime[5] = AseY + AsmY;
yPrime[6] = AmeX + AmsX;
yPrime[7] = AmeY + AmsY;


return yPrime; 
}

//-----------------------------------------------------------------------------------------------
// integrateCheck checks the integration termination conditions
//
//
// Inputs: 
//	double Xs, Ys, Xm, Ym - Current X and Y positions of s/c and Moon.
//	double clearance      - Distance in m defining how close the s/c is allowed to get to
//				the Moon 
// Output:
//	int stopFlag 	      - Flag used by rungeKutta which halts integration if equal to 1.
//
//----------------------------------------------------------------------------------------------
int integrateCheck(double Xs, double Ys, double Xm, double Ym, double clearance){

int stopFlag = 0;
double radEarth = 6371000; 	// [m] - Radius of Earth
double radMoon = 1737000; 	// [m] - Radius of Moon
double d_EM = 384403000;	  	      // [m] - initial distance from Earth to Moon;
double scEDist = sqrt(pow(Xs, 2) + pow(Ys, 2)); // [m] - Current distance to Earth from spacecraft
double scMDist = sqrt(pow((Xm-Xs), 2) + pow((Ym-Ys), 2)); // current distance to Moon from spacecraft


if (scEDist <= radEarth){
	stopFlag = 1;
}
else if (scEDist >= 1.5*d_EM){
	stopFlag = 1;
}
else if (scMDist <= radMoon + clearance){
	stopFlag = 1;
} 


return stopFlag;
}

//-----------------------------------------------------------------------------------------------
// rungeKutta integrates spacecraft location using Runge Kutta 1st order method
//
// Inputs:
//	int saveOn 	- Flag which tells rungeKutta to save each integration step to a file.
//	double *t	- Pointer to time - accessed from outside in order to omptimize for low
//			  est flight time.
//	double initStateArr[] - Array containing the inital state of spacecraft and moon
//	double clearance- Distance in m specifying how close the s/c is allowed to get to the 
//			  Moon. For use in checking stop condtions and writing output files.
//	double acc 	- Specifices grid accuracy, for use in writing output files.
// Output:
//	double *stateArr- Pointer to final state after integration. Technically, stateArr points
//			  to initStateArr, so initStateArr will be overwritten and must be re-
//			  initialized before attempting another integration.
//
//----------------------------------------------------------------------------------------------
double *rungeKutta(int saveOn, double *t, double initStateArr[],int obj,  double clearance, double acc){

int h = 10; // [s] - time step
*t = 0; // [s] - start time
int tf = 5*pow(10, 5); // [s] - end time

int arrSize  = 8; // State Array Size 
double tmpStateArr[8]; // tmpStateArr and yPrimeInp are temporary inputs which get overwritten and reused to get the rates of change at a given state. 
double yPrime1[8], yPrime2[8], yPrime3[8], yPrime4[8];
double *stateArr; // Set to the initial state, this will then be overwritten as integration iterates.
stateArr = initStateArr;

FILE *FID;
if (saveOn == 1){ // If the program should save to file
	char buffer[50];
	char *s = "./Results/Optimum_";
	snprintf(buffer, sizeof(char)*50, "%s%d_%.f_0p%.f.csv", s, obj, clearance, acc*10);
	FID = fopen(buffer, "w");
}

while (*t < tf){ // These are integration steps (y0, y1, y2, ..., yn)
	double *k1 = yPrimeFunc(stateArr, yPrime1);	
	
	// just in case we really really still want to try 4th order RK
	/*
	for (int i = 0; i < arrSize; i++){
		tmpStateArr[i] =  stateArr[1] + k1[i]*h/2; // Stepping forward to take a derivative
	}
	
	double *k2 = yPrimeFunc(tmpStateArr, yPrime2);
	for (int i = 0; i < arrSize; i++){
		tmpStateArr[i] = stateArr[i] + k2[i]*h/2; // Stepping forward to take a derivative
	}
	double *k3 = yPrimeFunc(tmpStateArr, yPrime3);
	
	for (int i = 0; i < arrSize; i++){
		tmpStateArr[i] = stateArr[i] + k3[i]*h; // Stepping forward to take a derivative. 
	}
	double *k4 = yPrimeFunc(tmpStateArr, yPrime4);
	
	*/

	// Implimenting 1ST ORDER RUNGE KUTTA (Euler's Method Essentially)
	for (int i = 0; i < arrSize; i++){
		stateArr[i] = stateArr[i] + h*(k1[i]);// + 2*k2[i] + 2*k3[i] + k4[i]); // Applying the derivative
	}


	if (saveOn == 1){ // Print current state to file if specified
		fprintf(FID, "%.2f, %f, %f, %f, %f, %f, %f, %f, %f\n", *t, stateArr[0], stateArr[1], stateArr[2], stateArr[3], stateArr[4], stateArr[5], stateArr[6], stateArr[7]);
	}

	// Checking stop conditions
	int stopFlag = integrateCheck(stateArr[0], stateArr[1], stateArr[2], stateArr[3], clearance);	
	if (stopFlag == 1){
		if (saveOn == 1){
			fclose(FID);
		}
		return stateArr;
	}

	*(t) +=h; // Iterate time step

}

if (saveOn == 1){ // Close file if it was opened
	fclose(FID);
}

return stateArr;
}


//-----------------------------------------------------------------------------------------------
// guessOpt optimizes a specified objective for a three body problem
//
//
// Inputs: 
//	double guess[]   - Array of size 2 which contains guess velocity and direction (angle)
//	double iter      - Iterator which defines how fine a search grid to create for search. The
//		           higher the iterator, the finer the search.
//	int obj  	 - Specifies either objective 1 or 2 which tells guessOpt what to optimize
//	double clearance - Distance in m defining how close the s/c can get to the Moon.
//	double acc 	 - Velocity in m/s which defines how fine a search should ultimately 
//			   performed by guessOpt.
// 	double *tol	 - Pointer to the current accuracy (or tolerance) which is updated with 
//			   each iteration of guessOpt.
// Output:
//	guess 		 - Pointer to current best guess, which is then fed back in as input for 
//			   next iteration of guessOpt
//
//----------------------------------------------------------------------------------------------
double guessOpt(double guess[], double iter, int obj, double clearance, double acc, double *tol){

// Initializing Constants and initial conditions
int numSearch = 10;
double pi = 4 * atan(1.0);
double G = 6.674*pow(10, -11); 		// [Nm^2/kg^2] - Gravitational Constant
double radEarth = 6371000; 	 	// [m] - Radius of Earth
double mE = 5.97219*pow(10, 24);	// [kg] - Mass of Earth 
double mM = 7.34767309*pow(10, 22); 	// [kg] - Mass of Moon

double theta0_s = 50; 		// [deg] - initial angle to spacecraft
double theta0_m = 42.5; 	// [deg] - initial angle to Moon

double d_ES = 338000000; 			      	// [m] - Initial distance from Earth to spacecraft
	double x_s0 = d_ES * cos(theta0_s * pi/180);  	// [m] - initial X position of spacecraft
	double y_s0 = d_ES * sin(theta0_s * pi/180);  	// [m] - initial y position of spacecraft

double d_EM = 384403000;			      	// [m] - initial distance from Earth to Moon;
	double x_m0 = d_EM * cos(theta0_m * pi/180);  	// [m] - initial x position of Moon
	double y_m0 = d_EM * sin(theta0_m * pi/180);  	// [m] - initial y position of Moon

double V_0s = 1000; 				      	// [m/s] Initial Velocity of spacecraft
	double Vx_s0 = V_0s * cos(theta0_s * pi/180); 	// [m/s] - initial x velocity of spacecraft
	double Vy_s0 = V_0s * sin(theta0_s * pi/180); 	// [m/s] - initial y velocity of spacecraft

double V_0m = sqrt(G*pow(mE, 2) / ((mE + mM)*d_EM));	// [m/s] - initial Velocity of Moon
	double Vx_m0 = V_0m * -sin(theta0_m * pi/180); 	// [m/s] - initial x velocity of Moon
	double Vy_m0 = V_0m * cos(theta0_m * pi/180); 	// [m/s] - initial y velocity of Moon

double initConds[8] = {x_s0, y_s0, x_m0, y_m0, Vx_s0, Vy_s0, Vx_m0, Vy_m0};

// Determining spacing & search ranges for iterating through velocities and angles
double velCenter = guess[0];
double velRange = 100*exp(-iter/2);
double vSpace = velRange / (numSearch - 1);
double vStrtGuess = velCenter - velRange/2;

double thetaCenter = guess[1];
double thetaRange = 2*pi*exp(-iter/2);
double tSpace = thetaRange / (numSearch - 1);
double tStrtGuess = thetaCenter - thetaRange/2;

*tol = velRange;

// Initialize opimization parameters
double vGuess, tGuess;
double minVal = 30 * 86400; // (More than 30 days to start - much larger than 100 m/s as well)
double tmpInitConds[8];
double t = 0; // Will be passed in as pointer to be accessed from inside and outside

// Iterating through grid looking for successful burns, successful burns are compared to previous
// best burns and replace the best burn if they are better than the previous best burn so that 
// the previous best burn is rewritten to be the current best burn which will then be compared
// to the next succesful burns in search for the next best burn. Why did you read all that?
for(int i = 0; i < numSearch; i++){
	vGuess = vStrtGuess + i * vSpace;
	for(int j = 0; j < numSearch; j++){
		tGuess = tStrtGuess + j * tSpace;
		// Increase in velocity in X and Y directions
		
		for (int i = 0; i<8; i++){
			tmpInitConds[i] = initConds[i];
		}
		tmpInitConds[4] = initConds[4] + vGuess * cos(tGuess); // Adding to spacecraft velocity
		tmpInitConds[5] = initConds[5] + vGuess * sin(tGuess); 
			
		// Result needs to be final position and flight time (called fTime in the loop below).
		double *finalState = rungeKutta(0, &t, tmpInitConds, obj, clearance, acc);	
		double finalDist = sqrt(pow(finalState[0], 2) + pow(finalState[1], 2));
		if (finalDist <= radEarth){ 	// We only care if it made it back to Earth
			if(obj == 1){ // If evaluating objective 1 (Min velocity)
				if (minVal > vGuess){	
					// Tolerance based on difference between previous and current best answer.
					minVal = vGuess;
	
					// Set new guess based on this best one
					guess[0] = vGuess;
					guess[1] = tGuess;
				}
			}
			else{ // Evaluating objective 2 (Min time)
				if (minVal > t){	
					//printf("\n\t\t Found new best time guess\n");	
					minVal = t;
					guess[0] = vGuess;
					guess[1] = tGuess;
				}
			} 
		}
			
	}
}
if (*tol < acc){ // write best case to file if accuracy has been reached
	for (int i = 0; i<8; i++){
		tmpInitConds[i] = initConds[i];
	}
	tmpInitConds[4] = initConds[4] + guess[0] * cos(guess[1]); // Adding to spacecraft velocity
	tmpInitConds[5] = initConds[5] + guess[0] * sin(guess[1]); 

	double* finalState = rungeKutta(1, &t, tmpInitConds, obj, clearance, acc);	
	printf("\nFinal Time = %.3f days \n", t/86400);
		
}


return minVal; 
}






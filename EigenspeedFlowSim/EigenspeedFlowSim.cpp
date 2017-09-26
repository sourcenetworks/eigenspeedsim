// EigenspeedFlowSim.cpp : Defines the entry point for the console application.
//
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <time.h>
#include "flowsim.h"

using namespace std;



int main(int argc, char *argv[])
{
	string temp;
	flowsim mySim;

	mySim.parseCommandLine(argc, argv);
	mySim.openBWCap();
	mySim.openInitMatrix();
	mySim.allocateArrays();
	mySim.simLoop();

	//TODO implement cleanup function
	mySim.cleanup();

	cin >> temp;

	return 0;


	
}
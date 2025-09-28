#include "Fit_Final.C"
#include "TString.h"

void AutoRun() {

	TString Output_Name = "Record";

	float MassBin[9] = { 7.5,17.5,27.5,37.5,47.5,57.5,67.5,77.5,1000 };

	for (int MassNumber = 1; MassNumber < 9; MassNumber++) { 
		Fit_Final(MassBin[MassNumber - 1], MassBin[MassNumber], MassNumber, Output_Name);
	}

}

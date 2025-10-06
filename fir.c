/*******************************************************************************
Vendor: Xilinx 
Associated Filename: fir.c
Purpose: Vivado HLS Tutorial Example 
Device: All 
Revision History: May 30, 2008 - initial release
                                                
*******************************************************************************
Copyright 2008 - 2012 Xilinx, Inc. All rights reserved. 

This file contains confidential and proprietary information of Xilinx, Inc. and 
is protected under U.S. and international copyright and other intellectual 
property laws.

DISCLAIMER
This disclaimer is not a license and does not grant any rights to the materials 
distributed herewith. Except as otherwise provided in a valid license issued to 
you by Xilinx, and to the maximum extent permitted by applicable law: 
(1) THESE MATERIALS ARE MADE AVAILABLE "AS IS" AND WITH ALL FAULTS, AND XILINX 
HEREBY DISCLAIMS ALL WARRANTIES AND CONDITIONS, EXPRESS, IMPLIED, OR STATUTORY, 
INCLUDING BUT NOT LIMITED TO WARRANTIES OF MERCHANTABILITY, NON-INFRINGEMENT, OR 
FITNESS FOR ANY PARTICULAR PURPOSE; and (2) Xilinx shall not be liable (whether 
in contract or tort, including negligence, or under any other theory of 
liability) for any loss or damage of any kind or nature related to, arising under 
or in connection with these materials, including for any direct, or any indirect, 
special, incidental, or consequential loss or damage (including loss of data, 
profits, goodwill, or any type of loss or damage suffered as a result of any 
action brought by a third party) even if such damage or loss was reasonably 
foreseeable or Xilinx had been advised of the possibility of the same.

CRITICAL APPLICATIONS
Xilinx products are not designed or intended to be fail-safe, or for use in any 
application requiring fail-safe performance, such as life-support or safety 
devices or systems, Class III medical devices, nuclear facilities, applications 
related to the deployment of airbags, or any other applications that could lead 
to death, personal injury, or severe property or environmental damage 
(individually and collectively, "Critical Applications"). Customer assumes the 
sole risk and liability of any use of Xilinx products in Critical Applications, 
subject only to applicable laws and regulations governing limitations on product 
liability. 

THIS COPYRIGHT NOTICE AND DISCLAIMER MUST BE RETAINED AS PART OF THIS FILE AT 
ALL TIMES.

*******************************************************************************/
#include "fir.h"
#include <stdio.h>

#define NULL_PTR NULL


//helper function for matrix multiplication
static float dotProduct(float* row, float* column, int length){
	float sum=0;
	for(int i=0; i<length; i++){
		sum += row[i] * column[i];
	}
	return sum;
}

//helper function to update weights
static void weightUpdate(float* w, float mu, float error, float* shift_reg, int length){
	for(int i=0; i<length; i++){
		w[i] += mu * error * shift_reg[i];
	}
}

//helper function for shift register
static void shift_insertion(float* shift_reg, float new_value, int length){
	for(int i=length-1; i>0; --i){
		shift_reg[i] = shift_reg[i-1];
	}
	shift_reg[0]= new_value;
}

int fir(float* y_in, float mu, float* ref, int nbTrain, float* output, int totalNumber){
	int result = -1;
	int const taps = 5;

	//Making sure the addresses passed are valid
	if(y_in != NULL_PTR && ref != NULL_PTR && output != NULL_PTR){
		float w[taps] = {0};
		float shift_reg[taps]= {0};
		float error[nbTrain];

		//least mean square algo
		for(int i=0; i<nbTrain; i++){
			shift_insertion(shift_reg, y_in[i], taps);
			output[i] = dotProduct(w,shift_reg,taps);
			error[i] = ref[i] - output[i];
			weightUpdate(w, mu, error[i], shift_reg, taps);
		}

		//equalization
		for(int j = nbTrain; j<totalNumber; j++){
			shift_insertion(shift_reg, y_in[j], taps);
			output[j] = dotProduct(w,shift_reg,taps);
			//error[j] = ref[j] - output[j];

		}

		result = 0; //success
	}

	else {
		// do nothing, it will return error
	}


	return result;

}

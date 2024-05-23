/*****************************

WRITTEN BY ISTVÁN MÁRTON

*****************************/

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define length 4096
#define NUM_OF_THREADS 16384
#define RANK_OF_NORM 20
#define NUM_OF_BITS 8 * sizeof(unsigned long long int)

__global__ void L1(int* d_mtx_to_vec, unsigned long long int steps, unsigned long long int steps_remainder, int *d_L1_vector, int *d_L1_strategy, int iLonger, int iShorter){ // This function calculates the L1 norm on the GPU.
	int i, l, logical;
	int temp[length], vect[NUM_OF_BITS], product, L1;
	unsigned long long int number, index, iMax, iMin, iNumofZeros, iNum_temp;

	index = blockIdx.x * blockDim.x + threadIdx.x; // Index of threads.

	iMax = (index + 1) *(steps) - 1; // This part calculates the minimal (iMin-th) and the maximal (iMax-th) word of the binary reflected Gray code for which the calculations must be performed by a given thread.
	iMin = index * (steps);
	if(index < (steps_remainder) ) iMax += index + 1;
	else iMax += (steps_remainder);
	if(index <= (steps_remainder)) iMin += index;
	else iMin += (steps_remainder);
	number = iMin;

	for(l=0; l < iLonger; l++) {temp[l] = d_mtx_to_vec[l];} // As the code can consider a row of the matrix with fixed sign, it considers the first row of the matrix with +1.
	product = 0;
	for(i = 1 ; iShorter > i; i++){
		iNum_temp = (unsigned long long int) 1 << i; // iNum_temp and iNumofZeros are coefficients to determine the number-th word of BRGC.
		iNumofZeros=(unsigned long long int) iNum_temp >> 1;		
		logical = ((number+ iNumofZeros)/iNum_temp) % 2; // logical can be 0 and 1. logical is the number-th word and i-th digit of the BRGC.
		vect[i] = (int) 2 * logical - 1; // vect is the possible strategy vector. It's elements consists of +1 and -1.
			if(vect[i] > 0){for(l=0; l < iLonger; l++){temp[l] += d_mtx_to_vec[i * iLonger + l]; }} // The code determines the vector-matrix multiplication belonging to the number-th word of the BRGC.
			else {for(l=0; l < iLonger; l++){temp[l] -= d_mtx_to_vec[i * iLonger + l]; }}				
	}
	for(l= 0; l < iLonger; l++) {product += abs(temp[l]);} // The code calculates the L1 value belonging to the number-th word of the BRGC.
	L1 = product; 
	for(l=1; l<iShorter; l++){d_L1_strategy[index * (iShorter - 1) + l - 1] = vect[l];} // The program stores the strategy vector belonging to the number-th BRGC word in the d_L1_strategy vector.

	for(number=iMin + 1; number <= iMax; number++){ //The code determines the BRGC words till number variable reaches iMax.
		product = 0;
		for(i = 1 ; iShorter > i; i++){
			iNum_temp = (unsigned long long int) 1 << i; // iNum_temp and iNumofZeros are coefficients to determine the j-th word of BRGC.
			iNumofZeros=(unsigned long long int) iNum_temp >> 1;
			if( ((number+ iNumofZeros) % iNum_temp) == 0 ) {vect[i]=-vect[i]; // The code calculates if there is a change in the i-th digit in the BRGC.
				if(vect[i] > 0){for(l=0; l < iLonger; l++){temp[l] += 2 * d_mtx_to_vec[i * iLonger + l]; }} // When the i-th digit is changed, the code changes the result of the vector-matrix multiplication. It only deals with the i-th row of the matrix.
				else {for(l=0; l < iLonger; l++){temp[l] -= 2 *d_mtx_to_vec[i * iLonger + l]; }}
			break;
		}
            		}
	for(l = 0; l < (iLonger ); l++) {product += abs(temp[l]);} // Calculates the number-th possible L1 value.
		if(product > L1) {L1 = product; // If the current possible L1 value is greater than the previous, it modifies both the value and the corresponding strategy vector as well.
			for(l=1; l<iShorter; l++){d_L1_strategy[index * (iShorter - 1) + l - 1] = vect[l];}
		}
    }
d_L1_vector[index] = L1; // Every thread writes the biggest found L1 value to the d_L1_vector. This vector will be copied to the host memory.
}

__global__ void L2(int* d_mtx_to_vec, unsigned long long int steps, unsigned long long int steps_remainder, int *d_L2_vector, int *d_L2_strategy, int iRows, int iCols){
	int i, l, logical;
	int temp_negative[length], temp_positive[length], vect[NUM_OF_BITS], product, L2;
	unsigned long long int number, index, iMax, iMin, iNumofZeros, iNum_temp;

	index = blockIdx.x * blockDim.x + threadIdx.x;

	iMax = (index + 1) *(steps) - 1;
	iMin = index * (steps);
	if(index < (steps_remainder) ) iMax += index + 1;
	else iMax += (steps_remainder);
	if(index <= (steps_remainder)) iMin += index;
	else iMin += (steps_remainder);
	number = iMin;

	for(l=0; l < iCols; l++) {temp_negative[l] = 0; temp_positive[l] = d_mtx_to_vec[l]; }
	product = 0;
	for(i = 1 ; iRows > i; i++){
		iNum_temp = (unsigned long long int) 1 << i;
		iNumofZeros=(unsigned long long int) iNum_temp >> 1;		
		logical = ((number+ iNumofZeros)/iNum_temp) % 2;
		vect[i] = (int) 2 * logical - 1;
			if(vect[i] > 0){for(l=0; l < iCols; l++){temp_positive[l] += d_mtx_to_vec[i * iCols + l]; }}
			else {for(l=0; l < iCols; l++){temp_negative[l] += d_mtx_to_vec[i * iCols + l]; }}				
	}

	for(l= 0; l < iCols; l++) {product += abs(temp_negative[l]) + abs(temp_positive[l]);}
	L2 = product;
	for(l=1; l<iRows; l++){d_L2_strategy[index * (iRows - 1) + l - 1] = vect[l];}

	for(number=iMin + 1; number <= iMax; number++){
	product = 0;
	for(i = 1 ; iRows > i; i++){
		iNum_temp = (unsigned long long int) 1 << i;
		iNumofZeros=(unsigned long long int) iNum_temp >> 1;
		if( ((number+ iNumofZeros) % iNum_temp) == 0 ) {vect[i]=-vect[i] ;				
			if(vect[i] > 0){for(l=0; l < iCols; l++){temp_positive[l] += d_mtx_to_vec[i * iCols + l]; temp_negative[l] -= d_mtx_to_vec[i * iCols + l]; }}
			else {for(l=0; l < iCols; l++){temp_positive[l] -= d_mtx_to_vec[i * iCols + l]; temp_negative[l] += d_mtx_to_vec[i * iCols + l]; }}
			break;
		}
	}
	
	for(l = 0; l < (iCols ); l++) {product += abs(temp_negative[l]) + abs(temp_positive[l]);}
	if(product > L2) {L2 = product;
		for(l=1; l<iRows; l++){d_L2_strategy[index * (iRows - 1) + l - 1] = vect[l];}
	}
    }
d_L2_vector[index] = L2;
}

__global__ void L3(int* d_mtx_to_vec, unsigned long long int steps, unsigned long long int steps_remainder, int *d_L3_vector, int *d_L3_strategy, int iRows, int iCols, unsigned long long int *d_iNumPower){
	int i, l, helper[6] = {0, 1, 2, 2, 1, 0};

	int temp_0[length], temp_1[length], temp_2[length], vect[NUM_OF_BITS + 1], product, L3 = 0, logical, temporary;
	unsigned long long int number, index, iMax, iMin, divide;

	index = blockIdx.x * blockDim.x + threadIdx.x;

	iMax = (index + 1) *(steps) - 1;
	iMin = index * (steps);
	if(index < (steps_remainder) ) iMax += index + 1;
	else iMax += (steps_remainder);
	if(index <= (steps_remainder)) iMin += index;
	else iMin += (steps_remainder);
	
	number = iMin;
	for(l=0; l < iCols; l++) {temp_0[l] = d_mtx_to_vec[(iRows-1) * iCols + l]; temp_1[l] = 0; temp_2[l] = 0;}
	product = 0;
	for(i = 0 ; (iRows - 1) > i; i++){
		logical = (number/d_iNumPower[i]) % 6; // Determines the ternary reflected Gray code (TRGC). d_iNumPower is a vector consisting of the power of 3. This vector was copied from the host to speed up the calculation of TRGC as the power of 3 does not need to be determined every time.
		vect[i] = helper[logical];
		switch(vect[i])
			{
			case 0:
			for(l=0; l < iCols; l++){temp_0[l] += d_mtx_to_vec[i * iCols + l]; }
			break;
			case 1:
			for(l=0; l < iCols; l++){temp_1[l] += d_mtx_to_vec[i * iCols + l]; }
			break;
			case 2:
			for(l=0; l < iCols; l++){temp_2[l] += d_mtx_to_vec[i * iCols + l]; }
			break;
			}				
	}
	
	for(l= 0; l < iCols; l++) {product += abs(temp_0[l]) + abs(temp_1[l]) + abs(temp_2[l]);}
	L3 = product;
	for(l=0; l<(iRows - 1); l++){d_L3_strategy[index * (iRows - 1) + l] = vect[l];}

	for(number=iMin + 1; number <= iMax; number++){
      	product = 0;
		for(i = 0 ; (iRows - 1) > i; i++){
			divide = number/d_iNumPower[i];
			logical = divide % 3; //printf("%d, ", logical);
			if(logical) {
				logical = divide % 6;
				temporary = helper[logical];
				if( (vect[i] == 0)  && (temporary == 1) ) {for(l=0; l < iCols; l++){temp_0[l] -= d_mtx_to_vec[i * iCols + l]; temp_1[l] += d_mtx_to_vec[i * iCols + l];}}
					else if((vect[i] == 1)  && (temporary == 2)) {for(l=0; l < iCols; l++){temp_1[l] -= d_mtx_to_vec[i * iCols + l]; temp_2[l] += d_mtx_to_vec[i * iCols + l];}}
					else if((vect[i] == 2)  && (temporary == 1)) {for(l=0; l < iCols; l++){temp_1[l] += d_mtx_to_vec[i * iCols + l]; temp_2[l] -= d_mtx_to_vec[i * iCols + l];}}
					else {for(l=0; l < iCols; l++){temp_0[l] += d_mtx_to_vec[i * iCols + l]; temp_1[l] -= d_mtx_to_vec[i * iCols + l];}}
					vect[i] = temporary;
					break;
			}				
		}
		
	for(l= 0; l < iCols; l++) {product += abs(temp_0[l]) + abs(temp_1[l]) + abs(temp_2[l]);}
		if(product > L3) {
			L3 = product;
			for(l=0; l<(iRows - 1); l++){d_L3_strategy[index * (iRows - 1) + l] = vect[l];}
		}
	}
d_L3_vector[index] = L3;
}

__global__ void Ln(int* d_mtx_to_vec, int* d_iHelper, unsigned long long int steps, unsigned long long int steps_remainder, int *d_Ln_vector, int *d_Ln_strategy, int iRows, int iCols, int n, unsigned long long int *d_iNumPower){
	int i, l;
	int temp[RANK_OF_NORM][length], vect[NUM_OF_BITS + 1], product, Ln = 0;
	unsigned long long int number, index, iMax, iMin, divide;

	int logical, temporary;
	index = blockIdx.x * blockDim.x + threadIdx.x;
	iMax = (index + 1) *(steps) - 1;
	iMin = index * (steps);
	if(index < (steps_remainder) ) iMax += index + 1;
	else iMax += (steps_remainder);
	if(index <= (steps_remainder)) iMin += index;
	else iMin += (steps_remainder);
	number = iMin;

	for(l=0; l < iCols; l++) {
		temp[0][l] = d_mtx_to_vec[(iRows-1) * iCols + l]; 
		for(i=1; i<n; i++){
			temp[i][l] = 0;
		}
	}
	
	product = 0;
	for(i = 0 ; (iRows - 1) > i; i++){
		logical = (number/d_iNumPower[i]) % (2*n); 
		vect[i] = d_iHelper[logical]; //d_iHelper is a vector consisting of the power of n. It helps to determine the words of the n-ary Gray code.
		for(l=0; l < iCols; l++){temp[vect[i]][l] += d_mtx_to_vec[i * iCols + l]; }				
	}

	for(l= 0; l < iCols; l++) {
		for(i=0; i < n; i++){
			product += abs(temp[i][l]);
		}
	}
	Ln = product;
	for(l=0; l<(iRows - 1); l++){d_Ln_strategy[index * (iRows - 1) + l] = vect[l];}

	for(number=iMin + 1; number <= iMax; number++){
		product = 0;
		for(i = 0 ; (iRows - 1) > i; i++){
			divide = number/d_iNumPower[i];
			logical = divide % n; //printf("%d, ", logical);
			if(logical) {
				logical = divide % (2*n);
				temporary = d_iHelper[logical];
				for(l=0; l < iCols; l++) {temp[vect[i]][l] -= d_mtx_to_vec[i * iCols + l]; temp[temporary][l] += d_mtx_to_vec[i * iCols + l];}
				vect[i] = temporary;
				break;
			}			
		}
		for(l= 0; l < iCols; l++) {
			for(i=0; i < n; i++){
				product += abs(temp[i][l]);
			}
		}
		if(product > Ln) {
			Ln = product;
			for(l=0; l<(iRows - 1); l++){d_Ln_strategy[index * (iRows - 1) + l] = vect[l];}
		}
	}
d_Ln_vector[index] = Ln;
}

int** mtx_read(int *iRows, int *iCols, char* fileName){
	int i = 0,j = 0, k = 0;
	int *row, **mtx, value;
	
	mtx = NULL;
	row = NULL;
	
	char g, cNum[256];
	
	FILE *fp;
	fp = fopen(fileName,"r");
	
	do{
		g = fgetc(fp);	
		if((((g - '0') < 10) && ((g - '0') >= 0)) || (g == 'e') || ( g == 'E') || (g == '.') || (g == '+') || (g == '-')) {cNum[i] = g; i++;}
		else {
			cNum[i] = '\0'; 
			if(cNum[0] != '\0') {sscanf(cNum, "%d", &value); j++; i = 0;  row = (int*) realloc(row, j * sizeof(int)); row[j-1] = value;}
			if( ((g == '\n') || (g == EOF)) && (j > 0)){*iCols = j; j = 0; k++; mtx = (int**) realloc(mtx, k * sizeof(int*)); mtx[k-1] = row; row = NULL;}
		}
		
	}while(!feof(fp));
	*iRows = k;
printf("rows: %d, cols: %d\n",*iRows, *iCols); 
	fclose(fp);
return mtx;
}

void fileN(char *fileName, char** argv, int *argc){
	int r;
	FILE *fp;
	fp = NULL;
	if((*argc) < 2) {
		while(fp == NULL){
			printf("Please give me a filename that exist within this directory: ");
			r = scanf("%s",fileName);
			if(r == 0) {printf("Something went wrong when a filename was typed!\n");}
			fp = fopen(fileName, "r");
		}
	}
	else if((*argc) < 3){
		sprintf(fileName,"%s", argv[1]);
		fp = fopen(fileName, "r");
		while(fp == NULL){
			printf("Please give me a filename that exist within this directory: ");
			r = scanf("%s",fileName);
			if(r == 0) {printf("Something went wrong when a filename was typed!\n");}
			fp = fopen(fileName, "r");
		}
	}
	else if((*argc) < 4){
		sprintf(fileName,"%s", argv[1]);
		fp = fopen(fileName, "r");
		if(fp == NULL) {
			sprintf(fileName,"%s", argv[2]);
			fp = fopen(fileName, "r");
			while(fp == NULL){
				printf("Please give me a filename that exist within this directory: ");
				r = scanf("%s",fileName);
				if(r == 0) {printf("Something went wrong when a filename was typed!\n");}
				fp = fopen(fileName, "r");
			}
		}
	}
	else {
		printf("Too many parameters!\n");
		exit(-1);
	}
	fclose(fp);
}

void nNumber(int* p, char** argv, int *argc){
	int sd, n = 0;
	char line[256];
	
	if((*argc) < 2) {
		printf("Please give an integer bigger than zero: ");
		fgets(line, sizeof(line), stdin);
		sd = sscanf(line, "%d", &n);
		while((sd == 0) || (n < 1)){
			if(sd == 0) {printf("The entry must be an integer! Please give an integer bigger than zero: ");}
			if(n < 1  && sd != 0) {printf("Please give an integer bigger than zero: ");}
			fgets(line, sizeof(line), stdin);
			sd = sscanf(line, "%d", &n);
			//printf("sd value: %lu\n", sizeof(line));
		}
	}
	else if((*argc) < 3){
		sd = sscanf(argv[1], "%d", &n);
		while((sd == 0) || (n < 1)){
			if(sd == 0) {printf("Please give an integer bigger than zero: ");}
			if(n < 1  && sd != 0) {printf("Please give an integer bigger than zero: ");}
			fgets(line, sizeof(line), stdin);
			sd = sscanf(line, "%d", &n);
			//printf("sd value: %lu\n", sizeof(line));
		}
	}
	else if((*argc) < 4){
		sd = sscanf(argv[1], "%d", &n);
		if((sd == 0) || (n < 1)) {sd = sscanf(argv[2], "%d", &n);}
		while((sd == 0) || (n < 1)){
			if(sd == 0) {printf("Please give an integer bigger than zero: ");}
			if(n < 1  && sd != 0) {printf("Please give an integer bigger than zero: ");}
			fgets(line, sizeof(line), stdin);
			sd = sscanf(line, "%d", &n);
			//printf("sd value: %lu\n", sizeof(line));
		}
	}
	else {
		printf("Too many parameters!\n");
		exit(-1);
	}

	if(n > RANK_OF_NORM) {printf("The order of the L norm is too big. Please increase the RANK_OF_NORM variable in the code to %d and compile and run it again!\n", n); exit(-1);}
	*p = n;
}

void calc_Lnorm(int* n, int* iRows, int* iCols, int** mtx){
	FILE *fp;
	char fileOutput[1024]; // The variable 'fileOutput' is the name of the file, the 
	int i, j, iMax, *mtx_to_vec, *d_mtx_to_vec, maxRows, *Ln_vector, *d_Ln_vector, Ln_max, *Ln_strategy, *d_Ln_strategy, num_ofBlock, num_ofThread; //i and j are the indices of the input matrix; iMax is the variable belonging to the strategy vector, the index of the strategy vector; mtx_to_vec: the input matrix is converted to a vector in the host; d_mtx_to_vec: the converted matrix in the device; maxRows: them maximal number of rows (in case of L1, the maximal number of rows or columns) of the matrix the program can deal with, this number is determined by the order of the L norm that should be calculated; Ln_vector and d_Ln_vector are the two vectors containing the possible L norms belonging to a given thread in the host and device respectively; Ln_max: The maximal possible value of the L norm in the host; Ln_strategy and d_Ln_strategy are the vector containing all the possible strategy vectors belonging to different threads; num_ofBlock is the number of blocks the program uses; num_ofThread is the number of thread in a block
	unsigned long long int steps, steps_remainder, Inner_num, copyNum; //Inner_num: the number of possible L value the device should calculate; copyNum: the number of possible L values that will be copied from the device to the host 
	cudaDeviceProp devProp; // devProp contains the number of cores a warp contain.
	cudaGetDeviceProperties(&devProp, 0);
	mtx_to_vec = (int*)calloc(*iRows * *iCols, sizeof(int)); // allocating memory for the mtx_to_vec variable

	if(NUM_OF_THREADS < 1) {printf("The  NUM_OF_THREADS variable must be greater than 1. Modify it, and compile again!\n"); exit(-1);}
	if(*n == 1){ // if the order of the L norm is 1 then this part of the code will be executed.
		int iShorter, iLonger; //iShorter is the number of rows or columns, whichever is less; iLonger is the number of rows or columns, whichever is bigger
		if( *iRows > *iCols ){ //In this if else sequence the code transposes the matrix if necessary and transform a matrix into a vector
			for(j = 0; j < *iCols; j++){
				for(i = 0; i < *iRows; i++){
					mtx_to_vec[j * *iRows + i] = mtx[i][j];
				}
			}
		}
		else{
			for(i = 0; i < *iRows; i++){
				for(j = 0; j < *iCols; j++){
					mtx_to_vec[i * *iCols + j] = mtx[i][j];
				}
			}
		}
		if(*iRows < *iCols) {iShorter = *iRows; iLonger = *iCols;}
		else {iShorter = *iCols; iLonger = *iRows;}
		if(iShorter > (NUM_OF_BITS)) {printf("Matrix is too big. The number of rows or columns can not be more than %lu.\n", NUM_OF_BITS); exit(-1);}
		if(iLonger > length) {printf("Matrix is too big. The length variable %d should be bigger or equal than %d.\n", length, iLonger); exit(-1);}
		cudaMalloc((void**)&d_mtx_to_vec, iShorter * iLonger * sizeof(int)); // Allocating memory for the matrix in the device.
		Inner_num = (unsigned long long int) 1 << (iShorter - 1); //The number of possible L values that the device should calculate is 2^(iShorter-1)
		copyNum = NUM_OF_THREADS > Inner_num ? Inner_num : NUM_OF_THREADS; // The possible number of L norms can not be more than the number of threads
		num_ofThread = copyNum < devProp.warpSize ? copyNum : devProp.warpSize; // Number of threads in a block can not be bigger than the number of warps.
		num_ofBlock = copyNum/num_ofThread; copyNum = num_ofBlock * num_ofThread; //The number of blocks the code uses.
		steps = Inner_num/copyNum; steps_remainder = Inner_num % copyNum;
		Ln_vector = (int*) malloc(copyNum * sizeof(int)); // The code allocates memory in the host for the possible L norms.
		Ln_strategy = (int*) malloc(copyNum * (iShorter - 1) * sizeof(int)); // The code allocates memory for the possible strategies belonging to L norms in the host.
		cudaMalloc((void**)&d_Ln_vector, copyNum * sizeof(int)); // The code allocates memory in the device for the possible L norms.
		cudaMalloc((void**)&d_Ln_strategy, copyNum * (iShorter - 1) * sizeof(int)); // The code allocates memory for the possible strategies belonging to L norms in the device.
		cudaMemcpy(d_mtx_to_vec, mtx_to_vec, iShorter * iLonger * sizeof(int), cudaMemcpyHostToDevice); // The matrix is copied from RAM to GPU memory.
		printf("num_ofBlock: %d, num_ofThread: %d\n",num_ofBlock,num_ofThread);
		L1<<<num_ofBlock,num_ofThread>>>(d_mtx_to_vec, steps, steps_remainder, d_Ln_vector, d_Ln_strategy, iLonger, iShorter); // The calculation of the L1 norm with GPU.
		cudaMemcpy(Ln_vector, d_Ln_vector, copyNum * sizeof(int), cudaMemcpyDeviceToHost); // Copy the possible L norm values from device to host.
		cudaMemcpy(Ln_strategy, d_Ln_strategy, copyNum * (iShorter - 1) * sizeof(int), cudaMemcpyDeviceToHost); // Copy the possible strategies belonging to L norm values from device to host.
		Ln_max = Ln_vector[0]; iMax = 0; // Determining the maximal element of the Ln_vector which is the L norm, and the index of the strategy vector as well.
		for(i = 1; i < copyNum; i++){ if(Ln_max < Ln_vector[i]) {Ln_max = Ln_vector[i]; iMax = i;}}

		fp = fopen("strategy_L1.txt", "w");// Print out the strategy vector to file.
		fprintf(fp,"1\n");
		for(i=0; i<(iShorter - 1); i++) {fprintf(fp, "%d\n", Ln_strategy[iMax * (iShorter - 1) + i]);}
		fclose(fp);
	}
	else if(*n == 2){ // if the order of the L norm is 2 then this part of the code will be executed.

			for(i = 0; i < *iRows; i++){
				for(j = 0; j < *iCols; j++){
					mtx_to_vec[i * *iCols + j] = mtx[i][j];
				}
			}

		if(*iRows > (NUM_OF_BITS)) {printf("Matrix is too big. The number of rows or columns can not be more than %lu.\n", NUM_OF_BITS); exit(-1);}
		if(*iCols > length) {printf("Matrix is too big. The length variable %d should be bigger or equal than %d.\n", length, *iCols); exit(-1);}
		cudaMalloc((void**)&d_mtx_to_vec, *iRows * *iCols * sizeof(int));
		Inner_num = (unsigned long long int) 1 << (*iRows - 1);
		copyNum = NUM_OF_THREADS > Inner_num ? Inner_num : NUM_OF_THREADS; // The possible number of L norms can not be more than the number of threads
		num_ofThread = copyNum < devProp.warpSize ? copyNum : devProp.warpSize; // Number of threads in a block can not be bigger than the number of warps.
		num_ofBlock = copyNum/num_ofThread; copyNum = num_ofBlock * num_ofThread; //The number of blocks the code uses.
		steps=Inner_num/copyNum; steps_remainder = Inner_num % copyNum;
		Ln_vector = (int*) malloc(copyNum * sizeof(int));
		Ln_strategy = (int*) malloc(copyNum * (*iRows - 1) * sizeof(int));
		cudaMalloc((void**)&d_Ln_vector, copyNum * sizeof(int));
		cudaMalloc((void**)&d_Ln_strategy, copyNum * (*iRows - 1) * sizeof(int));
		cudaMemcpy(d_mtx_to_vec, mtx_to_vec, *iRows * *iCols * sizeof(int), cudaMemcpyHostToDevice);
		printf("num_ofBlock: %d, num_ofThread: %d\n",num_ofBlock,num_ofThread);
		L2<<<num_ofBlock,num_ofThread>>>(d_mtx_to_vec, steps, steps_remainder, d_Ln_vector, d_Ln_strategy, *iRows, *iCols);
		cudaMemcpy(Ln_vector, d_Ln_vector, copyNum * sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(Ln_strategy, d_Ln_strategy, copyNum * (*iRows - 1) * sizeof(int), cudaMemcpyDeviceToHost);
		Ln_max = Ln_vector[0]; iMax = 0;
		for(i = 1; i < copyNum; i++){ if(Ln_max < Ln_vector[i]) {Ln_max = Ln_vector[i]; iMax = i;}}

		fp = fopen("strategy_L2.txt", "w");
		fprintf(fp,"1\n");
		for(i=0; i<(*iRows - 1); i++) {fprintf(fp, "%d\n", Ln_strategy[iMax * (*iRows - 1) + i]);}
		fclose(fp);
	}
	else if(*n == 3){ // if the order of the L norm is 3 then this part of the code will be executed.
		unsigned long long int *iNumPower, *d_iNumPower;

			for(i = 0; i < *iRows; i++){
				for(j = 0; j < *iCols; j++){
					mtx_to_vec[i * *iCols + j] = mtx[i][j];
				}
			}
		cudaMalloc((void**)&d_mtx_to_vec, *iRows * *iCols * sizeof(int));
		Inner_num = pow(3, *iRows - 1);
		copyNum = NUM_OF_THREADS > Inner_num ? Inner_num : NUM_OF_THREADS; // The possible number of L norms can not be more than the number of threads
		num_ofThread = copyNum < devProp.warpSize ? copyNum : devProp.warpSize; // Number of threads in a block can not be bigger than the number of warps.
		num_ofBlock = copyNum/num_ofThread; copyNum = num_ofBlock * num_ofThread; //The number of blocks the code uses.
printf("Inner_num: %llu, copyNum: %llu, num_ofThread: %d, num_ofBlock: %d", Inner_num, copyNum, num_ofThread, num_ofBlock);
		maxRows = (int) (floor (NUM_OF_BITS / log2(*n)) + 1);
		if( *iRows > maxRows) {printf("Matrix is too big. The number of rows can not be more than %d.\n", maxRows); exit(-1);}
		if(*iCols > length) {printf("Matrix is too big. The length variable %d should be bigger or equal than %d.\n", length, *iCols); exit(-1);}
		printf("NUM_OF_BITS: %lu, maxRows: %d\n", NUM_OF_BITS ,maxRows);
		steps=Inner_num/copyNum; steps_remainder = Inner_num % copyNum;
		Ln_vector = (int*) malloc(copyNum * sizeof(int));
		Ln_strategy = (int*) malloc(copyNum * (*iRows - 1) * sizeof(int));
		iNumPower = (unsigned long long int*) malloc(maxRows * sizeof(unsigned long long int));
		cudaMalloc((void**)&d_Ln_vector, copyNum * sizeof(int));
		cudaMalloc((void**)&d_Ln_strategy, copyNum * (*iRows - 1) * sizeof(int));
		cudaMalloc((void**)&d_iNumPower, (maxRows-1) * sizeof(unsigned long long int));
		iNumPower[0] = 1;
		for(i = 1; i < (maxRows-1); i++){iNumPower[i] = iNumPower[i-1] * 3; } // iNumPower is copied to the device memory to speed up the calculation of the ternary Gray code
		cudaMemcpy(d_iNumPower, iNumPower, (maxRows-1) * sizeof(unsigned long long int), cudaMemcpyHostToDevice);
		cudaMemcpy(d_mtx_to_vec, mtx_to_vec, *iRows * *iCols * sizeof(int), cudaMemcpyHostToDevice);
		printf("num_ofBlock: %d, num_ofThread: %d\n",num_ofBlock,num_ofThread);
		L3<<<num_ofBlock,num_ofThread>>>(d_mtx_to_vec, steps, steps_remainder, d_Ln_vector, d_Ln_strategy, *iRows, *iCols, d_iNumPower);
		cudaMemcpy(Ln_vector, d_Ln_vector, copyNum * sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(Ln_strategy, d_Ln_strategy, copyNum * (*iRows - 1) * sizeof(int), cudaMemcpyDeviceToHost);
		Ln_max = Ln_vector[0]; iMax = 0;
		for(i = 1; i < copyNum; i++){if(Ln_max < Ln_vector[i]) {Ln_max = Ln_vector[i]; iMax = i; }}

		FILE *fp;
		fp = fopen("strategy_L3.txt", "w");
		for(i=0; i<(*iRows - 1); i++) {fprintf(fp, "%d\n", Ln_strategy[iMax * (*iRows - 1) + i]);}
		fprintf(fp,"0\n");
		fclose(fp);

		free(iNumPower);
		cudaFree(d_iNumPower);
	}
	else{ // if the order of the L norm is bigger than 3, then this part of the code will be executed.
		int *iHelper, *d_iHelper;
		unsigned long long int *iNumPower, *d_iNumPower;

		for(i = 0; i < *iRows; i++){
			for(j = 0; j < *iCols; j++){
				mtx_to_vec[i * *iCols + j] = mtx[i][j];
			}
		}
		*n = ((*n < *iRows)? *n:*iRows);
		cudaMalloc((void**)&d_mtx_to_vec, *iRows * *iCols * sizeof(int));
		Inner_num = pow(*n, *iRows - 1);
		copyNum = NUM_OF_THREADS > Inner_num ? Inner_num : NUM_OF_THREADS; // The possible number of L norms can not be more than the number of threads
		num_ofThread = copyNum < devProp.warpSize ? copyNum : devProp.warpSize; // Number of threads in a block can not be bigger than the number of warps.
		num_ofBlock = copyNum/num_ofThread; copyNum = num_ofBlock * num_ofThread; //The number of blocks the code uses.

		maxRows = (int) (floor (NUM_OF_BITS / log2(*n)) + 1);
		if(*iRows > maxRows) {printf("Matrix is too big. The number of rows can not be more than %d.\n", maxRows); exit(-1);}
		if(*iCols > length) {printf("Matrix is too big. The length variable %d should be bigger or equal than %d.\n", length, *iCols); exit(-1);}
		iHelper = (int*) calloc(2 * *n, sizeof(int));
		for(i = 0; i < *n; i++){
			iHelper[i] = i;
			iHelper[2 * *n-i-1]=i;
		}
		printf("NUM_OF_BITS: %lu, maxRows: %d\n", NUM_OF_BITS ,maxRows);
		steps=Inner_num/copyNum; steps_remainder = Inner_num % copyNum;
		Ln_vector = (int*) malloc(copyNum * sizeof(int));
		Ln_strategy = (int*) malloc(copyNum * (*iRows - 1) * sizeof(int));
		iNumPower = (unsigned long long int*) malloc(maxRows * sizeof(unsigned long long int));
		cudaMalloc((void**)&d_Ln_vector, copyNum * sizeof(int));
		cudaMalloc((void**)&d_Ln_strategy, copyNum * (*iRows - 1) * sizeof(int));
		cudaMalloc((void**)&d_iNumPower, (maxRows-1) * sizeof(unsigned long long int));
		cudaMalloc((void**)&d_iHelper, 2 * *n * sizeof(int));
		iNumPower[0] = 1;
		for(i = 1; i < (maxRows-1); i++){iNumPower[i] = iNumPower[i-1] * *n; } //printf("iNumPower: %llu\n", iNumPower[maxRows -2]);
		cudaMemcpy(d_iNumPower, iNumPower, (maxRows-1) * sizeof(unsigned long long int), cudaMemcpyHostToDevice);
		cudaMemcpy(d_mtx_to_vec, mtx_to_vec, *iRows * *iCols * sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(d_iHelper, iHelper, 2 * *n * sizeof(int), cudaMemcpyHostToDevice);
		printf("num_ofBlock: %d, num_ofThread: %d\n",num_ofBlock,num_ofThread); 
		Ln<<<num_ofBlock,num_ofThread>>>(d_mtx_to_vec, d_iHelper, steps, steps_remainder, d_Ln_vector, d_Ln_strategy, *iRows, *iCols, *n, d_iNumPower);
		cudaMemcpy(Ln_vector, d_Ln_vector, copyNum * sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(Ln_strategy, d_Ln_strategy, copyNum * (*iRows - 1) * sizeof(int), cudaMemcpyDeviceToHost);
		Ln_max = Ln_vector[0]; iMax = 0;
		for(i = 1; i < copyNum; i++){if(Ln_max < Ln_vector[i]) {Ln_max = Ln_vector[i]; iMax = i; }}

		sprintf(fileOutput,"strategy_L%d.txt", *n);
		fp = fopen(fileOutput, "w");
		fprintf(fp, "#L%d is: %d\n", *n, Ln_max);
		for(i=0; i<(*iRows - 1); i++) {fprintf(fp, "%d\n", Ln_strategy[iMax * (*iRows - 1) + i]);}
		fprintf(fp,"0\n");
		fclose(fp);

		free(iNumPower);
		free(iHelper);
		cudaFree(d_iNumPower);
		cudaFree(d_iHelper);
	}
	
	printf("L%d is: %d\n", *n, Ln_max); // Write out the value of the L norm to the screen.
	
	free(Ln_vector); // Deallocates the vectors in the host memory.
	free(Ln_strategy);
	free(mtx_to_vec);
	
	cudaFree(d_Ln_vector); // Deallocates the vectors in the device memory.
	cudaFree(d_Ln_strategy);
	cudaFree(d_mtx_to_vec);
}

void mtx_free(int* iRows, int** mtx){
	int i;
	for(i = 0; i < *iRows; i++){
		free(mtx[i]);
	}
	free(mtx);
}

int main(int argc, char *argv[]){
	char fileName[1024]; //The 'fileName' variable contains the name of the file.
	int iRows, iCols, **mtx, n; // These variables are the number of rows and columns of the matrix, the matrix itself, and 'n' is the order of the L norm one wants to calculate.
	fileN(fileName, argv, &argc); // The 'fileN' function ensures that the filename containing the matrix exists within the working directory.
	nNumber(&n, argv, &argc); // The 'nNumber' function ensures that the order of the L norm is calculatable by the code.
	mtx = mtx_read(&iRows, &iCols, fileName); // The 'mtx_read' function reads the file containing the matrix.
	calc_Lnorm(&n, &iRows, &iCols, mtx); // The 'calc_Lnorm' calculates the n order L norm of the matrix.
	mtx_free(&iRows, mtx); //The 'mtx_free' function deallocates the memory for the pointer mtx.
	return 0;     
}

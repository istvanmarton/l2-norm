#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define length 4096
#define NUM_OF_THREADS 10752
#define NUM_OF_BITS 8 * sizeof(unsigned long long int)

__global__ void func(int* d_mtx_to_vec, unsigned long long int steps, unsigned long long int steps_remainder, int *d_L3_vector, int *d_L3_strategy, int iRows, int iCols, int n, unsigned long long int *d_iNumPower){
	int i, l;
	int temp_0[length], temp_1[length], temp_2[length], vect[NUM_OF_BITS + 1], product, L3 = 0;
	unsigned long long int number, index, iMax, iMin, iNum_temp;

	int logical;
	index = blockIdx.x * blockDim.x + threadIdx.x;

	iMax = (index + 1) *(steps) - 1;
	iMin = index * (steps);
	if(index < (steps_remainder) ) iMax += index + 1;
	else iMax += (steps_remainder);
	if(index <= (steps_remainder)) iMin += index;
	else iMin += (steps_remainder);
//	if(index == 0) printf("index: %llu, iMin: %llu, iMax: %llu\n", index, iMin, iMax);

	for(number=iMin; number <= iMax; number++){
		 for(l=0; l < iCols; l++) {temp_0[l] = d_mtx_to_vec[(iRows-1) * iCols + l]; temp_1[l] = 0; temp_2[l] = 0;}
		 product = 0;
			for(i = 0 ; (iRows - 1) > i; i++){
				iNum_temp = d_iNumPower[i];//pow(n,i+1);
				logical = (number/iNum_temp) % n; //printf("%d, ", logical);
				vect[i] = logical;
					switch(logical)
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
			//printf("\n");
			for(l= 0; l < iCols; l++) {product += abs(temp_0[l]) + abs(temp_1[l]) + abs(temp_2[l]);}
			if(product > L3) {
				L3 = product;
				for(l=0; l<(iRows - 1); l++){d_L3_strategy[index * (iRows - 1) + l] = vect[l];}
			}
	}
d_L3_vector[index] = L3;
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
	if((*argc) < 2) {
		do{
			printf("Please give me a filename: "); 
			r = scanf("%s",fileName);
		}while(r != 1);
	}
	else sprintf(fileName,"%s", argv[1]);

	fp = fopen(fileName, "r");
	if(fp == NULL) {
		do{
			printf("Please give me a filename that exist within this directory: ");
			r = scanf("%s",fileName);
			fp = fopen(fileName, "r");
		}while(fp == NULL);
	}
	fclose(fp);
}

int main(int argc, char *argv[]){
cudaDeviceProp devProp;
cudaGetDeviceProperties(&devProp, 0);
     char fileName[1024];
     fileN(fileName, argv, &argc);     
     int i, j, iMax, iRows, iCols, **mtx, *mtx_to_vec, *d_mtx_to_vec, n = 3, maxRows = (int) (floor (NUM_OF_BITS / log2(n)) + 1);
     printf("NUM_OF_BITS: %lu, maxRows: %d\n", NUM_OF_BITS ,maxRows);
     mtx = mtx_read(&iRows, &iCols, fileName);
     mtx_to_vec = (int*)calloc(iRows * iCols, sizeof(int));

		for(i = 0; i < iRows; i++){
			for(j = 0; j < iCols; j++){
				mtx_to_vec[i * iCols + j] = mtx[i][j];
			}
		}
	if( iRows > maxRows) {printf("Matrix is too big. The number of rows can not be more than %d.\n", maxRows); return 0;}
	if(iCols > length) {printf("Matrix is too big. The length variable %d should be bigger or equal than %d.\n", length, iCols); return 0;}
	cudaMalloc((void**)&d_mtx_to_vec, iRows * iCols * sizeof(int));
     	cudaMalloc((void**)&d_mtx_to_vec, iRows * iCols * sizeof(int));
     	unsigned long long int steps, steps_remainder, Inner_num = (unsigned long long int) pow(3, iRows - 1)/*1 << (iShorter - 1)*/, copyNum, *iNumPower, *d_iNumPower;
	copyNum = NUM_OF_THREADS > Inner_num ? Inner_num : NUM_OF_THREADS;
	int *L3_vector, *d_L3_vector, L3_max = 0, *L3_strategy, *d_L3_strategy, num_ofBlock = (int) ceil((float)copyNum/(float)devProp.warpSize), num_ofThread = copyNum < devProp.warpSize ? copyNum : devProp.warpSize;

	steps=Inner_num/copyNum; steps_remainder = Inner_num % copyNum;

	L3_vector = (int*) malloc(copyNum * sizeof(int));
	L3_strategy = (int*) malloc(copyNum * (iRows - 1) * sizeof(int));
	iNumPower = (unsigned long long int*) malloc(maxRows * sizeof(unsigned long long int));

	cudaMalloc((void**)&d_L3_vector, copyNum * sizeof(int));
	cudaMalloc((void**)&d_L3_strategy, copyNum * (iRows - 1) * sizeof(int));
	cudaMalloc((void**)&d_iNumPower, maxRows * sizeof(unsigned long long int));
     	for(i = 0; i < maxRows; i++){iNumPower[i] = pow(n, i);} //printf("%llu\n", iNumPower[maxRows-1]);
	cudaMemcpy(d_iNumPower, iNumPower, maxRows * sizeof(unsigned long long int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_mtx_to_vec, mtx_to_vec, iRows * iCols * sizeof(int), cudaMemcpyHostToDevice);
printf("num_ofBlock: %d, num_ofThread: %d\n",num_ofBlock,num_ofThread);
	 func<<<num_ofBlock,num_ofThread>>>(d_mtx_to_vec, steps, steps_remainder, d_L3_vector, d_L3_strategy, iRows, iCols, n, d_iNumPower);
	cudaMemcpy(L3_vector, d_L3_vector, copyNum * sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(L3_strategy, d_L3_strategy, copyNum * (iRows - 1) * sizeof(int), cudaMemcpyDeviceToHost);	 
	for(i = 0; i < copyNum; i++){
		if(L3_max < L3_vector[i]) {L3_max = L3_vector[i]; iMax = i;}
	}
FILE *fp;
fp = fopen("strategy_L3.txt", "w");
	for(i=0; i<(iRows - 1); i++) {fprintf(fp, "%d\n", L3_strategy[iMax * (iRows - 1) + i]);}
	fprintf(fp,"0\n");
fclose(fp);

	printf("L3 is: %d\n",L3_max);

	free(L3_vector);
	free(L3_strategy);
	free(mtx_to_vec);
	free(iNumPower);

	cudaFree(d_L3_vector);
	cudaFree(d_L3_strategy);
	cudaFree(d_mtx_to_vec);
	cudaFree(d_iNumPower);

     return 0;     
}

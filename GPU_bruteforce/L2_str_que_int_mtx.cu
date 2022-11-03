#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define length 128

__global__ void func(int* d_mtx_to_vec, unsigned long long int numThreads, unsigned long long int steps, unsigned long long int steps_remainder, int *d_L2_vector, int *d_L2_strategy, int iLonger, int iShorter){
	int i, l;
//printf("Itt meg jol vagyok.\n");

//	int *temp = new int[iLonger], product = 0 ,L1 = 0, *vect = new int[iShorter], elojel;
	int temp_L2_negative[length], temp_L2_positive[length], vect[length], product = 0, elojel, L2 = 0;
	unsigned long long int number, index, iMax, iMin, iNumofZeros, iNum_temp;
	int logical;
	index = blockIdx.x * blockDim.x + threadIdx.x;

//int iiiMax = 0;
	
	iMax = (index + 1) *(steps) - 1;
	iMin = index * (steps);
	if(index < (steps_remainder) ) iMax += index + 1;
	else iMax += (steps_remainder);
	if(index <= (steps_remainder)) iMin += index;
	else iMin += (steps_remainder);
//	printf("index: %llu, iMin: %llu, iMax: %llu\n", index, iMin, iMax);
		 number = iMin;
		 for(l=0; l < iLonger; l++) {temp_L2_negative[l] = 0; temp_L2_positive[l] = 0;}
		 product = 0;
			for(i = (iShorter - 1); i >=0 ; i--){
				iNum_temp = 1 << (i+1);
				iNumofZeros=iNum_temp>>1;
				logical = ((number+ iNumofZeros)/iNum_temp) % 2;
				logical = logical == 0 ? 0 : 1;
//				printf("\nhj: %u, j: %d\n", h[j], j); 
 //             printf("i: %d, j: %d\n", i,j);
//				printf("%d", logical == 0 ? 0 : 1);
				vect[i] = (int) 2 * logical - 1;
				for(l=0; l < iLonger; l++){
					if(vect[i] < 0) { temp_L2_negative[l] -= d_mtx_to_vec[i * iLonger + l];}
					else {		  temp_L2_positive[l] += d_mtx_to_vec[i * iLonger + l];}
				}
			}
			for(l= 0; l < iLonger; l++) {product += abs(temp_L2_positive[l]) + abs(temp_L2_negative[l]);}
			L2 = product;
			for(l=0; l<(iShorter); l++){d_L2_strategy[index * (iShorter) + l] = vect[l];}

     for(number=iMin + 1; number <= iMax; number++){
//		 printf("k: %llu\n",number);
		 product = 0; //szamlalo=0;

			for(i = (iShorter - 1); i >=0 ; i--){
				iNum_temp = 1 << (i+1);
				iNumofZeros=iNum_temp>>1;
				logical = ((number+ iNumofZeros)/iNum_temp) % 2;
				if(vect[i] != (2*logical - 1)) {elojel = (int) 2*logical - 1; 
					for(l=0; l < iLonger; l++){
						if(elojel < 0){
							temp_L2_positive[l] -= 2 * d_mtx_to_vec[i * iLonger + l];
							temp_L2_negative[l] += 2 * d_mtx_to_vec[i * iLonger + l];
						}
						else{
							temp_L2_positive[l] += 2 * d_mtx_to_vec[i * iLonger + l];
							temp_L2_negative[l] -= 2 * d_mtx_to_vec[i * iLonger + l];	
						}

					}
//				szamlalo++;
//				if(szamlalo > 1) printf("A szamlalo: %d\n", szamlalo);
				vect[i] = (int) 2 * logical - 1;
				break;
			}
								
//				printf("vector: %d\n", vect[i]);
//              printf("logical: %d, mask: %d\n", logical, mask << i);
            		}
	     for(l = 0; l < (iLonger ); l++) {product += abs(temp_L2_positive[l]) + abs(temp_L2_negative[l]);}
	     if(product > L2) {L2 = product;
		for(l=0; l<(iShorter); l++){d_L2_strategy[index * (iShorter) + l] = vect[l]; /*printf("%d\n", vect[k]);*/}
		}
//		 meret++;
//      printf("\n");
     }
d_L2_vector[index] = L2;
//printf("%max: %llu\n",iiiMax);
}

int** mtx_read(int *iRows, int *iCols, char* fileName){
printf("%s\n",fileName);
	int i = 0,j = 0, k = 0;
	int *row, **mtx, value;
	
	mtx = NULL; //(double**)malloc(sizeof(double*));
	row = NULL;
//	mtx[0] = row;
	
	char g, cNum[256];
	
	FILE *fp;
	fp = fopen(fileName,"r");
	
	do{
		
		g = fgetc(fp); 
		
//		
		if((((g - '0') < 10) && ((g - '0') >= 0)) || (g == 'e') || ( g == 'E') || (g == '.') || (g == '+') || (g == '-')) {cNum[i] = g;  /*printf("i: %d\tc: %c\n", i, cNum[i]);*/ i++;}
		else {
//			printf("%d\t%d\n",i,j);
			cNum[i] = '\0'; 
			if(cNum[0] != '\0') {i = 0; j++; sscanf(cNum, "%d", &value); /*printf("value: %d, col: %d, row: %d,\n", value, i,j);*/ row = (int*) realloc(row, j * sizeof(int)); row[j-1] = value; 
			if((g == '\n') || (g == EOF)){/*printf("j: %d\n",j);*/ *iCols = j; j = 0; k++; /*printf("k: %d\n",k);*/ mtx = (int**) realloc(mtx, k * sizeof(int*)); mtx[k-1] = row; /*printf("mtx: %d\t%d\t%d\n", mtx[k-1][0], mtx[k-1][1], mtx[k-1][2]);*/ /*free(row);*/ row = NULL;}
			}
		}
		
		
			
//		cNum[i] =  ? ;
		
	}while(!feof(fp));
//	printf("j: %d\n",j);
	*iRows = k;
printf("rows: %d, cols: %d\n",*iRows, *iCols); 
	fclose(fp);
//	free(row);
return mtx;
}

void fileN(char *fileName, char** argv, int *argc){
	if((*argc) < 2) {printf("Please give me a filename: "); scanf("%s",fileName);}
	else sprintf(fileName,"%s", argv[1]);//fileName = argv[1];
//	printf("fileName :%s", fileName);
}


int main(int argc, char *argv[]){
     char fileName[1024];
     fileN(fileName, argv, &argc);     
     printf("%s\n",fileName);
     int i, j, iMax, iRows, iCols, **mtx, *mtx_to_vec, *d_mtx_to_vec;
     mtx = mtx_read(&iRows, &iCols, fileName);
     mtx_to_vec = (int*)calloc(iRows * iCols, sizeof(int));

		for(i = 0; i < iRows; i++){
			for(j = 0; j < iCols; j++){
				mtx_to_vec[i * iCols + j] = mtx[i][j];
			}
		}

	cudaMalloc((void**)&d_mtx_to_vec, iRows * iCols * sizeof(int));
//for(i=0; i < iRows*iCols; i++) {printf("%d\n",mtx_to_vec[i]);}
     cudaMalloc((void**)&d_mtx_to_vec, iRows * iCols * sizeof(int));
     unsigned long long int numThreads, steps, steps_remainder, Inner_num = (unsigned long long int) 1 << (iRows), copyNum;
	 
//	 printf("numThread: %d\n", numThread);

	 numThreads = (unsigned long long int) 10752; steps=Inner_num/numThreads; steps_remainder = Inner_num % numThreads;// ( (int) Inner_num/(numThreads));

	int *L2_vector, *d_L2_vector, L2_max = 0, *L2_strategy, *d_L2_strategy;

	copyNum = numThreads > Inner_num ? Inner_num : numThreads;
	printf("copyNum: %llu\n", copyNum);
	L2_vector = (int*) malloc(copyNum * sizeof(int));
	L2_strategy = (int*) malloc(copyNum * (iRows) * sizeof(int));

	cudaMalloc((void**)&d_L2_vector, copyNum * sizeof(int));
	cudaMalloc((void**)&d_L2_strategy, copyNum * (iRows) * sizeof(int));
	
	cudaMemcpy(d_mtx_to_vec, mtx_to_vec, iRows * iCols * sizeof(int), cudaMemcpyHostToDevice);
printf("NUM of Blocks: %llu\n", (unsigned long long int) ceil((int) numThreads/32));
	 func<<<(unsigned long long int) ceil((int) numThreads/32),32>>>(d_mtx_to_vec, numThreads, steps, steps_remainder, d_L2_vector, d_L2_strategy, iCols, iRows);
	cudaMemcpy(L2_vector, d_L2_vector, copyNum * sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(L2_strategy, d_L2_strategy, copyNum * (iRows) * sizeof(int), cudaMemcpyDeviceToHost);
	 
	for(i = 0; i < copyNum; i++){
		if(L2_max < L2_vector[i]) {L2_max = L2_vector[i]; iMax = i;}
	}
FILE *fp;
fp = fopen("strategy_L2.txt", "w");	
	printf("The strategy is:\n");
	for(i=0; i<(iRows); i++) {printf("%d\n", L2_strategy[iMax * (iRows) + i]); fprintf(fp, "%d\n", L2_strategy[iMax * (iRows) + i]);}
fclose(fp);

//     printf("szam: %d\n", h >> 1);
//     printf("szam: %d\n", h >> 1);
	printf("L2 is: %d\n",L2_max);

	free(L2_vector);
	free(L2_strategy);
	free(mtx_to_vec);

	cudaFree(d_L2_vector);
	cudaFree(d_L2_strategy);
	cudaFree(d_mtx_to_vec);

     return 0;     
}

/* Compile with: g++ -O2 main_det8_class30.c -o det8_class30
Authors: C. Beierle, G. Leander  -- Feb 2020

This program constructs and returns all 8-bit APN permutations F s.t. AF = FA with A corresponding to class 30, i.e.,
A = B = 
[1 0|0 0|0 0|0 0]
[0 1|0 0|0 0|0 0]
[---+---+---+---]
[0 0|0 1|0 0|0 0]
[0 0|1 0|0 0|0 0]
[---+---+---+---]
[0 0|0 0|0 1|0 0]
[0 0|0 0|1 0|0 0]
[---+---+---+---]
[0 0|0 0|0 0|0 1]
[0 0|0 0|0 0|1 0] 
up to affine equivalence. It may return multiple representatives in one affine-equivalence class. The program is optimized for that particular class.
*/

#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <stdio.h>
#include <time.h>

// dimension of the APN permutation to find
#define  N 8

time_t start;
time_t runtime;

// the matrix A
int A[(1<<N)] = {0, 1, 2, 3, 8, 9, 10, 11, 4, 5, 6, 7, 12, 13, 14, 15, 32, 33, 34, 35, 40, 41, 42, 43, 36, 37, 38, 39, 44, 45, 46, 47, 16, 17, 18, 19, 24, 25, 26, 27, 20, 21, 22, 23, 28, 29, 30, 31, 48, 49, 50, 51, 56, 57, 58, 59, 52, 53, 54, 55, 60, 61, 62, 63, 128, 129, 130, 131, 136, 137, 138, 139, 132, 133, 134, 135, 140, 141, 142, 143, 160, 161, 162, 163, 168, 169, 170, 171, 164, 165, 166, 167, 172, 173, 174, 175, 144, 145, 146, 147, 152, 153, 154, 155, 148, 149, 150, 151, 156, 157, 158, 159, 176, 177, 178, 179, 184, 185, 186, 187, 180, 181, 182, 183, 188, 189, 190, 191, 64, 65, 66, 67, 72, 73, 74, 75, 68, 69, 70, 71, 76, 77, 78, 79, 96, 97, 98, 99, 104, 105, 106, 107, 100, 101, 102, 103, 108, 109, 110, 111, 80, 81, 82, 83, 88, 89, 90, 91, 84, 85, 86, 87, 92, 93, 94, 95, 112, 113, 114, 115, 120, 121, 122, 123, 116, 117, 118, 119, 124, 125, 126, 127, 192, 193, 194, 195, 200, 201, 202, 203, 196, 197, 198, 199, 204, 205, 206, 207, 224, 225, 226, 227, 232, 233, 234, 235, 228, 229, 230, 231, 236, 237, 238, 239, 208, 209, 210, 211, 216, 217, 218, 219, 212, 213, 214, 215, 220, 221, 222, 223, 240, 241, 242, 243, 248, 249, 250, 251, 244, 245, 246, 247, 252, 253, 254, 255};

// Inverse of A
int A_inv[(1<<N)] = {0, 1, 2, 3, 8, 9, 10, 11, 4, 5, 6, 7, 12, 13, 14, 15, 32, 33, 34, 35, 40, 41, 42, 43, 36, 37, 38, 39, 44, 45, 46, 47, 16, 17, 18, 19, 24, 25, 26, 27, 20, 21, 22, 23, 28, 29, 30, 31, 48, 49, 50, 51, 56, 57, 58, 59, 52, 53, 54, 55, 60, 61, 62, 63, 128, 129, 130, 131, 136, 137, 138, 139, 132, 133, 134, 135, 140, 141, 142, 143, 160, 161, 162, 163, 168, 169, 170, 171, 164, 165, 166, 167, 172, 173, 174, 175, 144, 145, 146, 147, 152, 153, 154, 155, 148, 149, 150, 151, 156, 157, 158, 159, 176, 177, 178, 179, 184, 185, 186, 187, 180, 181, 182, 183, 188, 189, 190, 191, 64, 65, 66, 67, 72, 73, 74, 75, 68, 69, 70, 71, 76, 77, 78, 79, 96, 97, 98, 99, 104, 105, 106, 107, 100, 101, 102, 103, 108, 109, 110, 111, 80, 81, 82, 83, 88, 89, 90, 91, 84, 85, 86, 87, 92, 93, 94, 95, 112, 113, 114, 115, 120, 121, 122, 123, 116, 117, 118, 119, 124, 125, 126, 127, 192, 193, 194, 195, 200, 201, 202, 203, 196, 197, 198, 199, 204, 205, 206, 207, 224, 225, 226, 227, 232, 233, 234, 235, 228, 229, 230, 231, 236, 237, 238, 239, 208, 209, 210, 211, 216, 217, 218, 219, 212, 213, 214, 215, 220, 221, 222, 223, 240, 241, 242, 243, 248, 249, 250, 251, 244, 245, 246, 247, 252, 253, 254, 255};

int pos_fix[32] = {0,1,2,3,12,13,14,15,48,49,50,51,60,61,62,63,192,193,194,195,204,205,206,207,240,241,242,243,252,253,254,255}; // the 32 fixpoints of A

int pos[112] = {4,5,6,7,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,52,53,54,55,64,65,66,67,68,69,70,71,72,73,74,  75,  76,  77,  78,  79,  80,  81,  82,  83,  84,  85,  86,  87,  88,  89, 90, 91, 92, 93,  94,  95,  96,  97,  98,  99,  100,  101,  102,  103,  104,  105, 106,  107,  108,  109,  110,  111,  112,  113,  114,  115,  116,  117,  118,  119,  120,  121,  122,  123,  124,  125,  126,  127,  196,  197,  198,  199,  208, 209,  210,  211,  212,  213,  214,  215,  216,  217,  218,  219,  220,  221,  222,  223,  244,  245,  246,  247}; // one element of each cycle within the cycle structure of A (without the fixpoints). The next free position x will be selected from this list of positions

int repr[(1<<N)] = {0, 1, 2, 3, 4, 5, 6, 7, 4, 5, 6, 7, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 16, 17, 18, 19, 24, 25, 26, 27, 20, 21, 22, 23, 28, 29, 30, 31, 48, 49, 50, 51, 52, 53, 54, 55, 52, 53, 54, 55, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 64, 65, 66, 67, 72, 73, 74, 75, 68, 69, 70, 71, 76, 77, 78, 79, 96, 97, 98, 99, 104, 105, 106, 107, 100, 101, 102, 103, 108, 109, 110, 111, 80, 81, 82, 83, 88, 89, 90, 91, 84, 85, 86, 87, 92, 93, 94, 95, 112, 113, 114, 115, 120, 121, 122, 123, 116, 117, 118, 119, 124, 125, 126, 127, 192, 193, 194, 195, 196, 197, 198, 199, 196, 197, 198, 199, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 208, 209, 210, 211, 216, 217, 218, 219, 212, 213, 214, 215, 220, 221, 222, 223, 240, 241, 242, 243, 244, 245, 246, 247, 244, 245, 246, 247, 252, 253, 254, 255}; // the cycle representatives of each element from 0-255.

// numbers of even Hamming weight up to 255 (needed for more efficient APN check)
int evens[127] = {3, 5, 6, 9, 10, 12, 15, 17, 18, 20, 23, 24, 27, 29, 30, 33, 34, 36, 39, 40, 43, 45, 46, 48, 51, 53, 54, 57, 58, 60, 63, 65, 66, 68, 71, 72, 75, 77, 78, 80, 83, 85, 86, 89, 90, 92, 95, 96, 99, 101, 102, 105, 106, 108, 111, 113, 114, 116, 119, 120, 123, 125, 126, 129, 130, 132, 135, 136, 139, 141, 142, 144, 147, 149, 150, 153, 154, 156, 159, 160, 163, 165, 166, 169, 170, 172, 175, 177, 178, 180, 183, 184, 187, 189, 190, 192, 195, 197, 198, 201, 202, 204, 207, 209, 210, 212, 215, 216, 219, 221, 222, 225, 226, 228, 231, 232, 235, 237, 238, 240, 243, 245, 246, 249, 250, 252, 255};

// the known classes of 5-bit APN permutations by Leander and Brinkmann
int APN[5][32] = {
{0,1,2,4,3,6,8,16,5,10,15,27,19,29,31,20,7,18,25,21,12,14,24,28,26,11,23,13,30,9,17,22},
{0,1,2,4,3,8,13,16,5,11,21,31,23,15,19,30,6,28,29,9,24,27,14,18,10,17,12,26,7,25,20,22},
{0,1,2,4,3,8,13,16,5,17,28,27,30,14,24,10,6,19,11,20,31,29,12,21,18,26,15,25,7,22,23,9},
{0,1,2,4,3,8,16,28,5,10,25,17,18,23,31,29,6,20,13,24,19,11,9,22,27,7,14,21,26,12,30,15},
{0,1,2,4,3,8,16,28,5,10,26,18,17,20,31,29,6,21,24,12,22,15,25,7,14,19,13,23,9,30,27,11}
};

int solutions;
long long iterations;

int sbox[(1<<N)]; // it stores the (partial) S-box that we are going to construct (undefined values represented by -1)
int sbox_DDT[(1<<N)][(1<<N)]; // it stores the (partial) DDT of sbox

void printArray(int A[1<<N]) {
	for (int i=0;i<(1<<N);i++)
		printf("%2x ",i);
	printf("\n");
	for (int i=0;i<(1<<N);i++)
		if (A[i]!=-1) printf("%2x ",A[i]);
			else printf("-- ");
	printf("\n\n");
}


void printArray2(int A[1<<N]) {
	printf("{0x%02x",A[0]);
	for (int i=1;i<(1<<N);i++)
		printf(",0x%02x",A[i]);
	printf("};\n");
}

void fprintArray2(FILE *fp, int A[1<<N]) {
	fprintf(fp, "{0x%02x",A[0]);
	for (int i=1;i<(1<<N);i++)
		fprintf(fp, ",0x%02x",A[i]);
	fprintf(fp, "};\n");
}


// checks if A does not contain -1
int isComplete(int A[1<<N]) {
    int val = A[(1<<N)-1];
    int i;
    if (val==-1) return 0;
    A[(1<<N)-1] = -1;
    for (i=0; A[i] != -1; i++);
    A[(1<<N)-1] = val;
	return i==((1<<N)-1);

}

// multiplicative order of x for matrix
int order_mat(int A[1<<N], int x) {
	int y=A[x];
	int i=1;
	while(y!=x) {
		y=A[y];
		i++;
	}
	return i;
}

// return 1 if y is not yet in the image of sbox. returns 0 if it is already taken (i.e., check whether y is contained in sbox or not). It only checks whether the cycle representative is already contained. The function does not work for y=0, since it is not included in the pos array
int isNotTaken(int y) {
	for (int i=0; i<112; i++) {
		if (repr[sbox[pos[i]]] == repr[y])
			return 0;
	}
	return 1;
}

// continue building DDT for new set point sbox[c]. Returns 0 if the S-box cannot be APN anymore (because one entry is set to 4), returns 1 otherwise
int addDDTInformation(int c) {
    int a;
    // for all input differences
	for (int i=0;i<127;i++) {
		a = evens[i];	// we only need to check those input differences corresponding to vectors with even Hamming Weight
                // if sbox[c + a] is set, increase DDT entry by 2
                if (sbox[c^a]!= -1) {
                    sbox_DDT[a][sbox[c]^sbox[c^a]]+=2;
                    // if set above 2, return that it cannot be APN anymore
                    if (sbox_DDT[a][sbox[c]^sbox[c^a]]>2) return 0;
                }

	}
	return 1;
}

// reducing DDT entries for new removed point sbox[c]. It reduces all entries that were increased during the APN check of addDDTInformation(c).
void removeDDTInformation(int c) {
	int a;
	for (int i=0;i<127;i++) {
		a = evens[i];
                if (sbox[c^a]!= -1) {
                    sbox_DDT[a][sbox[c]^sbox[c^a]]-=2;
                    if (sbox_DDT[a][sbox[c]^sbox[c^a]]==2) break; // important because addDDTInformation(c) did not increment any values above that point.
                }

	}
}

int maxDepth=0;

// recursive construction of the sbox.
void nextValue_with_filter(int depth, char* filename) {
	int xS,yS,j;
	iterations++;

	if (depth>maxDepth) maxDepth=depth; // to show maxDepth
	if ((time(NULL)-start)>=5) { // display the status every 5 seconds
		start=time(NULL);
		printf(" depth:%d\n",depth);
		printArray(sbox);
		printf("solutions so far:%d maxDepth:%d \n",solutions,maxDepth);
		printf("iterations:%lli  \n",iterations);
		printf("running time: %li sec\n", time(NULL)-runtime);
		fflush(stdout);
	}

	//complete
	if (isComplete(sbox)) {
			solutions++;
			printf("found a new apn permutation: #%d\n",solutions);
			printArray2(sbox);
			fflush(stdout);
			FILE *fp = fopen(filename, "a");
			if (fp == NULL)
            {
                printf("Error opening file!\n");
		fflush(stdout);
            }
            else
            {
                fprintArray2(fp, sbox);
            }
            fclose(fp);
		return;
	}

	//not complete

	int x=pos[depth];	// get next free position in the look-up table
	for (int y=(1<<N)-1;y>=0;y--) {	
		if (isNotTaken(y) && (order_mat(A,y)!=1)) { // exclude the fixpoints of A
			xS = x;
			yS = y;
			for (int i=0; i<1; i++) {	// we know that the order of x is 2
				// we set sbox[xS] = yS. Add this to the partial DDT information and check whether it can still be APN
				sbox[xS]=yS;
				if (!addDDTInformation(xS)) goto UNDO_WITH_FILTER; // undo if it can't be APN anymore
				xS=A[xS];
                    		yS=A[yS];
			}
			sbox[xS]=yS;
			if (!addDDTInformation(xS)) goto UNDO_WITH_FILTER; // undo if it can't be APN anymore
			 	
                    	nextValue_with_filter(depth+1, filename);
                      
			//undo the changes (i.e., delete the set points of sbox from this step and undo the changes in the DDT)
			UNDO_WITH_FILTER:
			// we remove the points in reverse order to undo the DDT changes correctly
			do { 
				removeDDTInformation(xS);
				sbox[xS] = -1;
				xS = A_inv[xS];
			} while (xS != A_inv[x]);
		}
	}
}

void test_for_matrix(char* filename) {
	if (N!=8) {
		printf("invalid N\n");
		exit(-1);
	}
    int xS, yS;
    solutions=0;
    iterations=0;
    FILE *fp = fopen(filename, "w");
    if (fp == NULL)
    {
            printf("Error opening file!\n");
    }
    fclose(fp);
    printf("Search S-boxes invariant for matrix: \n");
    printArray2(A);
    printf("\n\n\n");
    fflush(stdout);

    // (1) start with the first APN permutation on 5 bit
    printf("Fix next APN on fixpoints \n"); 
    fflush(stdout);
    memset(sbox_DDT,0,(1<<N)*(1<<N)*sizeof(int));
    for (int x=0;x<(1<<N);x++) {
        sbox[x]=-1;
    }

    // w.l.o.g., we can set sbox on the fixed points to the APN permutation in 5 bit
    for (int x=0; x<32;x++) {
	sbox[pos_fix[x]] = pos_fix[APN[0][x]];
	if (!addDDTInformation(pos_fix[x])) exit(0);
    }

    // start the recursive search
    nextValue_with_filter(0, filename);


    // (2) start with the first APN permutation on 5 bit
    printf("Fix next APN on fixpoints \n"); 
    fflush(stdout);
    memset(sbox_DDT,0,(1<<N)*(1<<N)*sizeof(int));
    for (int x=0;x<(1<<N);x++) {
        sbox[x]=-1;
    }

    // w.l.o.g., we can set sbox on the fixed points to the APN permutation in 5 bit
    memset(sbox_DDT,0,(1<<N)*(1<<N)*sizeof(int));
    for (int x=0; x<32;x++) {
	sbox[pos_fix[x]] = pos_fix[APN[1][x]];
	if (!addDDTInformation(pos_fix[x])) exit(0);
    }

    // start the recursive search
    nextValue_with_filter(0, filename);


    // (3) start with the first APN permutation on 5 bit
    printf("Fix next APN on fixpoints \n"); 
    fflush(stdout);
    memset(sbox_DDT,0,(1<<N)*(1<<N)*sizeof(int));
    for (int x=0;x<(1<<N);x++) {
        sbox[x]=-1;
    }

    // w.l.o.g., we can set sbox on the fixed points to the APN permutation in 5 bit
    for (int x=0; x<32;x++) {
	sbox[pos_fix[x]] = pos_fix[APN[2][x]];
	if (!addDDTInformation(pos_fix[x])) exit(0);
    }

    // start the recursive search
    nextValue_with_filter(0, filename);


    // (4) start with the first APN permutation on 5 bit
    printf("Fix next APN on fixpoints \n"); 
    fflush(stdout);
    memset(sbox_DDT,0,(1<<N)*(1<<N)*sizeof(int));
    for (int x=0;x<(1<<N);x++) {
        sbox[x]=-1;
    }

    // w.l.o.g., we can set sbox on the fixed points to the APN permutation in 5 bit
    for (int x=0; x<32;x++) {
	sbox[pos_fix[x]] = pos_fix[APN[3][x]];
	if (!addDDTInformation(pos_fix[x])) exit(0);
    }

    // start the recursive search
    nextValue_with_filter(0, filename);


    // (5) start with the first APN permutation on 5 bit
    printf("Fix next APN on fixpoints \n"); 
    fflush(stdout);
    memset(sbox_DDT,0,(1<<N)*(1<<N)*sizeof(int));
    for (int x=0;x<(1<<N);x++) {
        sbox[x]=-1;
    }

    // w.l.o.g., we can set sbox on the fixed points to the APN permutation in 5 bit
    for (int x=0; x<32;x++) {
	sbox[pos_fix[x]] = pos_fix[APN[4][x]];
	if (!addDDTInformation(pos_fix[x])) exit(0);
    }

    // start the recursive search
    nextValue_with_filter(0, filename);
}


int main(int argc, char* argv[])
{
	if (argc != 2) {
        	printf("Usage: %s <outfile> ...\n", argv[0]);
        	return -1;
    	}
	srand (time(NULL));
	start = time(NULL);
	runtime = time(NULL);
	
	
	test_for_matrix(argv[1]);
	printf("\nFinished. Total running time: %li sec\n", time(NULL)-runtime);
	exit(0);
}

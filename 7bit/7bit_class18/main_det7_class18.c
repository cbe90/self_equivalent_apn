/* Compile with: g++ -O2 main_det7_class18.c -o det7_class18
Authors: C. Beierle, G. Leander  -- Feb 2020

This program constructs and returns all 7-bit APN permutations F s.t. AF = FB with A, B corresponding to class 18, i.e.,
A = 
0 0 1 0 0 0 0
1 0 1 0 0 0 0
0 1 0 0 0 0 0
0 0 0 0 0 0 1
0 0 0 1 0 0 0
0 0 0 0 1 0 1
0 0 0 0 0 1 1,
B =
0 0 1 0 0 0 0
1 0 0 0 0 0 0
0 1 1 0 0 0 0
0 0 0 0 0 0 1
0 0 0 1 0 0 1
0 0 0 0 1 0 1
0 0 0 0 0 1 0,
up to affine equivalence. It may return multiple representatives in one affine-equivalence class. 
*/

#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <stdio.h>
#include <time.h>
#include "commuting_matrices_A.h"
#include "commuting_matrices_B.h"

// dimension of the APN permutation to find
#define  N 7

time_t start;
time_t runtime;

// the matrix A
int A[(1<<N)] = {0, 2, 4, 6, 3, 1, 7, 5, 16, 18, 20, 22, 19, 17, 23, 21, 32, 34, 36, 38, 35, 33, 39, 37, 48, 50, 52, 54, 51, 49, 55, 53, 64, 66, 68, 70, 67, 65, 71, 69, 80, 82, 84, 86, 83, 81, 87, 85, 96, 98, 100, 102, 99, 97, 103, 101, 112, 114, 116, 118, 115, 113, 119, 117, 104, 106, 108, 110, 107, 105, 111, 109, 120, 122, 124, 126, 123, 121, 127, 125, 72, 74, 76, 78, 75, 73, 79, 77, 88, 90, 92, 94, 91, 89, 95, 93, 40, 42, 44, 46, 43, 41, 47, 45, 56, 58, 60, 62, 59, 57, 63, 61, 8, 10, 12, 14, 11, 9, 15, 13, 24, 26, 28, 30, 27, 25, 31, 29};

// the matrix B
int B[(1<<N)] = {0, 2, 4, 6, 5, 7, 1, 3, 16, 18, 20, 22, 21, 23, 17, 19, 32, 34, 36, 38, 37, 39, 33, 35, 48, 50, 52, 54, 53, 55, 49, 51, 64, 66, 68, 70, 69, 71, 65, 67, 80, 82, 84, 86, 85, 87, 81, 83, 96, 98, 100, 102, 101, 103, 97, 99, 112, 114, 116, 118, 117, 119, 113, 115, 56, 58, 60, 62, 61, 63, 57, 59, 40, 42, 44, 46, 45, 47, 41, 43, 24, 26, 28, 30, 29, 31, 25, 27, 8, 10, 12, 14, 13, 15, 9, 11, 120, 122, 124, 126, 125, 127, 121, 123, 104, 106, 108, 110, 109, 111, 105, 107, 88, 90, 92, 94, 93, 95, 89, 91, 72, 74, 76, 78, 77, 79, 73, 75};

// inverse of B
int B_inv[(1<<N)] = {0, 6, 1, 7, 2, 4, 3, 5, 88, 94, 89, 95, 90, 92, 91, 93, 8, 14, 9, 15, 10, 12, 11, 13, 80, 86, 81, 87, 82, 84, 83, 85, 16, 22, 17, 23, 18, 20, 19, 21, 72, 78, 73, 79, 74, 76, 75, 77, 24, 30, 25, 31, 26, 28, 27, 29, 64, 70, 65, 71, 66, 68, 67, 69, 32, 38, 33, 39, 34, 36, 35, 37, 120, 126, 121, 127, 122, 124, 123, 125, 40, 46, 41, 47, 42, 44, 43, 45, 112, 118, 113, 119, 114, 116, 115, 117, 48, 54, 49, 55, 50, 52, 51, 53, 104, 110, 105, 111, 106, 108, 107, 109, 56, 62, 57, 63, 58, 60, 59, 61, 96, 102, 97, 103, 98, 100, 99, 101};

// an ordering of the indices mod 2^n. It orders cycle-wise with respect to the cycles of B
int lex[(1<<N)] = {0,1, 2, 4, 5, 7, 3, 6, 8, 16, 32, 64, 56, 112, 88, 9, 18, 36, 69, 63, 115, 94, 10, 20, 37, 71, 59, 118, 89, 11, 22, 33, 66, 60, 117, 95, 12, 21, 39, 67, 62, 113, 90, 13, 23, 35, 70, 57, 114, 92, 14, 17, 34, 68, 61, 119, 91, 15, 19, 38, 65, 58, 116, 93, 24, 48, 96, 120, 72, 40, 80, 25, 50, 100, 125, 79, 43, 86, 26, 52, 101, 127, 75, 46, 81, 27, 54, 97, 122, 76, 45, 87, 28, 53, 103, 123, 78, 41, 82, 29, 55, 99, 126, 73, 42, 84, 30, 49, 98, 124, 77, 47, 83, 31, 51, 102, 121, 74, 44, 85, 105, 106, 108, 109, 111, 107, 110,104};

int pos[18] = {1,8,9,10,11,12,13,14,15,24,25,26,27,28,29,30,31,105}; // one element of each cycle within the cycle structure of B (without 0 and 104). The next free position x will be selected from this list of positions

int repr[(1<<N)] = {0, 1, 1, 1, 1, 1, 1, 1, 8, 9, 10, 11, 12, 13, 14, 15, 8, 13, 9, 12, 10, 15, 11, 14, 24, 25, 26, 27, 28, 29, 30, 31, 8, 15, 13, 10, 9, 14, 12, 11, 24, 30, 31, 25, 29, 27, 26, 28, 24, 29, 25, 28, 26, 31, 27, 30, 8, 12, 11, 15, 14, 10, 13, 9, 8, 14, 15, 9, 13, 11, 10, 12, 24, 28, 27, 31, 30, 26, 29, 25, 24, 27, 30, 29, 31, 28, 25, 26, 88, 89, 89, 89, 89, 89, 89, 89, 24, 31, 29, 26, 25, 30, 28, 27, 8, 11, 14, 13, 15, 12, 9, 10, 8, 10, 12, 14, 11, 9, 15, 13, 24, 26, 28, 30, 27, 25, 31, 29}; // the cycle representatives of cycles in A of each element from 0-127.

// numbers of even Hamming weight up to 127 (needed for more efficient APN check)
int evens[63] = {3, 5, 6, 9, 10, 12, 15, 17, 18, 20, 23, 24, 27, 29, 30, 33, 34, 36, 39, 40, 43, 45, 46, 48, 51, 53, 54, 57, 58, 60, 63, 65, 66, 68, 71, 72, 75, 77, 78, 80, 83, 85, 86, 89, 90, 92, 95, 96, 99, 101, 102, 105, 106, 108, 111, 113, 114, 116, 119, 120, 123, 125, 126};

int solutions;
long long iterations;
int n_skipped;
int skipped[128];

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

void printSkippedArray(int A[128]) {
	for (int i=0;i<32;i++)
		printf("%d ",A[i]);
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

void printDDT(int A[1<<N][1<<N]) {
    for (int i=0; i<(1<<N); i++) {
        for (int j=0; j<(1<<N); j++) {
            printf("%d ",A[i][j]);
        }
        printf("\n");
    }
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


// returns 1 if A is the smallest permutation up to conjugation with matrices in commuting_matrices. We use the lexicographic order of the lookup tables.
int is_smallest_in_class(int A[1<<N]){
    int E;
    for (int i=0;i<N_COMMUTING_MATRICES_B;i++) {
        for (int j=0; j<N_COMMUTING_MATRICES_A;j++) {
            for (int x=0;x<(1<<N);x++){
                if (A[lex[x]]==-1) goto NEXT_MAT; // We cannot say whether A is smaller than current matrix because of the undefined value
                if (A[commuting_matrices_B[i][lex[x]]] == -1) goto NEXT_MAT; // We cannot say whether A is smaller than current matrix because of the undefined value
                else {
                    E = commuting_matrices_A[j][A[commuting_matrices_B[i][lex[x]]]];
                    if (E < A[lex[x]]) return 0; // E is smaller, so A cannot be the smallest in class
                    if (E > A[lex[x]]) goto NEXT_MAT; // A is smaller, so consider next matrix
                } // consider next position
            }
            NEXT_MAT:;
        }
    }
    // A is the smallest of all, so return 1
    return 1;
}

// return 1 if y is not yet in the image of sbox. returns 0 if it is already taken (i.e., check whether y is contained in sbox or not). It only checks whether the cycle representative is already contained. The function does not work for y=0 (and y=88), since 0 and 127 are not included in the pos array
int isNotTaken(int y) {
	for (int i=0; i<18; i++) {
		if (repr[sbox[pos[i]]] == repr[y])
			return 0;
	}
	return 1;
}

// continue building DDT for new set point sbox[c]. Returns 0 if the S-box cannot be APN anymore (because one entry is set to 4), returns 1 otherwise
int addDDTInformation(int c) {
    int a;
    // for all input differences
	for (int i=0;i<63;i++) {
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
	for (int i=0;i<63;i++) {
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
	if ((time(NULL)-start)>=300) { // display the status every 300 seconds
		start=time(NULL);
		printf(" depth:%d\n",depth);
		printArray(sbox);
		printf("solutions so far:%d maxDepth:%d \n",solutions,maxDepth);
		printf("iterations:%lli  skipped: %d \nDepth skipped: \n",iterations,n_skipped);
		printSkippedArray(skipped);
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
	for (int y=(1<<N)-1;y>=1;y--) {		// do not include y=0 and y=88 since those values are already set in the beginning. The function isNotTaken(y) would not work for 0.
	if (y != 88) {
		if (isNotTaken(y)) { // if the orders are different, it cannot be of the desired form (we know the order of each element is 7)
			xS = x;
			yS = y;
			for (int i=0; i<6; i++) {	// we know the order of each element is 7
				// we set sbox[xS] = yS. Add this to the partial DDT information and check whether it can still be APN
				sbox[xS]=yS;
				if (!addDDTInformation(xS)) goto UNDO_WITH_FILTER; // undo if it can't be APN anymore
				xS=B[xS];
                    		yS=A[yS];
			}
			sbox[xS]=yS;
			if (!addDDTInformation(xS)) goto UNDO_WITH_FILTER; // undo if it can't be APN anymore
			
                	if (depth <= 4) { // if it is too deep, it is faster to omit the check for smallest representative.
                        	if (is_smallest_in_class(sbox)) nextValue_with_filter(depth+1, filename);
                        	else {n_skipped++; skipped[depth]++;}
                	}
                	else { // omit check for smallest representative
                    		nextValue_with_filter(depth+1, filename);
                	}
            
			//undo the changes (i.e., delete the set points of sbox from this step and undo the changes in the DDT)
			UNDO_WITH_FILTER:
			// we remove the points in reverse order to undo the DDT changes correctly
			do { 
				removeDDTInformation(xS);
				sbox[xS] = -1;
				xS = B_inv[xS];
			} while (xS != B_inv[x]);
		}
	}
	}
}

void test_for_matrix(char* filename, int first) {
	if (N!=7) {
		printf("invalid N\n");
		exit(-1);
	}
    int xS, yS;
    solutions=0;
    iterations=0;
    memset(sbox_DDT,0,(1<<N)*(1<<N)*sizeof(int));
    for (int i=0; i<128; i++) skipped[i]=0;
    FILE *fp = fopen(filename, "w");
    if (fp == NULL)
    {
            printf("Error opening file!\n");
    }
    fclose(fp);
    printf("Search S-boxes invariant for matrices: \n");
    printArray2(B);
    printArray2(A);
    printf("\n\n\n");
    fflush(stdout);
    for (int x=0;x<(1<<N);x++) {
        sbox[x]=-1;
    }

    // w.l.o.g., we can set sbox[0]=0
    sbox[0]=0;
    if (!addDDTInformation(0)) exit(0);
    // w.l.o.g., we can set sbox[104]=88
    sbox[104]=88;
    if (!addDDTInformation(104)) exit(0);

    sbox[1]=first;
    if (!addDDTInformation(1)) exit(0);
    xS=B[1];
    yS=A[first];
    // set the orbits accordingly
    for (int i=0;i<6;i++) {
                    sbox[xS]=yS;
                    if (!addDDTInformation(xS)) exit(0);
                    xS=B[xS];
                    yS=A[yS];
    } 
    if (!is_smallest_in_class(sbox)) {
	exit(0);
    } 

    // start the recursive search
    nextValue_with_filter(1, filename);
}

int main(int argc, char* argv[])
{
	if (argc != 3) {
        	printf("Usage: %s <outfile> <SBox[1]> ...\n", argv[0]);
        	return -1;
    	}
	srand (time(NULL));
	start = time(NULL);
	runtime = time(NULL);
	
	
	if ((atoi(argv[2]) == 0) or (atoi(argv[2]) == 88)) {
		printf("No permutation ..\n");
		return -1;
	}
	test_for_matrix(argv[1],atoi(argv[2]));
	printf("\nFinished. Total running time: %li sec\n", time(NULL)-runtime);
	exit(0);
}

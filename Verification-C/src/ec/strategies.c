#include <stdlib.h>//malloc and free
#include <strategies.h>
//#include <stdio.h>

void optimised_strategies(unsigned int **strategies, float *cost, unsigned int n, const float mul_c, const float eval_c)
{
	/* Simple optimal strategies for dimension 1 isogenies. 
	Adapted to higher dimensions if we omit gluing and splitting, as shown in https://ia.cr/2023/508. 

	Returns all the strategies for isgogeny chains of length <=n.
	strategy[i] contains a strategy for a chain of length i+1.*/
	unsigned int b;
	float min_cost, cur_cost;

	strategies[0]=(unsigned int *)malloc(0);
	cost[0]=0;
	for(int i=1;i<n;i++){
		b=0;
		min_cost=cost[i-1]+cost[0]+(b+1)*mul_c+(i-b)*eval_c;
		for(int j=1;j<i;j++){
			cur_cost=cost[i-1-j]+cost[j]+(j+1)*mul_c+(i-j)*eval_c;
			if(cur_cost<min_cost){
				b=j;
				min_cost=cur_cost;
			}
		}
		cost[i]=min_cost;

		strategies[i]=(unsigned int *)malloc(i*sizeof(unsigned int));
		strategies[i][0]=b+1;
		for(int k=0;k<i-1-b;k++){
			strategies[i][k+1]=strategies[i-1-b][k];
		}
		for(int k=0;k<b;k++){
			strategies[i][k+i-b]=strategies[b][k];
		}
	}
}

void optimised_strategy(unsigned int *S,unsigned int n, const float mul_c, const float eval_c)
{
	// Returns one strategy for an isogeny chain of length n.
	const unsigned int cn=n;
	unsigned int *strategies[cn];
	float cost[cn];
	optimised_strategies(strategies,cost,n,mul_c,eval_c);

	for(int i=0;i<n-1;i++){
		S[i]=strategies[n-1][i];
	}

	for(int i=0;i<n;i++){
		free(strategies[i]);
	}
}

void optimised_strategies_with_first_eval(unsigned int **strategies_left, float *cost_left,unsigned int **strategies_right, 
	float *cost_right, unsigned int n, const float mul_c, const float eval_c, const float first_eval_c)
{
	/* Returns all strategies for isogeny chains of length <=n with constraint at the beginning 
	(no scalar multiplication on the first codomain). Suitable for higher dimensional isogenies
	as explained in https://arxiv.org/pdf/2407.15492, Appendix E.2.1. */

	unsigned int b;
	float min_cost, cur_cost;

	optimised_strategies(strategies_right,cost_right,n-2,mul_c,eval_c);

	strategies_left[0]=(unsigned int *)malloc(0);
	strategies_left[1]=(unsigned int *)malloc(sizeof(unsigned int));
	strategies_left[1][0]=1;
	cost_left[0]=0;
	cost_left[1]=mul_c+first_eval_c;
	for(int i=2;i<n;i++){
		b=0;
		min_cost=cost_left[i-1]+cost_right[0]+(b+1)*mul_c+(i-1-b)*eval_c+first_eval_c;
		for(int j=1;j<i-1;j++){
			cur_cost=cost_left[i-1-j]+cost_right[j]+(j+1)*mul_c+(i-1-j)*eval_c+first_eval_c;
			if(cur_cost<min_cost){
				b=j;
				min_cost=cur_cost;
			}
		}
		cost_left[i]=min_cost;

		strategies_left[i]=(unsigned int *)malloc(i*sizeof(unsigned int));
		strategies_left[i][0]=b+1;
		for(int k=0;k<i-1-b;k++){
			strategies_left[i][k+1]=strategies_left[i-1-b][k];
		}
		for(int k=0;k<b;k++){
			strategies_left[i][k+i-b]=strategies_right[b][k];
		}
	}
}

void optimised_strategy_with_first_eval(unsigned int *S, unsigned int n, const float mul_c, 
	const float eval_c, const float first_eval_c)
{
	/* Returns one strategy for isogeny chains of length n with constraint at the beginning 
	(no scalar multiplication on the first codomain). */

	const unsigned int cn=n;
	unsigned int *strategies_left[cn], *strategies_right[cn-2];
	float cost_left[cn], cost_right[cn-2];
	optimised_strategies_with_first_eval(strategies_left,cost_left,strategies_right,cost_right,n,mul_c,eval_c,first_eval_c);

	for(int i=0;i<n-1;i++){
		S[i]=strategies_left[n-1][i];
	}

	for(int i=0;i<n;i++){
		free(strategies_left[i]);
	}
	for(int i=0;i<n-2;i++){
		free(strategies_right[i]);
	}
}


void optimised_strategy_with_first_eval_and_splitting(unsigned int *S, unsigned int n, unsigned int m, const float mul_c, 
	const float eval_c, const float first_eval_c)
{
	const unsigned int cn=n, cm=m;
	unsigned int *strategies_left[cn];
	unsigned int *strategies_middle[cn-cm-1];
	float cost_left[cn], cost_middle[cn-cm-1];
	// strat_R_split[l] contains all strategies of length <=n-2 with constraint l+1 steps before the end with l\in{0,...,m-1}.
	unsigned int **strat_R_split[cm];
	float *cost_R_split[cm];
	unsigned int b, l;
	float min_cost, cur_cost;

	optimised_strategies_with_first_eval(strategies_left, cost_left, strategies_middle, cost_middle, n-m+1, mul_c, 
		eval_c, first_eval_c);

	// Filling in strat_R_split
	for(l=0;l<m;l++){
		strat_R_split[l]=(unsigned int **)malloc((n-m+l)*sizeof(void *));
		cost_R_split[l]=(float *)malloc((n-m+l)*sizeof(float));

		strat_R_split[l][0]=(unsigned int *)malloc(0);
		cost_R_split[l][0]=0;
		if(l==0){
			for(int i=1;i<n-m+l;i++){
				strat_R_split[l][i]=(unsigned int *)malloc(i*sizeof(unsigned int));
				b=0;
				min_cost=cost_middle[i-1-b]+cost_R_split[0][b]+(b+1)*mul_c+(i-b)*eval_c;
				for(int j=2;j<i;j++){
					cur_cost=cost_middle[i-1-j]+cost_R_split[0][j]+(j+1)*mul_c+(i-j)*eval_c;
					if(cur_cost<min_cost){
						min_cost=cur_cost;
						b=j;
					}
				}

				cost_R_split[l][i]=min_cost;

				strat_R_split[l][i][0]=b+1;
				for(int k=0; k<i-1-b;k++){
					strat_R_split[l][i][k+1]=strategies_middle[i-1-b][k];
				}
				for(int k=0; k<b;k++){
					strat_R_split[l][i][k+i-b]=strat_R_split[l][b][k];
				}
			}
		}
		else{
			for(int i=1;i<n-m+l;i++){
				strat_R_split[l][i]=(unsigned int *)malloc(i*sizeof(unsigned int));
				b=0;
				min_cost=cost_R_split[l-b-1][i-1-b]+cost_middle[b]+(b+1)*mul_c+(i-b)*eval_c;
				for(int j=1;j<i;j++){
					if(j<l){
						cur_cost=cost_R_split[l-j-1][i-1-j]+cost_middle[j]+(j+1)*mul_c+(i-j)*eval_c;
					}
					if(j>l){
						cur_cost=cost_middle[i-1-j]+cost_R_split[l][j]+(j+1)*mul_c+(i-j)*eval_c;
					}
					if(cur_cost<min_cost){
						min_cost=cur_cost;
						b=j;
					}
				}

				cost_R_split[l][i]=min_cost;

				strat_R_split[l][i][0]=b+1;
				if(b<l){
					for(int k=0; k<i-1-b;k++){
						strat_R_split[l][i][k+1]=strat_R_split[l-b-1][i-1-b][k];
					}
					for(int k=0; k<b;k++){
						strat_R_split[l][i][k+i-b]=strategies_middle[b][k];
					}
				}
				else{
					for(int k=0; k<i-1-b;k++){
						strat_R_split[l][i][k+1]=strategies_middle[i-1-b][k];
					}
					for(int k=0; k<b;k++){
						strat_R_split[l][i][k+i-b]=strat_R_split[l][b][k];
					}
				}
			}
		}
	}



	// Filling in the last element of strategies_left
	for(int i=n-m;i<n;i++){
		strategies_left[i]=(unsigned int *)malloc(i*sizeof(unsigned int));
		b=0;
		l=i-(n-m);
		min_cost=cost_left[i-1-b]+cost_R_split[l][b]+(b+1)*mul_c+(i-1-b)*eval_c+first_eval_c;
		for(int j=1;j<i;j++){
			if((j!=l)&&(j!=n-1)){
				cur_cost=cost_left[i-1-j]+cost_R_split[l][j]+(j+1)*mul_c+(i-1-j)*eval_c+first_eval_c;
			}
			if(cur_cost<min_cost){
				min_cost=cur_cost;
				b=j;
			}
		}
		cost_left[i]=min_cost;

		strategies_left[i][0]=b+1;
		for(int k=0;k<i-1-b;k++){
			strategies_left[i][k+1]=strategies_left[i-1-b][k];
		}
		for(int k=0;k<b;k++){
			strategies_left[i][k+i-b]=strat_R_split[l][b][k];
		}
	}

	//Filling in S
	for(int i=0;i<n;i++){
		S[i]=strategies_left[n-1][i];
	}


	// Cleaning-up memory
	for(int l=0;l<m;l++){
		free(cost_R_split[l]);
		for(int i=1;i<n-m+l;i++){
			free(strat_R_split[l][i]);
		}
		free(strat_R_split[l]);
	}
	for(int i=0;i<n;i++){
		free(strategies_left[i]);
	}
}

//int main(){
	//int n=25;
	//int tab[n];
	//unsigned int S0[127], S1[127], S2[127];
	//optimised_strategy(S0,128,1.0,1.0);
	//optimised_strategy_with_first_eval(S1, 128, 1.0, 1.0, 10.0);
	//optimised_strategy_with_first_eval_and_splitting(S2,128,4,1.0,1.0,10.0);
	//printf("Simple strategy:\n");
	//for(int i=0;i<127;i++){
		//printf("%i\n",S0[i]);
	//}
	//printf("Strategy with constraint at the beginning:\n");
	//for(int i=0;i<127;i++){
		//printf("%i\n",S1[i]);
	//}
	//printf("Strategy with constraint at the beginning and m steps before the end:\n");
	//for(int i=0;i<127;i++){
		//printf("%i\n",S2[i]);
	//}
//}
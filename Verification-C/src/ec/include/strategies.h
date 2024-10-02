#ifndef STRATEGIES_H
#define STRATEGIES_H

void optimised_strategies(unsigned int **strategies, float *cost, unsigned int n, const float mul_c, const float eval_c);
void optimised_strategy(unsigned int *S,unsigned int n, const float mul_c, const float eval_c);
void optimised_strategies_with_first_eval(unsigned int **strategies_left, float *cost_left,unsigned int **strategies_right, 
	float *cost_right, unsigned int n, const float mul_c, const float eval_c, const float first_eval_c);
void optimised_strategy_with_first_eval(unsigned int *S, unsigned int n, const float mul_c, 
	const float eval_c, const float first_eval_c);
void optimised_strategy_with_first_eval_and_splitting(unsigned int *S, unsigned int n, unsigned int m, const float mul_c, 
	const float eval_c, const float first_eval_c);

#endif
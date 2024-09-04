#include <stdio.h>
#include <string.h>
#include <mp.h>

bool test_div(){

	{// Prime pHD256
	digit_t pp1[4] = {0x0000000000000000, 0x4000000000000000, 0xa4382e87ff9dc589, 0x2827baebd5c8e56e};
	digit_t cofactor_3[4]={0x0000000000000000, 0x4000000000000000, 0x3, 0x0};
	digit_t power_3[4]={0x94fd9829d87f5079, 0xc5afe6ff302bcbf, 0x0, 0x0};

	digit_t pp1_div_power_3[4];

	mp_div(pp1_div_power_3,pp1,power_3,4);

	if(mp_compare(pp1_div_power_3,cofactor_3,4)!=0){
		return 0;
	}}

	{// Prime pHD384
	digit_t pp1[6] = {0x0, 0x0, 0x8000000000000000, 0x21d47fd4347f03cc, 0xf96d653c9a6b76df, 0x412508a24fc64c3};
	digit_t cofactor_3[6] = {0x0, 0x0, 0x8000000000000000, 0x0, 0x0, 0x0};
	digit_t power_3[6] = {0x43a8ffa868fe0799, 0xf2daca7934d6edbe, 0x824a11449f8c987, 0x0, 0x0, 0x0};

	digit_t pp1_div_power_3[6];

	mp_div(pp1_div_power_3,pp1,power_3,6);

	if(mp_compare(pp1_div_power_3,cofactor_3,6)!=0){
		return 0;
	}}

	{// Prime pHD512
	digit_t pp1[8] = {0x0, 0x0, 0x0, 0x0, 0x80429e2d58ebb467, 0x52e830945c1a0446, 0x5eec151f0a69c447, 0xa65f4ee938387923};
	digit_t cofactor_3[8] = {0x0, 0x0, 0x0, 0x0, 0x1f, 0x0, 0x0, 0x0};
	digit_t power_3[8] = {0xeb5cfcd82c28a2b9, 0x4cff3b5f9fdfce96, 0xb07b3a7cdf4dbc02, 0x55de9c5756d2d32, 0x0, 0x0, 0x0, 0x0};

	digit_t pp1_div_power_3[8];

	mp_div(pp1_div_power_3,pp1,power_3,8);

	if(mp_compare(pp1_div_power_3,cofactor_3,8)!=0){
		return 0;
	}}

	return 1;
}

int main(){
	bool b;
	b=test_div();

	if(b){
		printf("Test succeeded.\n");
	}
	else{
		printf("Test failed.\n");
	}

	return 0;
}
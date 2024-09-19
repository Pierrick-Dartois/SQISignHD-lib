#include <fp.h>
#include <fp2.h>
#include <test_utils.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void nqr_val(fp2_t *x){
	fp2_random_test(x);
	while(fp2_is_square(x)){
		fp2_random_test(x);
	}
}

void z_nqr_val(fp2_t *x){
	fp2_t xm1, one;
	fp2_set_one(&one);

	fp2_random_test(x);
	fp2_sub(&xm1,x,&one);

	while((!fp2_is_square(x))|(fp2_is_square(&xm1))){
		fp2_random_test(x);
		fp2_sub(&xm1,x,&one);
	}
}

uint16_t compute_cofactor(uint64_t * cofactor,const uint64_t * p){
	uint64_t one[NWORDS_FIELD];
	one[0]=1;
	for(int i=1;i<NWORDS_FIELD;i++){
		one[i]=0;
	}
	mp_add(cofactor,p,one,NWORDS_FIELD);
	uint16_t power_of_2=0;
	while(!(cofactor[0]&1)){
		mp_shiftr(cofactor,1,NWORDS_FIELD);
		power_of_2++;
	}
	return power_of_2;
}

uint16_t compute_cofactor_3g(uint64_t * cofactor,const uint64_t * p){
	uint64_t q[NWORDS_FIELD], r[NWORDS_FIELD], one[NWORDS_FIELD], three[NWORDS_FIELD];

	mp_set_small(one,1,NWORDS_FIELD);
	mp_add(cofactor,p,one,NWORDS_FIELD);

	mp_set_small(three,3,NWORDS_FIELD);
	uint16_t power_of_3=0;

	mp_div_with_remainder(q,r,cofactor,three,NWORDS_FIELD);
	while(mp_is_zero(r,NWORDS_FIELD)==1){
		mp_copy(cofactor,q,NWORDS_FIELD);
		mp_div_with_remainder(q,r,cofactor,three,NWORDS_FIELD);
		power_of_3++;
	}
	return power_of_3;
}

uint16_t compute_nbits(uint64_t * cofactor){
	uint16_t nbits=NWORDS_FIELD*RADIX;
	uint8_t is_one=0;
	uint64_t mask;
	for(int i=0;i<NWORDS_FIELD;i++){
		for(int j=0;j<RADIX;j++){
			// 1ULL otherwise, he cannot shift more than 32 bits... C is so wierd...
			mask=1ULL<<(RADIX-1-j);
			is_one=(mask&cofactor[NWORDS_FIELD-1-i])>>(RADIX-1-j);
			if(is_one){
				return nbits;
			}
			else{
				nbits--;
			}
		}
	}
	return nbits;
}

void write_cfile(const char filename[],const int num){
	FILE *fptr;
	fptr=fopen(filename,"w");

	fprintf(fptr,"%s\n","#include <stddef.h>");
	fprintf(fptr,"%s\n","#include <stdint.h>");
	fprintf(fptr,"%s\n","#include <tutil.h>");
	fprintf(fptr,"%s\n","#include <fp2.h>");
	fprintf(fptr,"%s\n\n","#include \"gf_constants.h\"");

	fprintf(fptr,"%s","const fp2_t NQR_TABLE[");
	fprintf(fptr,"%i",num);
	fprintf(fptr,"%s","] = {");

	fp2_t x;

	for(int i=0;i<num;i++){
		nqr_val(&x);
		//printf("%llu\n",x.re[0]);
		//printf("%llu\n",x.re[0]&1);
		fprintf(fptr,"%s","{{");
		for(int j=0;j<NWORDS_FIELD;j++){
			fprintf(fptr,"0x%01llx",x.re[j]);
			if(j<NWORDS_FIELD-1){
				fprintf(fptr,"%s",",");
			}
			else{
				fprintf(fptr,"%s","}, {");
			}
		}

		for(int j=0;j<NWORDS_FIELD;j++){
			fprintf(fptr,"0x%01llx",x.im[j]);
			if(j<NWORDS_FIELD-1){
				fprintf(fptr,"%s",",");
			}
			else if(i<num-1){
				fprintf(fptr,"%s","}}, ");
			}
			else{
				fprintf(fptr,"%s\n","}}};");
			}
		}
	}

	fprintf(fptr,"%s","const fp2_t Z_NQR_TABLE[");
	fprintf(fptr,"%i",num);
	fprintf(fptr,"%s","] = {");

	for(int i=0;i<num;i++){
		z_nqr_val(&x);
		fprintf(fptr,"%s","{{");
		for(int j=0;j<NWORDS_FIELD;j++){
			fprintf(fptr,"0x%01llx",x.re[j]);
			if(j<NWORDS_FIELD-1){
				fprintf(fptr,"%s",",");
			}
			else{
				fprintf(fptr,"%s","}, {");
			}
		}

		for(int j=0;j<NWORDS_FIELD;j++){
			fprintf(fptr,"0x%01llx",x.im[j]);
			if(j<NWORDS_FIELD-1){
				fprintf(fptr,"%s",",");
			}
			else if(i<num-1){
				fprintf(fptr,"%s","}}, ");
			}
			else{
				fprintf(fptr,"%s\n","}}};");
			}
		}
	}

	uint64_t cofactor[NWORDS_FIELD];
	uint16_t power_of_2, nbits;
	power_of_2=compute_cofactor(cofactor,p);
	nbits=compute_nbits(cofactor);
	
	uint8_t NWORDS_P_COFACTOR_FOR_2F;
	//Ceiling formula
	NWORDS_P_COFACTOR_FOR_2F=(nbits/RADIX)+!!(nbits%RADIX);

	fprintf(fptr,"\n%s","const digit_t p_cofactor_for_2f[");
	fprintf(fptr,"%u",NWORDS_P_COFACTOR_FOR_2F);
	fprintf(fptr,"%s","] = {");
	for(int i=0;i<NWORDS_P_COFACTOR_FOR_2F;i++){
		fprintf(fptr,"0x%01llx",cofactor[i]);
		if(i<NWORDS_P_COFACTOR_FOR_2F-1){
			fprintf(fptr,"%s",", ");
		}
		else{
			fprintf(fptr,"%s\n","};");
		}
	}

	uint64_t cofactor_3g[NWORDS_FIELD];
	uint16_t power_of_3, nbits_3;
	power_of_3=compute_cofactor_3g(cofactor_3g,p);
	nbits_3=compute_nbits(cofactor_3g);
	
	uint8_t NWORDS_P_COFACTOR_FOR_3G;
	//Ceiling formula
	NWORDS_P_COFACTOR_FOR_3G=(nbits_3/RADIX)+!!(nbits_3%RADIX);

	fprintf(fptr,"\n%s","const digit_t p_cofactor_for_3g[");
	fprintf(fptr,"%u",NWORDS_P_COFACTOR_FOR_3G);
	fprintf(fptr,"%s","] = {");
	for(int i=0;i<NWORDS_P_COFACTOR_FOR_3G;i++){
		fprintf(fptr,"0x%01llx",cofactor_3g[i]);
		if(i<NWORDS_P_COFACTOR_FOR_3G-1){
			fprintf(fptr,"%s",", ");
		}
		else{
			fprintf(fptr,"%s\n","};");
		}
	}

	fprintf(fptr,"%s","const uint16_t P_COFACTOR_FOR_2F_BITLENGTH = ");
	fprintf(fptr,"%hu",nbits);
	fprintf(fptr,"%s\n",";");

	fprintf(fptr,"%s","const uint16_t POWER_OF_2 = ");
	fprintf(fptr,"%hu",power_of_2);
	fprintf(fptr,"%s\n",";");

	fprintf(fptr,"%s","const uint16_t P_COFACTOR_FOR_3G_BITLENGTH = ");
	fprintf(fptr,"%hu",nbits_3);
	fprintf(fptr,"%s\n",";");

	fprintf(fptr,"%s","const uint16_t POWER_OF_3 = ");
	fprintf(fptr,"%hu",power_of_3);
	fprintf(fptr,"%s\n",";");
}

int main(){
	write_cfile(FILENAME,20);
	return 0;
}


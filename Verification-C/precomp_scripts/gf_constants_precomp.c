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

	while((~fp2_is_square(x))|(fp2_is_square(&xm1))){
		fp2_random_test(x);
		fp2_sub(&xm1,x,&one);
	}
}

void write_cfile(const char filename[],const int num){
	FILE *fptr;
	fptr=fopen(filename,"w");

	fprintf(fptr,"%s\n","#include <stddef.h>");
	fprintf(fptr,"%s\n","#include <stdint.h>");
	fprintf(fptr,"%s\n","#include <tutil.h>");
	fprintf(fptr,"%s\n","#include <fp2.h>");
	fprintf(fptr,"%s\n\n","#include \"gf_constants.h\"");

	fprintf(fptr,"%s","const const fp2_t NQR_TABLE[");
	fprintf(fptr,"%i",num);
	fprintf(fptr,"%s","] = {");

	fp2_t x;

	for(int i=0;i<num;i++){
		nqr_val(&x);
		fprintf(fptr,"%s","{{");
		for(int j=0;j<NWORDS_FIELD;j++){
			fprintf(fptr,"%llu",x.re[j]);
			if(j<NWORDS_FIELD-1){
				fprintf(fptr,"%s",",");
			}
			else{
				fprintf(fptr,"%s","}, {");
			}
		}

		for(int j=0;j<NWORDS_FIELD;j++){
			fprintf(fptr,"%llu",x.im[j]);
			if(j<NWORDS_FIELD-1){
				fprintf(fptr,"%s",",");
			}
			else{
				fprintf(fptr,"%s","}}, ");
			}
		}
	}
	fprintf(fptr,"%s\n","};");

	fprintf(fptr,"%s","const const fp2_t Z_NQR_TABLE[");
	fprintf(fptr,"%i",num);
	fprintf(fptr,"%s","] = {");

	for(int i=0;i<num;i++){
		z_nqr_val(&x);
		fprintf(fptr,"%s","{{");
		for(int j=0;j<NWORDS_FIELD;j++){
			fprintf(fptr,"%llu",x.re[j]);
			if(j<NWORDS_FIELD-1){
				fprintf(fptr,"%s",",");
			}
			else{
				fprintf(fptr,"%s","}, {");
			}
		}

		for(int j=0;j<NWORDS_FIELD;j++){
			fprintf(fptr,"%llu",x.im[j]);
			if(j<NWORDS_FIELD-1){
				fprintf(fptr,"%s",",");
			}
			else{
				fprintf(fptr,"%s","}}, ");
			}
		}
	}
	fprintf(fptr,"%s","};");

}

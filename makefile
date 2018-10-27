RCF : source/reformat.c source/contact_frequency.c makefile
	cc source/reformat.c -o source/reformat -lm
	cc source/contact_frequency.c -o source/contact_frequency -lm

clean :
	rm -f source/reformat source/contact_frequency

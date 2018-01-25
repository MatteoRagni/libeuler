newton:
	gcc -I. -g libnewton.c test/newton_test.c -llapacke  -llapack -lblas -lm -o newton_test

eulere:
	gcc -I. -g libnewton.c libeuler.c test/eulere_test.c -llapacke  -llapack -lblas -lm -o eulere_test

euleri:
	gcc -I. -g libnewton.c libeuler.c test/euleri_test.c -llapacke  -llapack -lblas -lm -o euleri_test

debug:
	gdb --tui ./test
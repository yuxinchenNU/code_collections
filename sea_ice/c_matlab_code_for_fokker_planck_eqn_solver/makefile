all: FPweno

runFPweno: FPweno
	FPweno FP.in

FPweno: FPweno.o
	gcc -o FPweno FPweno.o -framework Accelerate -lm

FPweno.o: FPweno.c
	gcc -c FPweno.c

clean:
	rm *.o FPweno


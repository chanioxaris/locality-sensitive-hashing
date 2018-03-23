OBJS = structs.o main.o functions.o metric_functions.o output_functions.o hashtable.o preprocessing.o search.o quicksort.o
EXEC = lsh
GCCC = gcc -c
GCCO = gcc -o
LIBS = -lm

$(EXEC): $(OBJS)
	$(GCCO) $(EXEC) $(OBJS) $(LIBS)

structs.o: structs.c
	$(GCCC) structs.c

main.o: main.c
	$(GCCC) main.c

metric_functions.o: metric_functions.c
	$(GCCC) metric_functions.c

output_functions.o: output_functions.c
	$(GCCC) output_functions.c

hashtable.o: hashtable.c
	$(GCCC) hashtable.c

functions.o: functions.c
	$(GCCC) functions.c

preprocessing.o: preprocessing.c
	$(GCCC) preprocessing.c

quicksort.o: quicksort.c
	$(GCCC) quicksort.c

search.o: search.c
	$(GCCC) search.c	

clean:
	rm -rf $(OBJS) $(EXEC)
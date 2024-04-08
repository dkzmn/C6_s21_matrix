FLAGS=-Wall -Wextra -Werror -std=c11# -fsanitize=address -pedantic
LIBS=-lcheck
OS:=$(shell uname -s)
ifeq ($(OS), Darwin)
	FLAGS += -D MACOS
else
	LIBS += -lsubunit -lrt -lm -lpthread
endif

rebuild: clean all

all: test gcov_report

s21_matrix.a: s21_matrix.c
	gcc $(FLAGS) -c s21_matrix.c -o s21_matrix.o
	ar rc s21_matrix.a s21_matrix.o
	ranlib s21_matrix.a

test: test_s21_matrix.c s21_matrix.a
	gcc $(FLAGS) -c test_s21_matrix.c -o test_s21_matrix.o
	gcc $(FLAGS) test_s21_matrix.o s21_matrix.a $(LIBS) -o test
	./test

gcov_report:
	gcc $(FLAGS) -c s21_matrix.c --coverage -o s21_matrix.o
	gcc $(FLAGS) -c test_s21_matrix.c -o test_s21_matrix.o
	gcc $(FLAGS) test_s21_matrix.o s21_matrix.o --coverage $(LIBS) -o test
	./test
	lcov -t "s21_matrix" -o s21_matrix.info -c -d .
	genhtml -o report s21_matrix.info
	rm -rf *.gcda *.gcno *.info

clean:
	rm -rf *.o *.a *.css *.html *.gcda *.gcno *.info ./report test

clf:
	clang-format -n --style=Google *.[ch]

cpk:
	cppcheck *.[ch]

style: clf cpk
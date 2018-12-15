FLAGS = -pedantic -Wall -Wextra -std=c++14 -g

all:launchsim
clean:rm -rf *.o

# Dependencies
main.o : rocket.hpp

%.o : %.cpp
	g++ $(FLAGS) -c $< -o $@

launchsim: main.o
	g++ $(FLAGS) $^ -o $@
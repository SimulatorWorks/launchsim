FLAGS = -pedantic -Wall -Wextra -std=c++14 -g

all:launchsim
clean:rm -rf *.o

# Dependencies
rocket.o : rocket.hpp
main.o : rocket.hpp

%.o : %.cpp
	g++ $(FLAGS) -c $< -o $@

launchsim: rocket.o main.o
	g++ $(FLAGS) $^ -o $@
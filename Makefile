FLAGS = -pedantic -Wall -Wextra -std=c++14

launchsim: main.cpp
	g++ -g $(FLAGS) main.cpp -o launchsim 
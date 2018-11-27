Fuzzy c-Means Library
===
[![Build Status](https://semaphoreci.com/api/v1/ahmad88me/fcm-cpp/branches/master/badge.svg)](https://semaphoreci.com/ahmad88me/fcm-cpp)
[![codecov](https://codecov.io/gh/ahmad88me/fcm-cpp/branch/master/graph/badge.svg)](https://codecov.io/gh/ahmad88me/fcm-cpp)

This project implements Fuzzy c-means from the original paper by James Bezdek. 


## Core advantages of this implementation

* The use of Eign library for matrices instead of fix n-dimentional arrays which are common in prototyping fuzzy c-means library.
* The use of dynamic sizes for matrices.
* Make the FCM as a class, with all the related variables inside, so no need for Global variables and it also allow having multiple instances of FCM.
* Add automated tests with googletest.
* Add coverage for the tests.

## To Install
`make install`

## To run the tests
`make test`

## To run the tests with Docker
`sh scripts/run_tests_with_docker.sh`

## Dependancies

If you are using ubuntu, there is a helper script 
`sh scripts/setup.sh`
*But we recommend using docker* as this script might not be well maintained

Some of the dependancies:
* eigen3
* git
* cmake
* googletests
* g++
* lcov
* curl


#! /bin/bash
./configure --with-blaslib=veclib --with-lapacklib=veclib LDFLAGS=-L/opt/local/lib/ CXXFLAGS=-I/opt/local/include/ CXX=g++ --prefix=$PROJECT_TDA

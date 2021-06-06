#!/bin/bash

BINDIR := ../bin
DEPDIR := .d
$(shell mkdir -p $(DEPDIR) >/dev/null)
$(shell mkdir -p $(BINDIR) >/dev/null)
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td
POSTCOMPILE = [ -f $(DEPDIR)/$*.d ] && diff -q $(DEPDIR)/$*.Td $(DEPDIR)/$*.d || mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d

INCLUDE := -I../include

SRC := $(wildcard *.cpp) ## all cpp files
CXXFLAGS = -g $(INCLUDE) -pipe -Wall -Wno-sign-compare -Wno-unused-variable -Wno-unused-but-set-variable -std=c++14 -O3
LDFLAGS = 

OS = $(shell uname -s)
ifeq ($(OS), Linux)
  CXX = g++
  CXXFLAGS += -fopenmp
  LDFLAGS += -fopenmp
endif
ifeq ($(OS), Darwin)
  CXX = g++
endif

LDFLAGS += -ldl

.DEFAULT_GOAL := simulate

%.o : %.cpp
%.o : %.cpp $(DEPDIR)/%.d
	$(CXX) $(DEPFLAGS) $(CXXFLAGS) -c $< -o $@
	$(POSTCOMPILE)

$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d

simulate : simulate.o system.o body.o math_func.o customize.o
	$(CXX) -o $(BINDIR)/$@ $^ $(LDFLAGS)

.PHONY: clean
clean:
	-rm *.o $(BINDIR)/simulate

-include $(patsubst %,$(DEPDIR)/%.d,$(basename $(SRC)))

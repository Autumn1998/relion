#atom.cpp  atom.h  main.cpp  Makefile  mpicc.sh  mrcfile_atom2.cpp  mrcfile_atom.h  run_ASART1.sh

RM := rm -rf


LIBS :=

O_SRCS := 
CPP_SRCS := 
C_UPPER_SRCS := 
C_SRCS := 
S_UPPER_SRCS := 
OBJ_SRCS := 
CU_SRCS := 
ASM_SRCS := 
CXX_SRCS := 
C++_SRCS := 
CC_SRCS := 
CU_DEPS := 
OBJS := 
C++_DEPS := 
C_DEPS := 
CC_DEPS := 
CPP_DEPS := 
EXECUTABLES := 
CXX_DEPS := 
C_UPPER_DEPS := 


CPP_SRCS += \
main.cpp \
mrcfile_atom2.cpp \
atom.cpp \



OBJS += \
main.o \
mrcfile_atom2.o \
atom.o \



# Tool invocations
TXBRSIRT: $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: MPIC++ Linker'		
	mpic++  $(OBJS)  -lm   -lrt -lpthread -o TXBRSIRT
	@echo 'Finished building target: $@'
	@echo ' '

%.o: %.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: MPIC++ Compiler'
	mpic++ -O0 -g3 -Wall -c -fmessage-length=0 -lm  -o $@ $<    
	@echo 'Finished building: $<'
	@echo ' '




# Other Targets
clean:
	-$(RM) $(CU_DEPS)$(OBJS)$(C_DEPS)$(CC_DEPS)$(CPP_DEPS)$(EXECUTABLES)$(CXX_DEPS)$(C_UPPER_DEPS) TXBRSIRT
	-@echo ' '



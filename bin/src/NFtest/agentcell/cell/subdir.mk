################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/NFtest/agentcell/cell/cell.cpp \
../src/NFtest/agentcell/cell/environment.cpp \
../src/NFtest/agentcell/cell/util.cpp 

OBJS += \
./src/NFtest/agentcell/cell/cell.o \
./src/NFtest/agentcell/cell/environment.o \
./src/NFtest/agentcell/cell/util.o 

CPP_DEPS += \
./src/NFtest/agentcell/cell/cell.d \
./src/NFtest/agentcell/cell/environment.d \
./src/NFtest/agentcell/cell/util.d 


# Each subdirectory must supply rules for building sources it contributes
src/NFtest/agentcell/cell/%.o: ../src/NFtest/agentcell/cell/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/include/c++/4.8.5 -O3 -c -fmessage-length=0 -std=c++11 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



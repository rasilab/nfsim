################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/NFreactions/reactantLists/reactantList.cpp \
../src/NFreactions/reactantLists/reactantTree.cpp 

OBJS += \
./src/NFreactions/reactantLists/reactantList.o \
./src/NFreactions/reactantLists/reactantTree.o 

CPP_DEPS += \
./src/NFreactions/reactantLists/reactantList.d \
./src/NFreactions/reactantLists/reactantTree.d 


# Each subdirectory must supply rules for building sources it contributes
src/NFreactions/reactantLists/%.o: ../src/NFreactions/reactantLists/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/include/c++/4.8.5 -O3 -c -fmessage-length=0 -std=c++11 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



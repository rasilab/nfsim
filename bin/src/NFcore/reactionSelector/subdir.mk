################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/NFcore/reactionSelector/directSelector.cpp \
../src/NFcore/reactionSelector/logClassSelector.cpp 

OBJS += \
./src/NFcore/reactionSelector/directSelector.o \
./src/NFcore/reactionSelector/logClassSelector.o 

CPP_DEPS += \
./src/NFcore/reactionSelector/directSelector.d \
./src/NFcore/reactionSelector/logClassSelector.d 


# Each subdirectory must supply rules for building sources it contributes
src/NFcore/reactionSelector/%.o: ../src/NFcore/reactionSelector/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/include/c++/4.8.5 -O3 -c -fmessage-length=0 -std=c++11 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



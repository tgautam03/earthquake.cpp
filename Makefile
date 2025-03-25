CC = g++ -std=c++11
HOST_COMPILE_FLAG = -c
CPU_OPTIMIZE = -O3 -Ofast -march=native -funroll-loops -ffast-math -flto
#######################
##### Dependencies ####
#######################
# Data Structures
MatrixFP32.o: src/MatrixFP32.cpp
	@$(CC) $(HOST_COMPILE_FLAG) $(CPU_OPTIMIZE) src/MatrixFP32.cpp -o build/MatrixFP32.o

Tensor3dFP32.o: src/Tensor3dFP32.cpp
	@$(CC) $(HOST_COMPILE_FLAG) $(CPU_OPTIMIZE) src/Tensor3dFP32.cpp -o build/Tensor3dFP32.o

Tensor4dFP32.o: src/Tensor4dFP32.cpp
	@$(CC) $(HOST_COMPILE_FLAG) $(CPU_OPTIMIZE) src/Tensor4dFP32.cpp -o build/Tensor4dFP32.o

# Utils
utils.o: src/utils.cpp
	@$(CC) $(HOST_COMPILE_FLAG) $(CPU_OPTIMIZE) src/utils.cpp -o build/utils.o

#######################
###### Wave Sim #######
#######################
1d_sim: MatrixFP32.o utils.o src/Wave1D.cpp sim_1d.cpp 
	@$(CC) $(CPU_OPTIMIZE) build/MatrixFP32.o build/utils.o src/Wave1D.cpp sim_1d.cpp -o bin/sim_1d.out

3d_sim: Tensor4dFP32.o Tensor3dFP32.o utils.o src/Wave3D.cpp sim_3d.cpp 
	@$(CC) $(CPU_OPTIMIZE) build/Tensor4dFP32.o build/Tensor3dFP32.o build/utils.o src/Wave3D.cpp sim_3d.cpp -o bin/sim_3d.out

3d_sim_slices: Tensor4dFP32.o Tensor3dFP32.o utils.o src/Wave3D.cpp sim_3d_slices.cpp 
	@$(CC) $(CPU_OPTIMIZE) build/Tensor4dFP32.o build/Tensor3dFP32.o build/utils.o src/Wave3D.cpp sim_3d_slices.cpp -o bin/sim_3d_slices.out

##########################
# Clean executable files #
##########################
clean: 
	@echo "Removing object files..."
	rm bin/*.out build/*.o

clean_img:
	@echo "Removing frames..."
	rm viz/frames/*.png
FC = gfortran

NC_INC = -I/usr/lib64/gfortran/modules
NC_LIB = -lnetcdff

SRC_DIR = Src
BIN_DIR = Bin

MOD_FLAG = -J$(SRC_DIR)

TARGET = $(BIN_DIR)/run_model.exe

SRCS = $(SRC_DIR)/mod_io.f90 $(SRC_DIR)/main.f90 $(SRC_DIR)/mod_advect.f90 $(SRC_DIR)/mod_physics.f90
OBJS = $(SRCS:.f90=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	@mkdir -p $(BIN_DIR)
	$(FC) -o $@ $^ $(NC_LIB)

$(SRC_DIR)/main.o: $(SRC_DIR)/mod_io.o $(SRC_DIR)/mod_advect.o $(SRC_DIR)/mod_physics.o

%.o: %.f90
	$(FC) -c $< -o $@ $(NC_INC) $(MOD_FLAG)

clean:
	rm -f $(SRC_DIR)/*.o $(SRC_DIR)/*.mod $(TARGET)

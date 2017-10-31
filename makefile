# simple make file
SOURCES=PES_GP_Symm.f90 
PRODUCT=3Ar_PES.out


all: $(PRODUCT)

$(PRODUCT) : $(SOURCES)
	gfortran -o $(PRODUCT) $(SOURCES)

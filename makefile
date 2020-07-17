FC	= ifort
F_OPT	= -mcmodel=medium -shared-intel -fpic

.SUFFIXES: .f90 .o

OBJ_F	= marker_contourOut13_analyze.o bmp_out.o #tif_out.o

TARGET_F	= marker


$(TARGET_F):$(OBJ_F)
	$(FC) $(F_OPT) -o $(TARGET_F) $(OBJ_F)

# .f90.o:
# 	$(FC) -c free=f90 $(F_OPT) $<

.f90.o:
	$(FC) -c $(F_OPT) $<

.f.o:
	$(FC) -c $(F_OPT) $<

clean:
	rm -f $(TARGET_F) *.o *.mod

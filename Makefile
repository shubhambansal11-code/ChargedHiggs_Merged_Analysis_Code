SDIR  = ./
IDIR  = ./
LDIR  = ./
ROOTSYS   = $(shell echo $ROOTSYS)
MAKEFLAGS = --no-print-directory -r -s
INCLUDE = $(shell root-config --cflags)
LIBS    = $(shell root-config --libs) -lTMVA -lMLP -lTreePlayer -lMinuit


BINS7 = Event-Loop
OBJS = 
all: $(BINS7)


$(BINS7): % : main_batch.C
	@echo -n "Building $@ ... "
	$(CXX) $(CCFLAGS) $< -I$(IDIR) $(INCLUDE) TH1Fs/TH1Fs.C main/EventLoop.C utilis/NeutrinoBuilder.C utilis/Chi2_minimization.C -o execute_V2.exe $(LIBS)
	@echo "Done"

clean:
	rm -f $(BINS) 

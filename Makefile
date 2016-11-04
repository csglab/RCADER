MAKE = make		#change this line if you are using a different GNU make software

dirFASTAtoRF = ./src/FASTAtoRF
dirRCADER = ./src/RCADER
dirExamples= ./examples

all: Extract MK_dir CC_FASTAtoRF CC_RCADER RM_objectFiles

Extract:
	gunzip $(dirExamples)/*/*.gz

MK_dir:
	mkdir -p ./bin

CC_FASTAtoRF: $(dirFASTAtoRF)/Makefile
	$(MAKE) -C $(dirFASTAtoRF)
	
CC_RCADER: $(dirRCADER)/Makefile
	$(MAKE) -C $(dirRCADER)

RM_objectFiles:
	rm -f $(dirFASTAtoRF)/*.o $(dirRCADER)/*.o 

clean:
	rm -r -f $(dirFASTAtoRF)/*.o $(dirRCADER)/*.o ./bin

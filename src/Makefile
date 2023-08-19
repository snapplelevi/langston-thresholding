# Copy of Makefile.win, but this file should target Unix-like machines and build
# the dependencies for the package


EXTERNALDIR = $(PACKAGEDIR)/external
BUILDDIR = $(PACKAGEDIR)/build
IGRAPHDIR = igraph-0.8.0
IGRAPHLIBS = igraph-libs
INCLUDEDIR = $(PACKAGEDIR)/include
LIBDIR = $(PACKAGEDIR)/lib

all: igraph_extract 

.PHONY: igraph_extract
igraph_extract:
# Extract igraph files from the tar file
	tar -zxvf ../external/$(IGRAPHDIR).tar.gz -C ../external
# Extract linked libraries from the tar file   
	tar -zxvf ../lib/$(IGRAPHLIBS).tar.gz -C ../lib/
	echo "IS THIS THE DEFAULT BEING RUN???"

.PHONY: igraph
igraph:
	mkdir -p $(BUILDDIR)
	mkdir -p $(INCLUDEDIR)
	cd $(EXTERNALDIR) && \
	tar -zxvf $(IGRAPHDIR).tar.gz && \
	cd $(IGRAPHDIR) && \
	./configure --libdir=$(LIBDIR) --prefix=$(BUILDDIR) && \
	$(MAKE) MAKEFLAGS= && \
	$(MAKE) install && \
	cp $(BUILDDIR)/include/igraph/* $(INCLUDEDIR)/igraph/

.PHONY: clean
clean:
		rm ./src/*.o
		rm ./src/thresholding.dll
		rm -r ./external/$(IGRAPHDIR)
PKG_CPPFLAGS = -I../external/igraph-0.8.0/include -I../external/alglib
PKG_LIBS = -L../lib -ligraph-0.8.0

EXTERNALDIR = $(PWD)/external
BUILDDIR = $(PWD)/build
IGRAPHDIR = igraph-0.8.0
INCLUDEDIR = $(PWD)/include
LIBDIR = $(PWD)/lib

all: igraph

.PHONY: igraph
igraph:
	mkdir -p $(BUILDDIR)
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
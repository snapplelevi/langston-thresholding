

BUILDIR = $(PWD)/build
IGRAPHDIR = igraph-0.8.0
INCLUDEDIR = $(PWD)/include
LIBDIR = $(PWD)/lib

all: igraph

.PHONY: igraph
igraph:
	echo "TEST"
	@mkdir -p $(BUILDDIR)
	cd $(EXTERNALDIR) && \
	tar zxvf $(IGRAPHDIR).tar.gz && \
	cd $(IGRAPHDIR) && \
	./configure --libdir=$(LIBDIR) --prefix=$(BUILDDIR) && \
	$(MAKE) MAKEFLAGS= && \
	$(MAKE) install && \
	cp $(BUILDDIR)/include/igraph/* $(INCLUDEDIR)/igraph/

.PHONY: clean
clean:
		rm ./src/*.o
		rm ./src/thresholding.dll
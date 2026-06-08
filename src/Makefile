# src/Makefile
# R CMD INSTALL includes: Makeconf, [makevars_site], THIS_FILE, [makevars_user]
# so $(SHLIB_EXT), $(ALL_CPPFLAGS), $(ALL_CXXFLAGS), $(ALL_LIBS), $(SHLIB_LINK)
# are all defined by Makeconf before this file is parsed.

PKG_CPPFLAGS = -DTMBAD_FRAMEWORK
PKG_CXXFLAGS = $(SHLIB_OPENMP_CXXFLAGS)
PKG_LIBS     = $(SHLIB_OPENMP_CXXFLAGS)

SHLIB_LDFLAGS = $(SHLIB_CXXLDFLAGS)

CXX_LINK = $(SHLIB_CXXLD) $(SHLIB_LDFLAGS) $(LIBR0) $(LDFLAGS)

all: gllvm$(SHLIB_EXT) gllvm_HO$(SHLIB_EXT)

gllvm$(SHLIB_EXT): gllvm.o
	$(CXX_LINK) -o $@ $< $(ALL_LIBS)

gllvm_HO$(SHLIB_EXT): gllvm_HO.o
	$(CXX_LINK) -o $@ $< $(ALL_LIBS)

shlib-clean:
	@rm -f gllvm.o gllvm_HO.o

.PHONY: all shlib-clean

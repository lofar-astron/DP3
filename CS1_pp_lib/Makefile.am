SUBDIRS=src test include

ACLOCAL_AMFLAGS = -I autoconf_share

pkgextdir     = $(prefix)/config/$(PACKAGE)
pkgext_DATA   = pkgext pkgextcppflags pkgextcxxflags pkgextldflags

DISTCHECK_CONFIGURE_FLAGS= \
      --with-common=$(prefix) \
      --with-aips++ 

EXTRA_DIST = \
      Makefile.common \
      CS1_IDPPP.spec \
      autoconf_share/compiletool

include $(top_srcdir)/Makefile.common

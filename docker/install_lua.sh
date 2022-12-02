#!/bin/bash

# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Script to install lua
set -euo pipefail

pushd /tmp

echo "Downloading Lua"
curl -fsSLO https://www.lua.org/ftp/lua-${LUA_VERSION}.tar.gz
tar zxf lua-${LUA_VERSION}.tar.gz

echo "Building & installing Lua"
pushd lua-${LUA_VERSION}

# Modify Makefiles to build shared library, as explained in
# https://stackoverflow.com/questions/20848275/compiling-lua-create-so-files
patch -p1 <<'EOF'
diff -Naur lua-5.3.6.orig/Makefile lua-5.3.6/Makefile
--- lua-5.3.6.orig/Makefile	2022-11-17 11:19:00.069298984 +0100
+++ lua-5.3.6/Makefile	2022-11-17 11:32:10.759562052 +0100
@@ -6,6 +6,10 @@
 # Your platform. See PLATS for possible values.
 PLAT= none

+# Lua version and release.
+V= 5.3
+R= $V.6
+
 # Where to install. The installation starts in the src and doc directories,
 # so take care if INSTALL_TOP is not an absolute path. See the local target.
 # You may want to make INSTALL_LMOD and INSTALL_CMOD consistent with
@@ -41,18 +45,14 @@
 # What to install.
 TO_BIN= lua luac
 TO_INC= lua.h luaconf.h lualib.h lauxlib.h lua.hpp
-TO_LIB= liblua.a
+TO_LIB= liblua.a liblua.so.$(R)
 TO_MAN= lua.1 luac.1

-# Lua version and release.
-V= 5.3
-R= $V.6
-
 # Targets start here.
 all:	$(PLAT)

 $(PLATS) clean:
-	cd src && $(MAKE) $@
+	cd src && $(MAKE) $@ V=$(V) R=$(R)

 test:	dummy
 	src/lua -v
@@ -63,6 +63,8 @@
 	cd src && $(INSTALL_DATA) $(TO_INC) $(INSTALL_INC)
 	cd src && $(INSTALL_DATA) $(TO_LIB) $(INSTALL_LIB)
 	cd doc && $(INSTALL_DATA) $(TO_MAN) $(INSTALL_MAN)
+	ln -s $(INSTALL_LIB)/liblua.so.$(R) $(INSTALL_LIB)/liblua.so.$(V)
+	ln -s $(INSTALL_LIB)/liblua.so.$(R) $(INSTALL_LIB)/liblua.so

 uninstall:
 	cd src && cd $(INSTALL_BIN) && $(RM) $(TO_BIN)
diff -Naur lua-5.3.6.orig/src/Makefile lua-5.3.6/src/Makefile
--- lua-5.3.6.orig/src/Makefile	2022-11-17 11:19:00.118298380 +0100
+++ lua-5.3.6/src/Makefile	2022-11-17 11:30:42.499648924 +0100
@@ -7,7 +7,7 @@
 PLAT= none

 CC= gcc -std=gnu99
-CFLAGS= -O2 -Wall -Wextra -DLUA_COMPAT_5_2 $(SYSCFLAGS) $(MYCFLAGS)
+CFLAGS= -O2 -Wall -Wextra -DLUA_COMPAT_5_2 $(SYSCFLAGS) $(MYCFLAGS) -fPIC
 LDFLAGS= $(SYSLDFLAGS) $(MYLDFLAGS)
 LIBS= -lm $(SYSLIBS) $(MYLIBS)

@@ -29,6 +29,7 @@
 PLATS= aix bsd c89 freebsd generic linux macosx mingw posix solaris

 LUA_A=	liblua.a
+LUA_SO= liblua.so
 CORE_O=	lapi.o lcode.o lctype.o ldebug.o ldo.o ldump.o lfunc.o lgc.o llex.o \
 	lmem.o lobject.o lopcodes.o lparser.o lstate.o lstring.o ltable.o \
 	ltm.o lundump.o lvm.o lzio.o
@@ -43,7 +44,7 @@
 LUAC_O=	luac.o

 ALL_O= $(BASE_O) $(LUA_O) $(LUAC_O)
-ALL_T= $(LUA_A) $(LUA_T) $(LUAC_T)
+ALL_T= $(LUA_A) $(LUA_T) $(LUAC_T) $(LUA_SO)
 ALL_A= $(LUA_A)

 # Targets start here.
@@ -59,6 +60,11 @@
 	$(AR) $@ $(BASE_O)
 	$(RANLIB) $@

+$(LUA_SO): $(CORE_O) $(LIB_O)
+	$(CC) -shared -ldl -Wl,-soname,$(LUA_SO).$(V) -o $@.$(R) $? -lm $(MYLDFLAGS)
+	ln -sf $(LUA_SO).$(R) $(LUA_SO).$(V)
+	ln -sf $(LUA_SO).$(R) $(LUA_SO)
+
 $(LUA_T): $(LUA_O) $(LUA_A)
 	$(CC) -o $@ $(LDFLAGS) $(LUA_O) $(LUA_A) $(LIBS)
EOF

make linux -j${THREADS}
make install

popd

# Clean up the files from the build
rm -rf lua-${LUA_VERSION} lua-${LUA_VERSION}.tar.gz

popd
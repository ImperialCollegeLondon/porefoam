#------------------------------------------------------------------------------
# =========                 |
# \\      /  F ield         | foam-extend: Open Source CFD
#  \\    /   O peration     | Version:     4.1
#   \\  /    A nd           | Web:         http://www.foam-extend.org
#    \\/     M anipulation  | For copyright notice see file Copyright
#------------------------------------------------------------------------------
# License
#     This file is part of foam-extend.
#
#     foam-extend is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation, either version 3 of the License, or (at your
#     option) any later version.
#
#     foam-extend is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.
#
# Script
#     Makefile
#
# Description
#     Makefile used to compile wmake utilities.
#
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# The Makefile use a POSIX shell
#------------------------------------------------------------------------------

SHELL      = /bin/sh

#------------------------------------------------------------------------------
# Initialise suffices list
#------------------------------------------------------------------------------

.SUFFIXES:
.SUFFIXES: .o


#------------------------------------------------------------------------------
# set compilation and dependency building rules
#------------------------------------------------------------------------------

GENERAL_RULES = $(WM_DIR)/rules/General
RULES         = $(WM_DIR)/rules/$(WM_ARCH)$(WM_COMPILER)
BIN           = $(WM_DIR)/bin/$(WM_ARCH)$(WM_COMPILER)

include $(RULES)/general
include $(RULES)/$(WM_LINK_LANGUAGE)


#------------------------------------------------------------------------------
# targets
#------------------------------------------------------------------------------

all:  $(BIN)/dirToString $(BIN)/wmkdep

clean:
	rm -f $(BIN)/dirToString $(BIN)/wmkdep 2>/dev/null


$(BIN)/dirToString: dirToString.c
	@mkdir -p $(BIN)
	$(cc) $(cFLAGS) dirToString.c -o $(BIN)/dirToString

$(BIN)/wmkdep: wmkdep.l
	@mkdir -p $(BIN)
	flex wmkdep.l || printf "please install flex for this to work:\n  sudo apt install flex # in Ubuntu\n\n"
	$(cc) $(cFLAGS) lex.yy.c -o $(BIN)/wmkdep  || printf "flex not found, ignore wmkdep"
	@rm -f lex.yy.c 2>/dev/null


#------------------------------------------------------------------------------

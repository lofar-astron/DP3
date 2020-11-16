// Register.h: Register AOFlag steps in DPPP
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
/// @file
/// @brief Register AOFlag steps in DPPP
/// @author Ger van Diepen

#ifndef DPPP_AOFLAGGER_REGISTER_H
#define DPPP_AOFLAGGER_REGISTER_H

/// Define the function (without name mangling) to register the 'constructor'.
extern "C" {
void register_aoflagger();
}

#endif

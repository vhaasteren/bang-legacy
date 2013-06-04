/* filefunctions.h -- Helpful functions to work with file-access
 *
 * Rutger van Haasteren 19 March 2008 haasteren@strw.leidenuniv.nl
 *
 * Copyright (C) 2006-2008 Rutger van Haasteren.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 * */

#include <sys/types.h>
#include <dirent.h>
#include <vector>
#include <string>
#include <stdio.h>

#ifndef __FILEFUNCTIONS_H__
#define __FILEFUNCTIONS_H__

using namespace std;

int get_allfiles_dir(string dir, vector<string> &files);
int get_basenames_dir(string dir, string extension, vector<string> &names);
bool FileExists(string strFilename);

// __FILEFUNCTIONS_H__
#endif

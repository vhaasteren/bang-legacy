/* filefunctions.cpp -- Helpful functions to work with file-access
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

#include "filefunctions.h"
#include <algorithm>
#include <sys/stat.h> 

/* This function reads all files is a directory, and stores the filenames in an
 * array.
 *
 * Input: dir: The directory name
 *
 * Output: files: Array of strings filled with the filenames
 */
int get_allfiles_dir(string dir, vector<string> &files) {
    DIR *dp;
    struct dirent *dirp;
    if((dp = opendir(dir.c_str())) == NULL) {
        return 1;
    }

    while ((dirp = readdir(dp)) != NULL) {
        files.push_back(string(dirp->d_name));
    }
    closedir(dp);
    return 0;
} // get_allfiles_dir

/* This function returns all filenames of a certain extension in a directory,
 * and returns the base-names (without the extension) in an array.
 *
 * Input: dir: The directory name
 *
 * Input: extension: The extension
 *
 * Output: names: Array of strings filled with the file basenames
 */
int get_basenames_dir(string dir, string extension, vector<string> &names) {
  vector<string> files = vector<string>();
  string basename;

  get_allfiles_dir(dir, files);

  for(int i=0; i<files.size(); i++) {
    if( (int) files[i].rfind('.') >= 0) {
      if( files[i].compare(files[i].rfind('.'), files[i].length() - files[i].rfind('.'), extension) == 0) {
	// The extension matches

	basename = files[i].substr(0, files[i].rfind('.'));
	names.push_back(basename);

      } // if files
    } // if files
  } // for i

  // Sort the names alphabetically using qsort
  sort(names.begin(), names.end());
  return 0;
} // get_basenames_dir

/* This function returns true if the file in question exists. If the file
 * does not exist false is returned
 * */
bool FileExists(string strFilename) {
  struct stat stFileInfo;
  bool blnReturn;
  int intStat;

  // Attempt to get the file attributes
  intStat = stat(strFilename.c_str(),&stFileInfo);
  if(intStat == 0) {
    // We were able to get the file attributes
    // so the file obviously exists.
    blnReturn = true;
  } else {
    // We were not able to get the file attributes.
    // This may mean that we don't have permission to
    // access the folder which contains this file. If you
    // need to do that level of checking, lookup the
    // return values of stat which will give you
    // more details on why stat failed.
    blnReturn = false;
  }
  return(blnReturn);
}

/*
int main()
{
    string dir = string(".");
    vector<string> basenames = vector<string>();

    get_basenames_dir(dir, string(".par"), basenames);

    printf("Size of basenames: %i\n", basenames.size());

    for (int i=0; i<basenames.size(); i++) {
      printf("basename[%i]: %s\n", i, basenames[i].c_str());
    }
    return 0;
}
*/

// stacker, Python module for stacking of interferometric data.
// Copyright (C) 2014  Lukas Lindroos
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.  //
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "DataIO.h"
#include "DataIOFits.h"
#include "Chunk.h"

using std::string;

const double pi = M_PI;

DataIOFits::DataIOFits(const char* infilename, const char* outfilename,
    pthread_mutex_t* mutex)
{
  this->mutex = *mutex;

  // Initiate status variables
  fitsioStatus = 0, anynull = 0;

  fits_open_file(&infile, infilename, READONLY, &fitsioStatus);
  if(fitsioStatus != 0)
  {
    string message;
    message += "Failed to open infile: ";
    message += infilename;

//     fileException e(fileException::OPEN, message);
    throw fileException(fileException::OPEN, message);
//     throw e;
  }

  if(outfilename != NULL && strlen(outfilename) > 0)
  {
    fits_open_file(&outfile, outfilename, READWRITE, &fitsioStatus);
    if(fitsioStatus != 0)
    {
      string message;
      message += "Failed to open outfile: ";
      message += outfilename;

      throw fileException( fileException::OPEN, message);
    }
  }
  else
  {
    outfile = NULL;
  }

  string object;
  try {
    object = readKeywordStr("OBJECT");
  }
  catch(keynameError e){
    throw fileException(fileException::HEADER_INFO_MISSING, "Could not find OBJECT keyword");
  }

  if(object == "MULTI")
    throw fileException(fileException::MOSAIC, "Mosaic uvfits files not supported");

  nfields = 1;
  x_phase_centre = new float[nfields];
  y_phase_centre = new float[nfields];

  float ra, dec;

  // TODO: Follow aips memo 102, read CRVAL first.
  try {
    if(nfields == 1) { 
      ra = readKeywordFloat("OBSRA");
      dec = readKeywordFloat("OBSDEC");
      x_phase_centre[0] = ra/180.*pi;
      y_phase_centre[0] = ra/180.*pi;
    } 
  }
  catch(keynameError e){
    throw fileException(fileException::HEADER_INFO_MISSING, 
                        "No field information present in fits file.");
  }

  std::cout << "(RA, DEC): (" << ra << ", " << dec << ")" << std::endl;
}

float DataIOFits::readKeywordFloat(const char* keyword) {

  char keyname[FLEN_KEYWORD];
  char comment[FLEN_COMMENT];

  float value;
  int status;

  strcpy(keyname, keyword);

  fits_read_key_flt(infile, keyname, &value, comment, &status);
  if(status != 0) {
    throw keynameError();
  }

  return value;
}

long DataIOFits::readKeywordLong(const char* keyword) {

  char keyname[FLEN_KEYWORD];
  char comment[FLEN_COMMENT];
  long value;
  int status;

  strcpy(keyname, keyword);

  fits_read_key_lng(infile, keyname, &value, comment, &status);
  if(status != 0) {
    throw keynameError();
  }

  return value;
}

string DataIOFits::readKeywordStr(const char* keyword) {

  char keyname[FLEN_KEYWORD];
  char comment[FLEN_COMMENT];
  char value[FLEN_VALUE];
  int status;

  strcpy(keyname, keyword);

  fits_read_key_str(infile, keyname, value, comment, &status);
  if(status != 0) {
    throw keynameError();
  }

  return string(value);
}


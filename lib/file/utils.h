/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by J-Donald Tournier, 23/07/09.

    This file is part of MRtrix.

    MRtrix is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MRtrix is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MRtrix.  If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef __file_ops_h__
#define __file_ops_h__

#include <stdint.h>
#include <sys/types.h>
#include <sys/stat.h>
//#include <sys/sendfile.h>
#include <fcntl.h>

#include "debug.h"
#include "mrtrix.h"
#include "types.h"
#include "file/path.h"
#include "file/overwrite.h"

#define TMPFILE_ROOT "mrtrix-tmp-"
#define TMPFILE_ROOT_LEN 11


namespace MR
{
  namespace File
  {

    namespace
    {
      inline char random_char ()
      {
        char c = rand () % 62;
        if (c < 10) return c+48;
        if (c < 36) return c+55;
        return c+61;
      }

      inline int do_create (const std::string& filename, int64_t size)
      {
        int fid = open(filename.c_str(), O_CREAT | O_RDWR | O_TRUNC, 0644);
        if (fid < 0) 
          throw Exception ("error creating file \"" + filename + "\": " + strerror (errno));

        if (size) size = ftruncate(fid, size);
        if (size) 
          throw Exception ("WARNING: cannot resize file \"" + filename + "\": " + strerror (errno));

        return fid;
      }
    }




    inline void create (ConfirmOverwrite& confirm_overwrite, const std::string& filename, int64_t size = 0)
    {
      confirm_overwrite (filename);
      close (do_create (filename, size));
    }

    inline void create (const std::string& filename, int64_t size = 0)
    { 
      ConfirmOverwrite::single_file (filename);
      close (do_create (filename, size));
    }



    inline void resize (const std::string& filename, int64_t size)
    {
      debug ("resizing file \"" + filename + "\" to " + str (size) + "...");

      int fd = open(filename.c_str(), O_RDWR, 0644);
      if (fd < 0)
        throw Exception ("error opening file \"" + filename + "\" for resizing: " + strerror (errno));
      int status = ftruncate(fd, size);
      close (fd);
      if (status)
        throw Exception ("cannot resize file \"" + filename + "\": " + strerror (errno));
    }




    inline bool is_tempfile (const std::string& name, const char* suffix = NULL)
    {
      if (Path::basename (name).compare (0, TMPFILE_ROOT_LEN, TMPFILE_ROOT)) return false;
      if (suffix) if (!Path::has_suffix (name, suffix)) return false;
      return true;
    }




    inline std::string create_tempfile (int64_t size = 0, const char* suffix = NULL)
    {
      debug ("creating temporary file of size " + str (size));

      std::string filename (TMPFILE_ROOT"XXXXXX.");
      if (suffix) filename += suffix;

      int fid;
      do {
        for (int n = 0; n < 6; n++)
          filename[TMPFILE_ROOT_LEN+n] = random_char();
      }
      while ( (fid = open(filename.c_str(), O_CREAT | O_RDWR | O_EXCL, 0644)) < 0);

      int status = size ? ftruncate(fid, size) : 0;
      close (fid);
      if (status) throw Exception ("cannot resize file \"" + filename + "\": " + strerror (errno));
      return filename;
    }


    inline void mkdir (const std::string& folder) 
    {
      if (::mkdir (folder.c_str(), 0777))
        throw Exception ("error creating folder \"" + folder + "\": " + strerror (errno));
    }
    

    inline void unlink (const std::string& file) 
    {
      if (::unlink (file.c_str()))
        throw Exception ("error deleting file \"" + file + "\": " + strerror (errno));;
    }

    inline void rmdir (const std::string& folder, bool recursive = false)
    {
      if (recursive) {
        Path::Dir dir (folder);
        std::string entry;
        while ((entry = dir.read_name()).size()) {
          std::string path = Path::join (folder, entry);
          if (Path::is_dir (path))
            rmdir (path, true);
          else 
            unlink (path);
        }
      }
      debug ("deleting folder \"" + folder + "\"...");
      if (::rmdir (folder.c_str()))
        throw Exception ("error deleting folder \"" + folder + "\": " + strerror (errno));
    }


  }
}

#endif


//# BlobStreamable.h: Interface for classes that can be streamed using blobs.
//#
//# Copyright (C) 2006
//# ASTRON (Netherlands Institute for Radio Astronomy)
//# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//#
//# This file is part of the LOFAR software suite.
//# The LOFAR software suite is free software: you can redistribute it and/or
//# modify it under the terms of the GNU General Public License as published
//# by the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The LOFAR software suite is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//#
//# $Id: BlobStreamable.h 26359 2013-09-04 13:58:06Z loose $

#ifndef LOFAR_BLOB_BLOBSTREAMABLE_H
#define LOFAR_BLOB_BLOBSTREAMABLE_H

// \file
// Interface for classes that can be streamed using blobs.

//# Includes
#include "../Common/Singleton.h"
#include "../Common/ObjectFactory.h"

namespace DP3
{
  //# Forward Declarations.
  class BlobIStream;
  class BlobOStream;

  // \addtogroup Blob
  // @{

  // Classes that want to be blob-streamable must implement the methods
  // declared in this interface.
  class BlobStreamable
  {
  public:
    // Destructor.
    virtual ~BlobStreamable() {}

    // Create a new BlobStreamable object by deserializing the contents of
    // the blob input stream \a bis. The first element in the input stream
    // \a bis shall be a string containing the class type of the object to
    // be constructed.
    static BlobStreamable* deserialize(BlobIStream& bis);

    // Serialize \c *this by writing its contents into the blob output
    // stream \a bos.
    void serialize(BlobOStream& bos) const;

  protected:
    // Return the class type of \c *this as a string.
    virtual const std::string& classType() const = 0;

  private:
    // Read the contents from the blob input stream \a bis into \c *this.
    virtual void read(BlobIStream& bis) = 0;

    // Write the contents of \c *this into the blob output stream \a bos.
    virtual void write(BlobOStream& bos) const = 0;

    // The private methods must be accessible for DH_BlobStreamable. 
    friend class DH_BlobStreamable;
  };

  // Factory that can be used to generate new BlobStreamable objects.
  // The factory is defined as a singleton.
  typedef Singleton< ObjectFactory< BlobStreamable*(), std::string > >
  BlobStreamableFactory;

  // @}
    
} // namespace LOFAR

#endif

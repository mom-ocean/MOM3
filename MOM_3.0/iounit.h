c====================== include file "iounit.h" ========================
c
c     i/o units and related variables
c
c     taum1disk = disk pointer for tau-1 latitude rows
c     taudisk   = disk pointer for tau   latitude rows
c     taup1disk = disk pointer for tau+1 latitude rows
c     kflds     = disk unit used for two dimensional fields
c     latdisk   = disk units for latitude rows (alternately pointed to 
c                by taum1disk, taudisk, and taup1disk)
c              
c     iodoc  = unit for documentation
c     iostab = unit for stability testing
c     iotim  = unit for time means
c     iotim1 = scratch disk (SSD) unit for accumulating time means
c     ionew1 = unit for reading sponge layer data
c     ionew2 = mirror unit of sponge layer data on SSD
#ifndef orlanski
# if defined obc_south || defined obc_north
c     ionew3 = unit for reading obc data (T,S) 
c     ionew4 = mirror unit of obc data (T,S) on SSD
c     ionew7 = unit for reading obc data (psi)
c     ionew8 = mirror unit of obc data (psi) on SSD
# endif
# if defined obc_west || defined obc_east
c     ionew5 = unit for reading obc data (T,S) 
c     ionew6 = mirror unit of obc data (T,S) on SSD
c     ionew9 = unit for reading obc data (psi)
c     ionew10= mirror unit of obc data (psi) on SSD
# endif
#endif
c     
c     for the following, a control # < 0 implies that unformatted data
c     will be written to a unit selected by the i/o manager "iomngr.F"
c     and given a hardwired name (grep getunit *.F to see names) 
c     and formatted data (to stdout) will be written. if a # > 0 and
c      # <> stdout, only unformatted data will be written.
c
c     iotavg = control # for tracer averages
c     iotmb  = control # for writing tracer meridional budget.
c     iotrmb = control # for term balances for tracer and momentum
c     ioglen = control # for writing global energetics integrals
c     iovmsf = control # for writing meridional stream function
c     iogyre = control # for writing gyre transport.
c     ioprxz = control # for writing x-z sections from latitudes
c     ioext  = control # for writing external mode (stream function)
c     iodsp  = control # for writing diagnostic surface pressure
c     iotsi  = control # for writing time step integrals
c     ioxbt  = control # for writing time averaged xbt data
c     iozmbc = control # for writing zonal mean surf boundary conditions
c
      integer taum1disk, taudisk, taup1disk
      common /iounit/ taum1disk, taudisk, taup1disk
      common /iounit/ latdisk(3), kflds
      common /iounit/ iodoc, iostab, iotavg, iotmb, iotrmb
      common /iounit/ iotim, iotim1
      common /iounit/ ioglen, iovmsf, iogyre, ioprxz, ioext, iodsp
      common /iounit/ iotsi, iozmbc, ionew1, ionew2, ioxbt
#ifndef orlanski
# if defined obc_south || defined obc_north
      common /iounit/ ionew3, ionew4, ionew7, ionew8
# endif
# if defined obc_west || defined obc_east
      common /iounit/ ionew5, ionew6, ionew9, ionew10
# endif
#endif
c
c     iotext  = 80 character text string for describing the details
c                of the next unformatted data record.
c     expnam  = 60 character text string for the experiment name
c
      character iotext*80, expnam*60
      common /iounitc/ iotext, expnam
c
c     when writing unformatted data records in MOM, each data record is
c     preceded by a header record which was written as:
c     write(unit) stamp, iotext, expnam
c     where stamp is a 32 character specification of the model date &
c     time corresponding to the time step when the data was written and
c     iotext is a 80 character description of what is in the
c     data record and how it is to be read. expnam is a 60 character
c     experiment name which shows which experiment wrote the data.
c     this makes it easy to decipher any unformatted output from the
c     model by using a program similar to the following: 
c
c      program decifr
c
c-----------------------------------------------------------------------
c      decifer an unformatted file from MOM by showing the header
c      records. the file needs to copied to file "fort.21"
c-----------------------------------------------------------------------
c
c      character*32 stamp
c      character*80 iotext
c      character*60 expnam
c
c      iounit = 21
c      rewind iounit
c      do n=1,100000
c
c        read the header record
c
c        read (iounit, end=110) stamp, iotext, expnam
c        write (*,'(1x,a32,1x,a80)') stamp, iotext
c
c        skip the data record
c
c        read (iounit)
c      enddo
c110   continue
c      write (*,*) " => end of file on fort.",iounit
c      stop
c      end
c
c     note: all unformatted diagnostic MOM data is handled this way.
c     to insure that data is read properly, verify that arrays are
c     dimensioned correctly by comparing the listed variables against
c     those in the *.h files. (grep -i -n "variable" *.h) Also, most
c     data from MOM is written IEEE 32bit so it is read directly by
c     most workstations. However, when trying to read these IEEE files
c     on the CRAY, they must be assigned IEEE before being read.
c     Some diagnostic data is averaged over time before being written.
c     In these cases, the time "stamp" refers to the last time step
c     at the end of the averaging period. An averaging interval is
c     also written as part of the data. Averaging periods = zero 
c     indicate instantaneous data.
c

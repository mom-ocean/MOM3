module read_ic_mod
  parameter (maxlatt = 256) ! maximum length of an attribute string

  type, public  	      ::  nc_variable
   sequence
   character(len=128)	      :: name
   integer		      :: ixtype
   integer		      :: ndim
   integer		      :: dimids(4)  ! assume four grid dimensions
   integer		      :: len   (4)  ! assume four grid dimensions
   character(len=80)	      :: units
  end type nc_variable

  type(nc_variable), allocatable, public, dimension(:) :: icvar

  integer, public, dimension(4)  :: start,count  ! assume four grid dimensions
  integer, public			  :: istatus
  integer, public			  :: ncid,ndims,nvars
  integer, public			  :: ixtype
  integer , allocatable 		  :: len_dim(:)
  character(len=128), allocatable	  :: name_dim(:)
  character(len=80)             	  :: time_name,time_unit
  character(len=14), dimension(6), public :: xtype_name
  data xtype_name /'NF_BYTE','NF_CHAR','NF_SHORT','NF_INT',&
                       'NF_FLOAT','NF_DOUBLE'/
! data buffers
#if defined _CRAY
  real(kind=8), allocatable, public	  :: data4 (:,:,:,:)
#else
  real(kind=4), allocatable, public	  :: data4 (:,:,:,:)
#endif
  real(kind=8), allocatable, public	  :: data8 (:,:,:,:)
!
! it is assumed that MOM always uses missing=-1e34 in netCDF in/output
! using -1e33 here to avoid problems with single and double precision data
  real, public :: missing
  data missing /-1.0e33/
!
#ifdef read_ic 
!     the following code is used for input of ic data in setocn.F
!-----------------------------------------------------------------------
!
!     read initial model fields from a netCDF file prepared by user 
!     to keep things easy the following assumptions are made:
!     - input ic_file covers i=2...imt-1, j=2...jmt-1 like snapshots
!       (ic input may be read from a snapshot)
!     - 3dim data are defined on a 4dim grid x,y,z,t (2dim on x,y,t),
!       hence 4dim buffers and access vectors (data4, data8, start, count)
!     - the x,y,z-axes are equal to model grid (size_check imt,jmt,km)
!     - ic data units are equal to model units (ic units shown in results)
!     - ic data time is specified by the time step counter ictimestep
!       (ic data time is shown but NOT checked against model ic_time)
!     - ic_file contains all necessary 3dim and/or 2dim data in any order
!       (if required data for tracers are not present the program stops, 
!       if surface elevation and/or velocity data are absent, these
!       fields are initialized with zero)
!     - netCDF variable names must be equal to trname(nt) for tracers,
!       'u','v' for 3dim velocities, 'eta', 'surf_press' for sea surface
!       elevation or pressure, and 'ubar','vbar' for barotropic velocities
!     - input data may have single or double precision
!     - missing=-1e34 is used for blanking land in netCDF files
!	    
!-----------------------------------------------------------------------
  integer, public :: ic3dim,ic2dim,ictimestep 
  integer, allocatable, public :: id3dim(:),id2dim(:)
! default ic input file name and time step (user input by namelist.icdata) 
  character(len=40), public :: ic_file
  data ic_file /'ic.dta.nc'/
  data ictimestep /1/	    

  contains

      subroutine open_ic
  
#include "netcdf.inc"  
#include "param.h"
!!#include "docnam.h" !! has not the Fortran90 format (! instead of c)
      character*12 trname
      common /docnam/ trname(nt)
      
      integer ioic
      namelist /icdata/ ic_file,ictimestep ! specify user input there

      call getunit (ioic, 'namelist.icdata','formatted sequential rewind')
      read  (ioic, icdata)
      write (stdout,icdata)
      call relunit (ioic)
      call getunit (ioic, 'document.dta','formatted sequential append')
      write (ioic,icdata)
      call relunit (ioic)
!
!-----------------------------------------------------------------------
!     open the netCDF ic_file 
!-----------------------------------------------------------------------
!
      write(stdout,*) 'Open File ',trim(ic_file)
      istatus=NF_OPEN(ic_file,0,ncid)
      if (istatus.ne.nf_noerr) then
        write(stderr,*) 'Error => failed opening ic file ',ic_file
        write(stderr,*) nf_strerror(istatus)
        call abort()
      endif
      write(stdout,*) 'Reading initial conditions from file ',trim(ic_file) &
        	     ,' opened with istatus = ',istatus,' ncid is ',ncid
!
!-----------------------------------------------------------------------
!     get the number of dimensions and variables
!-----------------------------------------------------------------------
!
      istatus=NF_INQ_NDIMS(ncid, ndims)
      if (istatus.ne.nf_noerr) then
        write(stderr,*) 'Error => failed reading '
        write(stderr,*) nf_strerror(istatus)
        call abort()
      endif
      istatus=NF_INQ_NVARS(ncid, nvars)
      if (istatus.ne.nf_noerr) then
        write(stderr,*) 'Error => failed reading '
        write(stderr,*) nf_strerror(istatus)
        call abort()
      endif
      write(stdout,*) 'ic file ',trim(ic_file) &
           ,' has',ndims,' dimensions for',nvars,' variables'
!
!-----------------------------------------------------------------------
!     get information about netCDF dimensions (4 axes x,y,z,t)
!-----------------------------------------------------------------------
!
       allocate(name_dim(ndims),len_dim(ndims))
       do id=1,ndims
         istatus=NF_INQ_DIM(ncid, id, name_dim(id), len_dim(id))
# ifdef debug_setocn
         write(stdout,*) 'Dimension(',id,'): ', trim(name_dim(id)) &
	                ,' Length: ', len_dim(id)
# endif
       enddo
!
!-----------------------------------------------------------------------
!     check with model dimensions
!-----------------------------------------------------------------------
!
      call size_check (len_dim(1)+2, len_dim(2)+2, len_dim(3) &
        ,'open_ic (read_ic)', 'stop')
      call size_check (len_dim(4)+2, len_dim(5)+2, len_dim(6) &
        ,'open_ic (read_ic)', 'stop')
!           
!-----------------------------------------------------------------------
!     get names and types of netCDF variables
!-----------------------------------------------------------------------
!
      allocate(icvar(1:nvars))
      do id=1,nvars
        istatus=NF_INQ_varname(ncid,id,icvar(id)%name)
        istatus=NF_INQ_vartype(ncid,id,icvar(id)%ixtype)
	icvar(id)%units=' '
        istatus=NF_GET_ATT_text(ncid,id,'units',icvar(id)%units)
# ifdef debug_setocn
        write(stdout,*) 'Variable(',id,'): ', trim(icvar(id)%name) &
        	      ,' ',xtype_name(icvar(id)%ixtype) &
        	      ,' ',trim(icvar(id)%units)
# endif
      enddo
!     
!-----------------------------------------------------------------------
!     get number, id, and length of dimensions
!-----------------------------------------------------------------------
!
      do id=ndims+1,nvars
     	icvar(id)%len(:) = 1
     	istatus=NF_INQ_varndims(ncid,id,icvar(id)%ndim)
     	istatus=NF_INQ_vardimid(ncid,id,icvar(id)%dimids)
     	do n=1,icvar(id)%ndim
     	  icvar(id)%len(n) = len_dim(icvar(id)%dimids(n))
     	enddo
# ifdef debug_setocn
     	write(stdout,*) 'Variable(',id,'): ',trim(icvar(id)%name) &
     	  ,' number of dimensions :', icvar(id)%ndim
     	write(stdout,*) 'Variable(',id,'): ',trim(icvar(id)%name) &
     	  ,' dimension id :',(icvar(id)%dimids(n),n=1,icvar(id)%ndim)
     	write(stdout,*) 'Variable(',id,'): ',trim(icvar(id)%name) &
     	  ,' dimensions: ',(icvar(id)%len(n),n=1,4)
# endif
      enddo
!
!-----------------------------------------------------------------------
!     assign variable id to input arrays, variable names must fit exactly
!-----------------------------------------------------------------------
!
      allocate(id3dim(nt+2))
      do n=1,nt
     	id3dim(n)=-1
     	do id=ndims+1,nvars
     	  if (trim(trname(n)).eq.trim(icvar(id)%name)) id3dim(n)=id
     	enddo
     	if (id3dim(n).lt.1) then
     	  write(stderr,*) 'Error => no ic input found for tracer ' &
               ,trim(trname(n)),' in file ',trim(ic_file)
     	  call abort()
     	else
     	  write(stdout,*) 'get tracer ',trname(n) &
     	       ,' from variable(',id3dim(n),'): ',trim(icvar(id3dim(n))%name) &
     	       ,', units: ',trim(icvar(id3dim(n))%units)
     	endif
      enddo
      id3dim(nt+1)=-1
      id3dim(nt+2)=-1
      do id=ndims+1,nvars
     	if (trim(icvar(id)%name).eq.'u') id3dim(nt+1)=id
     	if (trim(icvar(id)%name).eq.'v') id3dim(nt+2)=id
      enddo
      if (id3dim(nt+1).lt.1 .or. id3dim(nt+2).lt.1) then
     	write(stdout,*) 'Warning => no ic velocities initialize (u,v) with zero'
     	ic3dim=nt
      else
     	ic3dim=nt+2
     	write(stdout,*) 'get ic uhor from variable(',id3dim(nt+1),'): ' &
     	     ,trim(icvar(id3dim(nt+1))%name)
     	write(stdout,*) 'get ic vhor from variable(',id3dim(nt+2),'): ' &
     	     ,trim(icvar(id3dim(nt+2))%name)
      endif
# ifdef free_surf_on
!
      allocate(id2dim(nt+2))
      id2dim(:)=-1
      do id=ndims+1,nvars
#  ifdef explicit_free_surface
     	if (trim(icvar(id)%name).eq.'eta') id2dim(1)=id
#  else
     	if (trim(icvar(id)%name).eq.'surf_press') id2dim(1)=id
#  endif
     	if (trim(icvar(id)%name).eq.'ubar') id2dim(2)=id
     	if (trim(icvar(id)%name).eq.'vbar') id2dim(3)=id
      enddo
      if (id2dim(1).lt.1) then
        ic2dim=0
        write(stdout,*) 'Warning => no ic surface data, initialise with zero'
        else
        ic2dim=1
        write(stdout,*) 'get surface data from variable(',id2dim(1),'): ' &
             ,trim(icvar(id2dim(1))%name) &
	     ,', units: ',trim(icvar(id2dim(1))%units)
	if (id2dim(2).lt.1 .or. id2dim(3).lt.1) then
	  write(stdout,*) 'Warning => no ic velocities,' &
               ,' initialize (ubar,vbar) with zero'
        else
	  ic2dim=ic2dim+2
          write(stdout,*) 'get ic ubar from variable(',id2dim(2),'): ' &
               ,trim(icvar(id2dim(2))%name) &
	     ,', units: ',trim(icvar(id2dim(2))%units)
          write(stdout,*) 'get ic vbar from variable(',id2dim(3),'): ' &
               ,trim(icvar(id2dim(3))%name) &
	     ,', units: ',trim(icvar(id2dim(3))%units)
        endif
      endif
# endif
!
!-----------------------------------------------------------------------
!     allocate the data buffers for slabs (imt,km)
!-----------------------------------------------------------------------
!
      isreal4=0
      isreal8=0
      do id=ndims+1,nvars
     	if (icvar(id)%ixtype.eq.NF_FLOAT) then
     	  isreal4=isreal4+1
     	elseif (icvar(id)%ixtype.eq.NF_DOUBLE) then
     	  isreal8=isreal8+1
     	else
     	  write(stdout,*) '=> ic data must be float or double' &
     	       ,' found ',icvar(id)%ixtype,' for ',icvar(id)%name
     	  call abort()
     	endif		  
      enddo
      if (isreal4.gt.0) then
        allocate(data4(2:imtm1,1,km,1))
     	write(stdout,*) 'ic real4 data buffer allocated'
      endif
      if (isreal8.gt.0) then
        allocate(data8(2:imtm1,1,km,1))
     	write(stdout,*) 'ic real8 data buffer allocated'
      endif
!
!-----------------------------------------------------------------------
!     get time of input data
!-----------------------------------------------------------------------
!
      if (ictimestep.gt.len_dim(idtime)) then
      endif
      idtime=-1
      do id=1,nvars
     	if (trim(icvar(id)%name).eq.'Time') idtime=id
      enddo
      if (idtime.gt.0) then
        time_name=' '
	time_unit=' '
        istatus=NF_GET_ATT_text(ncid,idtime,'long_name',time_name)
        istatus=NF_GET_ATT_text(ncid,idtime,'units',time_unit)
        if (icvar(idtime)%ixtype.eq.NF_FLOAT) then
          istatus=NF_GET_VAR1_real(ncid,idtime,ictimestep,timeic)
        else
          istatus=NF_GET_VAR1_double(ncid,idtime,ictimestep,timeic)
        endif
        write(stdout,*) 'ic data time name: ',trim(time_name)
        write(stdout,*) 'ic data time is',timeic,' ',trim(time_unit)
        if (ictimestep.gt.len_dim(idtime)) then
     	  write(stdout,*) '=> ictimestep',ictimestep &
	       ,' exeeds number of timesteps',len_dim(idtime) &
	       ,' in input file ',trim(ic_file),' see namelist.icdata'
	  call abort()
        endif
      else
     	write(stdout,*) '=> no ic time found in ',trim(ic_file)
	ictimestep=1
      endif   
!
      return
      end subroutine open_ic

      subroutine close_ic
!      
#include "stdunits.h"
!      
      istatus=NF_CLOSE(ncid)
      if (allocated(data4)) deallocate(data4)
      if (allocated(data8)) deallocate(data8)
      deallocate(icvar,name_dim,len_dim,id3dim)
      if (allocated(id2dim)) deallocate(id2dim)
# ifdef debug_setocn
      write(stdout,*) 'File ',trim(ic_file) &
           ,' closed with istatus = ',istatus
# endif
      return
      end subroutine close_ic
#endif /* read_ic */
end module read_ic_mod 

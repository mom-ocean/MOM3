      write(stdout,*) " doing topog user inputs for test_case_A:"

c     Close off Panama
      kmt(92,46) = 0
      kmt(93,45) = 0
      kmt(94,45) = 0
      kmt(94,44) = 0
      kmt(95,44) = 0
c     Deepen ocean point in Drake Passage
      kmt(99,13) = 15
c     Sink Iceland and deepen surrounding points
c     as per jrt analysis
      kmt(115,65) = 12
      kmt(116,65) = 14
      kmt(117,65) = 14
      kmt(118,65) = 17
      kmt(119,65) = 17
      kmt(120,65) = 18
c
      kmt(113,64) = 12
      kmt(114,64) = 13
      kmt(115,64) = 13
      kmt(116,64) = 14
      kmt(117,64) = 15
      kmt(118,64) = 15
      kmt(119,64) = 16
c
      kmt(109,63) =  8
      kmt(110,63) = 14
      kmt(111,63) = 15
      kmt(112,63) = 15
      kmt(113,63) = 14
      kmt(114,63) = 13
      kmt(115,63) = 13
      kmt(116,63) = 12
      kmt(117,63) = 12
      kmt(118,63) = 12
      kmt(119,63) = 15
c
      kmt(113,62) = 16
      kmt(118,62) = 14
      kmt(119,62) = 14
      kmt(120,62) = 13
      kmt(121,62) = 13
      kmt(122,62) = 11
      kmt(  1,62) = 13
      kmt(  2,62) = 11
      kmt(  3,62) =  8
c
      kmt(119,61) = 13
      kmt(120,61) = 13
      kmt(121,61) = 10
      kmt(122,61) =  8
      kmt(  1,61) = 10
      kmt(  2,61) =  8
c
      kmt(116,60) = 15
      kmt(117,60) = 15
      kmt(119,60) = 13

c     Add the Bering Strait
      do j=63,64
         do i=64,65
           kmt(i,j) = 3
         enddo
      enddo

c     Deepen Spitzbergen
      kmt(7,68) = 11
      kmt(8,68) = 11
c     Add points to New Zealand
      kmt(59,17) = 0
      kmt(59,18) = 0
      kmt(60,18) = 0
c     Delete Great Lakes
      do j=60,64
         do i=92,96
            if( kmt(i,j) .ne.0 ) then
               kmt(i,j) = 0
            endif
         enddo
      enddo
c     Deepen Florida Straits
      kmt(94,51) = 0
      do j=49,52
         do i=94,95
            if( kmt(i,j) .gt. 0 .and. kmt(i,j).lt.12) then
               write(stdout,*) 
     *              ' topog ',i,j,' deepen from',kmt(i,j),' to 12'
               kmt(i,j) = 12
            endif
         enddo
      enddo
c     Open Channel thru Baffin Bay
*Eli: close it and make this an inland-sea mixing straits:
c     kmt(96,69) = 7
c     kmt(97,69) = 7
c     kmt(98,69) = 7
c     kmt(97,68) = 7
c     kmt(98,68) = 7
c     Deepen one point in Indonesian Straights
      kmt(45,35) = 10
c     Deepen Gibraltar for cross mixing
      kmt(119,54 ) = 10
      kmt(119,53 ) = 10

c     Keep japan sea as separate from Pacific ocean, connected by diffusion only:
c     add a point to Japan:
      kmt(47,53)=0
c     eliminate two isolated cells south-east of Korea:
      kmt(43,53)=0
      kmt(42,53)=0
c     eliminate a land point that sticks westward out of Indonesia(?):
      kmt(41,42)=18
c     sink two islands south-east of japan (Taiwan?):
      kmt(43,44) = 20
      kmt(42,47) = 21
c     prevent isolated point near sea of japan:
      kmt(44,53) = 13
c     Enlarge hudson bay and connect it to Labrador Sea as inland sea:
      do i=92,95
         kmt(i,62) = 3
      end do
c     Make sure connection to lab sea is to a velocity point by streightening
c     the coast line there:
      do i=99,100
         kmt(i,62) = 0
      end do
      do i=92,95
         kmt(i,61) = 3
      end do
      kmt(94,60) = 3
      kmt(95,60) = 3
      kmt(94,63) = 3
      kmt(95,63) = 3
      kmt(94,64) = 3
      kmt(95,64) = 3

c     Straighten the coastline of antarctica:
      do i=23,28
         kmt(i,9) = 0
      end do
      do i=8,10
         kmt(i,8) = 0
      end do

c     Make baltic an inland sea:
      kmt(4,60) = 0

c     Deepen points north west of baltic to prevent sharp topographic
c     gradients and digging there:
      do i=2,4
         do j=60,64
            if (kmt(i,j).ne.0.and.kmt(i,j).lt.13) kmt(i,j)=13
         end do
      end do
      do i=121,121
         do j=60,64
            if (kmt(i,j).ne.0.and.kmt(i,j).lt.13) kmt(i,j)=13
         end do
      end do

c     disconnect Madagaskar from Africa:
      kmt(15,26) = 16
      kmt(16,26) = 16

c     Enlarge Caspean Sea:
      kmt(18,55) = 3
      kmt(19,54) = 3

c     Enlarge black sea and connect to Mediterranean as inland sea:
      kmt(12,57) = 14

c     Enlarge Red Sea and connect as inland sea to indian ocean:
      kmt(15,48) = 9
      kmt(14,50) = 9
      kmt(14,48) = 9
      kmt(17,46) = 9
      kmt(18,45) = 14

c     Fill the Arabian Sea:
      kmt(19,50) = 0

c     Add a couple of points to Antarctica in the Drake Passage:
      kmt(99,9) = 0
      kmt(100,9) = 0
c     fill a diffusion-only row of single T-cells south of australia
c     on the antarctica coastline:
      kmt(57,6) = 0
      kmt(58,6) = 0

c     connect some "islands" in the indonesian throughflow to the larger land
c     masses there, leaving only a single passage:
      do j=33,35
         kmt(41,j) = 0
      end do
      kmt(37,31) = 0
      kmt(37,34) = 0
c     Eliminate some single point dips in the indonesian coastline:
      kmt(36,38) = 0
      kmt(48,26) = 0
      kmt(48,27) = 0
      kmt(42,34) = 0

c     sink spiesbergen:
      kmt(20,66) = 4

c     sink another arctic "island":
      kmt(83,66) = 3

c     eliminate single dips in arctic ocean coastline:
      kmt(17,64) = 0
      kmt(41,66) = 0

c     straighten Egyptian coastline in Mediterranean:
      kmt(10,52) = 0
c     another East Mediterranean fix to avoind an isolated cell:
      kmt(13,52) = 13

c     Without the following change the east med is connected with a
c     diffusion point only, and the free surface there is not changing
c     with the west med:
      kmt(10,52) = 14
      kmt(9,52)  = 14

c     All topography shallower than level 3 deepened to level 4
      do j=1,jmt
         do i=1,imt
            if( kmt(i,j) .le. 2 .and. kmt(i,j).ne.0 ) then
               kmt(i,j) = 3
            endif
         enddo
      enddo
c     Option fill_isolated_cells is overridden below using
c     -Dfill_isolated_dont_change_landmask in case it tries to change
c     the landmask (that is, the land outline).  Sometime this is
c     desirable, when it tries to eliminate marginal seas and such.  But
c     sometimes the overriding is not desirable, when
c     fill_isolated_cells tries to eliminate ocean points that are
c     isolated and not essential.  So, after running the model, check
c     the printout file to see exactly where fill_isolated_cells was
c     overridden, and in cases where this was not a good idea make the
c     changes that fill_isolated cells would try to do:
      kmt(3,59)=0
      kmt(48,60)=0
      kmt(49,60)=0
      kmt(35,45)=0
      kmt(37,36)=0
      kmt(37,37)=0
      kmt(37,38)=0
      kmt(48,28)=0
      kmt(49,28)=0
      kmt(91,66)=0
      kmt(92,66)=0
      kmt(96,66)=0
      kmt(98,62)=0
      kmt(104,65)=0
C     In the med, add the following T point, to prevent an isolated
C     cell near the straits of Gibraltar:
      kmt(121,55)=8

# ifdef bbl_mask
c     mask bbl in the Bering straits:
      mask_t_bbl(64:65,63) = 0
# endif

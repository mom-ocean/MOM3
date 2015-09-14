#if defined rs6000
c       This file is needed for AIX or other UNIX versions, if the 
c       entry points in libs of the functions specified below have  
c       C underscores. In this case the linker quits with unresolved 
c       references.  
c       Enable -Drs6000 and compile this file together with the other 
c       files of MOM.
c       Martin Schmidt    e-mail martin.schmidt@io-warnemuende.de       
c
        subroutine flush(lunit)
        integer*4 lunit
        call flush_(lunit)
        return
        end

#else
        subroutine libpatch
        return
        end
#endif

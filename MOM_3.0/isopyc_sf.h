c======================== include file "isopyc_sf.h" ===================
c
c     statement functions
c
      drodxe(i,k,j,ip) =    drho(i+ip,k,j,1)*ddxt(i,k,j,1) + 
     &                      drho(i+ip,k,j,2)*ddxt(i,k,j,2) 
      drodze(i,k,j,ip,kr) = drho(i+ip,k,j,1)*ddzt(i+ip,k-1+kr,j,1) + 
     &                      drho(i+ip,k,j,2)*ddzt(i+ip,k-1+kr,j,2)
c
      drodyn(i,k,j,jq) =    drho(i,k,j+jq,1)*ddyt(i,k,j,1) + 
     &                      drho(i,k,j+jq,2)*ddyt(i,k,j,2) 
      drodzn(i,k,j,jq,kr) = drho(i,k,j+jq,1)*ddzt(i,k-1+kr,j+jq,1) + 
     &                      drho(i,k,j+jq,2)*ddzt(i,k-1+kr,j+jq,2)
c
      drodxb(i,k,j,ip,kr) = drho(i,k+kr,j,1)*ddxt(i-1+ip,k+kr,j,1) + 
     &                      drho(i,k+kr,j,2)*ddxt(i-1+ip,k+kr,j,2) 
      drodyb(i,k,j,jq,kr) = drho(i,k+kr,j,1)*ddyt(i,k+kr,j-1+jq,1) + 
     &                      drho(i,k+kr,j,2)*ddyt(i,k+kr,j-1+jq,2) 
      drodzb(i,k,j,kr) =    drho(i,k+kr,j,1)*ddzt(i,k,j,1) + 
     &                      drho(i,k+kr,j,2)*ddzt(i,k,j,2)
c
# ifdef full_tensor
      drodye(i,k,j,ip,jq) = drho(i+ip,k,j,1)*ddyt(i+ip,k,j-1+jq,1) + 
     &                      drho(i+ip,k,j,2)*ddyt(i+ip,k,j-1+jq,2) 
      drodxn(i,k,j,ip,jq) = drho(i,k,j+jq,1)*ddxt(i-1+ip,k,j+jq,1) + 
     &                      drho(i,k,j+jq,2)*ddxt(i-1+ip,k,j+jq,2) 
# endif

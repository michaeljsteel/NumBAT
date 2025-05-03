#include "numbat_decl.h"

!     This subroutine takes an array of complex numbers, computes their squared magnitudes,
!     and sorts the indices of these magnitudes in descending order using a hybrid quicksort-insertion
!     sort algorithm.
!     The sorted indices are stored in the indx array, which can be used to
!     reorder the original array cor for other purposes where sorted order is required.

SUBROUTINE find_eigvals_order(n, arr, indx, order, nberr)

   use alloc
   use numbatmod

   type(NBError) nberr

   integer(8) n, indx(n), order
   complex(8) arr(n)

   integer(8), parameter :: M=7
   integer(8), parameter :: NSTACK=50

   integer(8) i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
   double precision a
   double precision, dimension(:), allocatable :: arr_0
   character(len=EMSG_LENGTH) :: emsg


!write(*,*) 'feo', n, arr
   call double_nalloc_1d(arr_0, n, 'arr0', nberr); RET_ON_NBERR(nberr)


   do  j=1,n
      indx(j) = j
   enddo

   if (order .eq. 0) then
      do j=1,n   ! sort largest real part first
         arr_0(j) = -1*dble(arr(j)**2)
      enddo

   else if (order .eq. 1) then
      do j=1,n   ! sort smallest  first
         arr_0(j) = abs(arr(j))
      enddo

   else
      write(emsg, *) 'Unknown order in find_eigvals_order'
      call nberr%set(NBERR_UNKNOWN_SORT_ORDER, emsg)
      return
   end if


   jstack=0
   l=1
   ir=n
1  if(ir-l.lt.M)then
      do 13 j=l+1,ir
         indxt=indx(j)
         a = arr_0(indxt)
         do 12 i=j-1,l,-1
            if( arr_0(indx(i)) .le. a)goto 2
            indx(i+1)=indx(i)
12       continue
         i=l-1
2        indx(i+1)=indxt
13    continue

      if(jstack.eq.0)return

      ir=istack(jstack)
      l=istack(jstack-1)
      jstack=jstack-2
   else
      k=(l+ir)/2
      itemp=indx(k)
      indx(k)=indx(l+1)
      indx(l+1)=itemp

      if( arr_0(indx(l)) .gt. arr_0(indx(ir))) then
         itemp=indx(l)
         indx(l)=indx(ir)
         indx(ir)=itemp
      endif

      if( arr_0(indx(l+1)) .gt. arr_0(indx(ir))) then
         itemp=indx(l+1)
         indx(l+1)=indx(ir)
         indx(ir)=itemp
      endif

      if( arr_0(indx(l)) .gt. arr_0(indx(l+1))) then
         itemp=indx(l)
         indx(l)=indx(l+1)
         indx(l+1)=itemp
      endif

      i=l+1
      j=ir
      indxt=indx(l+1)
      a = arr_0(indxt)

3     continue
      i=i+1
      if( arr_0(indx(i)) .lt. a) goto 3

4     continue
      j=j-1
      if( arr_0(indx(j)) .gt. a) goto 4

      if(j.lt.i) goto 5

      itemp=indx(i)
      indx(i)=indx(j)
      indx(j)=itemp
      goto 3

5     indx(l+1)=indx(j)
      indx(j)=indxt
      jstack=jstack+2

      if(jstack.gt.NSTACK) then
         write(emsg,*) 'NSTACK too small in indexx'
         call nberr%set(NBERR_SORT_STACK_TOO_SMALL, emsg)
         return
      endif

      if(ir-i+1.ge.j-l)then
         istack(jstack)=ir
         istack(jstack-1)=i
         ir=j-1
      else
         istack(jstack)=j-1
         istack(jstack-1)=l
         l=i
      endif

   endif
   goto 1

end

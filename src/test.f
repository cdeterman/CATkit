

      subroutine test(matrix, M, N)
      
      IMPLICIT NONE
      
      integer, intent(in)  ::  M,N
      double precision, intent(inout)   ::  matrix(M,N)
      integer   ::  i,j
      
      print *, "hello world"
      
      do i=1,M
        do j=1,N
            print *, matrix(i,j)
        end do
      end do
      
      do 10 i=1,M
        do 20 j=1,N
            print *, i
            print *, j
            matrix(i,j) = 2.0D0
   20   end do
   10 end do
      
      do i=1,M
        do j=1,N
            print *, matrix(i,j)
        end do
      end do
      
      end
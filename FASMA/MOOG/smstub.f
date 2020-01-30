      subroutine sm_window(i1,i2,i3,i4,i5,i6)
      integer i1,i2,i3,i4,i5,i6
      end subroutine

      subroutine sm_box(i1,i2,i3,i4)
      integer i1,i2,i3,i4
      end subroutine

      subroutine sm_location(i1,i2,i3,i4)
      integer i1,i2,i3,i4
      end subroutine


      subroutine sm_ltype(i1)
      integer i1
      end subroutine

      subroutine sm_grid(i1,i2)
      integer i1,i2
      end subroutine


      subroutine sm_putlabel(i1,c1)
      integer i1
      character c1(:)
      end subroutine

      subroutine sm_label(c1)
      character*80 c1
      end subroutine

      subroutine sm_defvar(c1,c2)
      character*80 c1,c2
      end subroutine

      subroutine sm_ctype(c1)
      character*80 c1
      end subroutine


      subroutine sm_xlabel(c1)
      character*80 c1
      end subroutine

      subroutine sm_ylabel(c1)
      character*80 c1
      end subroutine

      subroutine sm_device(c1)
      character*80 c1
      end subroutine


      subroutine sm_limits(x1,x2,x3,x4)
      real x1,x2,x3,x4
      end subroutine

      subroutine sm_ticksize(x1,x2,x3,x4)
      real x1,x2,x3,x4
      end subroutine

      subroutine sm_draw(x1,x2)
      real x1,x2
      end subroutine

      subroutine sm_lweight(x1)
      real x1
      end subroutine

      subroutine sm_expand(x1)
      real x1
      end subroutine

      subroutine sm_relocate(x1,x2)
      real x1,x2
      end subroutine

      subroutine sm_curs(x1,x2,i1)
      real x1,x2
      integer i1
      end subroutine


      subroutine sm_angle(x1)
      real x1
      end subroutine

      subroutine sm_ptype(x1,i1)
      integer i1
      real*4 x1(:)
      end subroutine

      subroutine sm_points(x1,x2,i1)
      integer i1
      real*4 x1(:), x2(:)
      end subroutine

      subroutine sm_connect(x1,x2,i1)
      integer i1
      real*4 x1(:), x2(:)
      end subroutine

      subroutine sm_histogram(x1,x2,i1)
      integer i1
      real*4 x1(:), x2(:)
      end subroutine
      
      subroutine sm_graphics()
      end subroutine

      subroutine sm_gflush()
      end subroutine

      subroutine sm_erase()
      end subroutine

      subroutine sm_alpha()
      end subroutine

      subroutine sm_hardcopy()
      end subroutine

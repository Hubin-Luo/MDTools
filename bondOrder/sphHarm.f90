        FUNCTION sphHarm(l,m,theta,phi)
!=====================================================
!spherical harmonics up to l=10 using the table from
!wikipedia.
!theta and phi are in rad.
!                                          Hubin Luo
!=====================================================
        IMPLICIT NONE
!
        COMPLEX(8) :: sphHarm
        INTEGER, INTENT(IN) :: l,m
        REAL(8), INTENT(IN) :: theta, phi
        REAL(8) :: pi,sin2,sin3,sin4,sin5,sin6,sin7,sin8,sin9,sin10
        REAL(8) :: cos2,cos3,cos4,cos5,cos6,cos7,cos8,cos9,cos10
        COMPLEX(8),PARAMETER :: im=(0.d0, 1.d0)
!
        IF(l<0 .OR. ABS(m)>l) THEN
           WRITE(*,*) "Wrong integers for spherical harmonics!!"
           RETURN
        END IF
!
        pi=ACOS(-1.d0)
        sin2=SIN(theta)*SIN(theta)
        sin3=SIN(theta)*SIN(theta)*SIN(theta)
        sin4=SIN(theta)*SIN(theta)*SIN(theta)*SIN(theta)
        sin5=SIN(theta)*SIN(theta)*SIN(theta)*SIN(theta)*SIN(theta)
        sin6=SIN(theta)*SIN(theta)*SIN(theta)*SIN(theta)*SIN(theta)*SIN(theta)
        sin7=SIN(theta)*SIN(theta)*SIN(theta)*SIN(theta)*SIN(theta)*SIN(theta)*SIN(theta)
        sin8=SIN(theta)*SIN(theta)*SIN(theta)*SIN(theta)*SIN(theta)*SIN(theta)*SIN(theta)*SIN(theta)
        sin9=SIN(theta)*SIN(theta)*SIN(theta)*SIN(theta)*SIN(theta)*SIN(theta)*SIN(theta)*SIN(theta)*SIN(theta)
        sin10=SIN(theta)*SIN(theta)*SIN(theta)*SIN(theta)*SIN(theta)*SIN(theta)*SIN(theta)*SIN(theta)*SIN(theta)*SIN(theta)
        cos2=COS(theta)*COS(theta)
        cos3=COS(theta)*COS(theta)*COS(theta)
        cos4=COS(theta)*COS(theta)*COS(theta)*COS(theta)
        cos5=COS(theta)*COS(theta)*COS(theta)*COS(theta)*COS(theta)
        cos6=COS(theta)*COS(theta)*COS(theta)*COS(theta)*COS(theta)*COS(theta)
        cos7=COS(theta)*COS(theta)*COS(theta)*COS(theta)*COS(theta)*COS(theta)*COS(theta)
        cos8=COS(theta)*COS(theta)*COS(theta)*COS(theta)*COS(theta)*COS(theta)*COS(theta)*COS(theta)
        cos9=COS(theta)*COS(theta)*COS(theta)*COS(theta)*COS(theta)*COS(theta)*COS(theta)*COS(theta)*COS(theta)
        cos10=COS(theta)*COS(theta)*COS(theta)*COS(theta)*COS(theta)*COS(theta)*COS(theta)*COS(theta)*COS(theta)*COS(theta)
!
        SELECT CASE(l)
          CASE(0)
            sphHarm=0.5*SQRT(1.d0/pi)
            RETURN
          CASE(1)
            SELECT CASE(m)
              CASE(-1)
                sphHarm=0.5*SQRT(1.5/pi)*EXP(-1.d0*im*phi)*SIN(theta)
                RETURN
              CASE(0)
                sphHarm=0.5*SQRT(3.d0/pi)*COS(theta)
                RETURN
              CASE(1)
                sphHarm=-0.5*SQRT(1.5/pi)*EXP(im*phi)*SIN(theta)
                RETURN
            END SELECT
          CASE(2)
            SELECT CASE(m)
              CASE(-2)
                sphHarm=0.25*SQRT(7.5/pi)*EXP(-2.d0*im*phi)*sin2
                RETURN
              CASE(-1)
                sphHarm=0.5*SQRT(7.5/pi)*EXP(-1.d0*im*phi)*SIN(theta)*COS(theta)
                RETURN
              CASE(0)
                sphHarm=0.25*SQRT(5.d0/pi)*(3.d0*cos2-1.d0)
                RETURN
              CASE(1)
                sphHarm=-0.5*SQRT(7.5/pi)*EXP(im*phi)*SIN(theta)*COS(theta)
                RETURN
              CASE(2)
                sphHarm=0.25*SQRT(7.5/pi)*EXP(2.d0*im*phi)*sin2
                RETURN
            END SELECT
          CASE(3)
            SELECT CASE(m)
              CASE(-3)
                sphHarm=0.125*SQRT(35.d0/pi)*EXP(-3.d0*im*phi)*sin3
                RETURN
              CASE(-2)
                sphHarm=0.25*SQRT(52.5/pi)*EXP(-2.d0*im*phi)*sin2*COS(theta)
                RETURN
              CASE(-1)
                sphHarm=0.125*SQRT(21.d0/pi)*EXP(-1.d0*im*phi)*SIN(theta)*(5.d0*cos2-1.d0)
                RETURN
              CASE(0)
                sphHarm=0.25*SQRT(7.d0/pi)*(5.d0*cos3-3.d0*COS(theta))
                RETURN
              CASE(1)
                sphHarm=-0.125*SQRT(21.d0/pi)*EXP(im*phi)*SIN(theta)*(5.d0*cos2-1.d0)
                RETURN
              CASE(2)
                sphHarm=0.25*SQRT(52.5/pi)*EXP(2.d0*im*phi)*sin2*COS(theta)
                RETURN
              CASE(3)
                sphHarm=-0.125*SQRT(35.d0/pi)*EXP(3.d0*im*phi)*sin3
                RETURN
            END SELECT
          CASE(4)
            SELECT CASE(m)
              CASE(-4)
                sphHarm=0.1875*SQRT(17.5/pi)*EXP(-4.d0*im*phi)*sin4
                RETURN
              CASE(-3)
                sphHarm=0.375*SQRT(35.d0/pi)*EXP(-3.d0*im*phi)*sin3*COS(theta)
                RETURN
              CASE(-2)
                sphHarm=0.375*SQRT(2.5/pi)*EXP(-2.d0*im*phi)*sin2*(7.d0*cos2-1.d0)
                RETURN
              CASE(-1)
                sphHarm=0.375*SQRT(5.d0/pi)*EXP(-1.d0*im*phi)*SIN(theta)*(7.d0*cos3-3.d0*COS(theta))
                RETURN
              CASE(0)
                sphHarm=0.1875*SQRT(1.d0/pi)*(35.d0*cos4-30.d0*cos2+3.d0)
                RETURN
              CASE(1)
                sphHarm=-0.375*SQRT(5.d0/pi)*EXP(im*phi)*SIN(theta)*(7.d0*cos3-3.d0*COS(theta))
                RETURN
              CASE(2)
                sphHarm=0.375*SQRT(2.5/pi)*EXP(2.d0*im*phi)*sin2*(7.d0*cos2-1.d0)
                RETURN
              CASE(3)
                sphHarm=-0.375*SQRT(35.d0/pi)*EXP(3.d0*im*phi)*sin3*COS(theta)
                RETURN
              CASE(4)
                sphHarm=0.1875*SQRT(17.5/pi)*EXP(4.d0*im*phi)*sin4
                RETURN
            END SELECT
          CASE(5)
            SELECT CASE(m)
              CASE(-5)
                sphHarm=0.09375*SQRT(77.d0/pi)*EXP(-5.d0*im*phi)*sin5
                RETURN
              CASE(-4)
                sphHarm=0.1875*SQRT(192.5/pi)*EXP(-4.d0*im*phi)*sin4*COS(theta)
                RETURN
              CASE(-3)
                sphHarm=0.03125*SQRT(385.d0/pi)*EXP(-3.d0*im*phi)*sin3*(9.d0*cos2-1.d0)
                RETURN
              CASE(-2)
                sphHarm=0.125*SQRT(577.5/pi)*EXP(-2.d0*im*phi)*sin2*(3.d0*cos3-COS(theta))
                RETURN
              CASE(-1)
                sphHarm=0.0625*SQRT(82.5/pi)*EXP(-1.d0*im*phi)*SIN(theta)*(21.d0*cos4-14.d0*cos2+1.d0)
                RETURN
              CASE(0)
                sphHarm=0.0625*SQRT(11.d0/pi)*(63.d0*cos5-70.d0*cos3+15.d0*COS(theta))
                RETURN
              CASE(1)
                sphHarm=-0.0625*SQRT(82.5/pi)*EXP(im*phi)*SIN(theta)*(21.d0*cos4-14.d0*cos2+1.d0)
                RETURN
              CASE(2)
                sphHarm=0.125*SQRT(577.5/pi)*EXP(2.d0*im*phi)*sin2*(3.d0*cos3-COS(theta))
                RETURN
              CASE(3)
                sphHarm=-0.03125*SQRT(385.d0/pi)*EXP(3.d0*im*phi)*sin3*(9.d0*cos2-1.d0)
                RETURN
              CASE(4)
                sphHarm=0.1875*SQRT(192.5/pi)*EXP(4.d0*im*phi)*sin4*COS(theta)
                RETURN
              CASE(5)
                sphHarm=-0.09375*SQRT(77.d0/pi)*EXP(5.d0*im*phi)*sin5
                RETURN
            END SELECT
          CASE(6)
            SELECT CASE(m)
              CASE(-6)
                sphHarm=0.015625*SQRT(3303.d0/pi)*EXP(-6.d0*im*phi)*sin6
                RETURN
              CASE(-5)
                sphHarm=0.09375*SQRT(1001.d0/pi)*EXP(-5.d0*im*phi)*sin5*COS(theta)
                RETURN
              CASE(-4)
                sphHarm=0.09375*SQRT(45.5/pi)*EXP(-4.d0*im*phi)*sin4*(11.d0*cos2-1.d0)
                RETURN
              CASE(-3)
                sphHarm=0.03125*SQRT(1365.d0/pi)*EXP(-3.d0*im*phi)*sin3*(11.d0*cos3-3.d0*COS(theta))
                RETURN
              CASE(-2)
                sphHarm=0.015625*SQRT(1365.d0/pi)*EXP(-2.d0*im*phi)*sin2*(33.d0*cos4-18.d0*cos2+1.d0)
                RETURN
              CASE(-1)
                sphHarm=0.0625*SQRT(136.5/pi)*EXP(-1.d0*im*phi)*SIN(theta)*(33.d0*cos5-30.d0*cos3+5.d0*COS(theta))
                RETURN
              CASE(0)
                sphHarm=0.03125*SQRT(13.d0/pi)*(231.d0*cos6-315.d0*cos4+105.d0*cos2-5.d0)
                RETURN
              CASE(1)
                sphHarm=-0.0625*SQRT(136.5/pi)*EXP(im*phi)*SIN(theta)*(33.d0*cos5-30.d0*cos3+5.d0*COS(theta))
                RETURN
              CASE(2)
                sphHarm=0.015625*SQRT(1365.d0/pi)*EXP(2.d0*im*phi)*sin2*(33.d0*cos4-18.d0*cos2+1.d0)
                RETURN
              CASE(3)
                sphHarm=-0.03125*SQRT(1365.d0/pi)*EXP(3.d0*im*phi)*sin3*(11.d0*cos3-3.d0*COS(theta))
                RETURN
              CASE(4)
                sphHarm=0.09375*SQRT(45.5/pi)*EXP(4.d0*im*phi)*sin4*(11.d0*cos2-1.d0)
                RETURN
              CASE(5)
                sphHarm=-0.09375*SQRT(1001.d0/pi)*EXP(5.d0*im*phi)*sin5*COS(theta)
                RETURN
              CASE(6)
                sphHarm=0.015625*SQRT(3003.d0/pi)*EXP(6.d0*im*phi)*sin6
                RETURN
            END SELECT
          CASE(7)
            SELECT CASE(m)
              CASE(-7)
                sphHarm=0.046875*SQRT(357.5/pi)*EXP(-7.d0*im*phi)*sin7
                RETURN
              CASE(-6)
                sphHarm=0.046875*SQRT(5005.d0/pi)*EXP(-6.d0*im*phi)*sin6*COS(theta)
                RETURN
              CASE(-5)
                sphHarm=0.046875*SQRT(192.5/pi)*EXP(-5.d0*im*phi)*sin5*(13.d0*cos2-1.d0)
                RETURN
              CASE(-4)
                sphHarm=0.09375*SQRT(192.5/pi)*EXP(-4.d0*im*phi)*sin4*(13.d0*cos3-3.d0*COS(theta))
                RETURN
              CASE(-3)
                sphHarm=0.046875*SQRT(17.5/pi)*EXP(-3.d0*im*phi)*sin3*(143.d0*cos4-66.d0*cos2+3.d0)
                RETURN
              CASE(-2)
                sphHarm=0.046875*SQRT(35.d0/pi)*EXP(-2.d0*im*phi)*sin2*(143.d0*cos5-110.d0*cos3+15.d0*COS(theta))
                RETURN
              CASE(-1)
                sphHarm=0.015625*SQRT(52.5/pi)*EXP(-1.d0*im*phi)*SIN(theta)*(429.d0*cos6-495.d0*cos4+135.d0*cos2-5.d0)
                RETURN
              CASE(0)
                sphHarm=0.03125*SQRT(15.d0/pi)*(429.d0*cos7-693.d0*cos5+315.d0*cos3-35.d0*COS(theta))
                RETURN
              CASE(1)
                sphHarm=-0.015625*SQRT(52.5/pi)*EXP(im*phi)*SIN(theta)*(429.d0*cos6-495.d0*cos4+135.d0*cos2-5.d0)
                RETURN
              CASE(2)
                sphHarm=0.046875*SQRT(35.d0/pi)*EXP(2.d0*im*phi)*sin2*(143.d0*cos5-110.d0*cos3+15.d0*COS(theta))
                RETURN
              CASE(3)
                sphHarm=-0.046875*SQRT(17.5/pi)*EXP(3.d0*im*phi)*sin3*(143.d0*cos4-66.d0*cos2+3.d0)
                RETURN
              CASE(4)
                sphHarm=0.09375*SQRT(192.5/pi)*EXP(4.d0*im*phi)*sin4*(13.d0*cos3-3.d0*COS(theta))
                RETURN
              CASE(5)
                sphHarm=-0.046875*SQRT(192.5/pi)*EXP(5.d0*im*phi)*sin5*(13.d0*cos2-1.d0)
                RETURN
              CASE(6)
                sphHarm=0.046875*SQRT(5005.d0/pi)*EXP(6.d0*im*phi)*sin6*COS(theta)
                RETURN
              CASE(7)
                sphHarm=-0.046875*SQRT(357.5/pi)*EXP(7.d0*im*phi)*sin7
                RETURN
            END SELECT
          CASE(8)
            SELECT CASE(m)
              CASE(-8)
                sphHarm=0.01171875*SQRT(6077.5/pi)*EXP(-8.d0*im*phi)*sin8
                RETURN
              CASE(-7)
                sphHarm=0.046875*SQRT(6077.5/pi)*EXP(-7.d0*im*phi)*sin7*COS(theta)
                RETURN
              CASE(-6)
                sphHarm=0.0078125*SQRT(7293.d0/pi)*EXP(-6.d0*im*phi)*sin6*(15.d0*cos2-1.d0)
                RETURN
              CASE(-5)
                sphHarm=0.046875*SQRT(8508.5/pi)*EXP(-5.d0*im*phi)*sin5*(5.d0*cos3-COS(theta))
                RETURN
              CASE(-4)
                sphHarm=0.0234375*SQRT(654.5/pi)*EXP(-4.d0*im*phi)*sin4*(65.d0*cos4-26.d0*cos2+1.d0)
                RETURN
              CASE(-3)
                sphHarm=0.015625*SQRT(9817.5/pi)*EXP(-3.d0*im*phi)*sin3*(39.d0*cos5-26.d0*cos3+3.d0*COS(theta))
                RETURN
              CASE(-2)
                sphHarm=0.0234375*SQRT(595.d0/pi)*EXP(-2.d0*im*phi)*sin2*(143.d0*cos6-143.d0*cos4+33.d0*cos2-1.d0)
                RETURN
              CASE(-1)
                sphHarm=0.046875*SQRT(8.5/pi)*EXP(-1.d0*im*phi)*SIN(theta)*(715.d0*cos7-1001.d0*cos5+385.d0*cos3-35.d0*COS(theta))
                RETURN
              CASE(0)
                sphHarm=0.00390625*SQRT(17.d0/pi)*(6435.d0*cos8-12012.d0*cos6+6930.d0*cos4-1260.d0*cos2+35.d0)
                RETURN
              CASE(1)
                sphHarm=-0.046875*SQRT(8.5/pi)*EXP(im*phi)*SIN(theta)*(715.d0*cos7-1001.d0*cos5+385.d0*cos3-35.d0*COS(theta))
                RETURN
              CASE(2)
                sphHarm=0.0234375*SQRT(595.d0/pi)*EXP(2.d0*im*phi)*sin2*(143.d0*cos6-143.d0*cos4+33.d0*cos2-1.d0)
                RETURN
              CASE(3)
                sphHarm=-0.015625*SQRT(9817.5/pi)*EXP(3.d0*im*phi)*sin3*(39.d0*cos5-26.d0*cos3+3.d0*COS(theta))
                RETURN
              CASE(4)
                sphHarm=0.0234375*SQRT(654.5/pi)*EXP(4.d0*im*phi)*sin4*(65.d0*cos4-26.d0*cos2+1.d0)
                RETURN
              CASE(5)
                sphHarm=-0.046875*SQRT(8508.5/pi)*EXP(5.d0*im*phi)*sin5*(5.d0*cos3-COS(theta))
                RETURN
              CASE(6)
                sphHarm=0.0078125*SQRT(7293.d0/pi)*EXP(6.d0*im*phi)*sin6*(15.d0*cos2-1.d0)
                RETURN
              CASE(7)
                sphHarm=-0.046875*SQRT(6077.5/pi)*EXP(7.d0*im*phi)*sin7*COS(theta)
                RETURN
              CASE(8)
                sphHarm=0.01171875*SQRT(6077.5/pi)*EXP(8.d0*im*phi)*sin8
                RETURN
            END SELECT
          CASE(9)
            SELECT CASE(m)
              CASE(-9)
                sphHarm=0.001953125*SQRT(230945.d0/pi)*EXP(-9.d0*im*phi)*sin9
                RETURN
              CASE(-8)
                sphHarm=0.01171875*SQRT(115472.5/pi)*EXP(-8.d0*im*phi)*sin8*COS(theta)
                RETURN
              CASE(-7)
                sphHarm=0.005859375*SQRT(13585.d0/pi)*EXP(-7.d0*im*phi)*sin7*(17.d0*cos2-1.d0)
                RETURN
              CASE(-6)
                sphHarm=0.0078125*SQRT(40755.d0/pi)*EXP(-6.d0*im*phi)*sin6*(17.d0*cos3-3.d0*COS(theta))
                RETURN
              CASE(-5)
                sphHarm=0.01171875*SQRT(2717.d0/pi)*EXP(-5.d0*im*phi)*sin5*(85.d0*cos4-30.d0*cos2+1.d0)
                RETURN
              CASE(-4)
                sphHarm=0.0234375*SQRT(47547.5/pi)*EXP(-4.d0*im*phi)*sin4*(17.d0*cos5-10.d0*cos3+COS(theta))
                RETURN
              CASE(-3)
                sphHarm=0.00390625*SQRT(21945.d0/pi)*EXP(-3.d0*im*phi)*sin3*(221.d0*cos6-195.d0*cos4+39.d0*cos2-1.d0)
                RETURN
              CASE(-2)
                sphHarm=0.0234375*SQRT(1045.d0/pi)*EXP(-2.d0*im*phi)*sin2*(221.d0*cos7-273.d0*cos5+91.d0*cos3-7.d0*COS(theta))
                RETURN
              CASE(-1)
                sphHarm=0.01171875*SQRT(47.5/pi)*EXP(-1.d0*im*phi)*SIN(theta)* &
                      & (2431.d0*cos8-4004.d0*cos6+2002.d0*cos4-308.d0*cos2+7.d0)
                RETURN
              CASE(0)
                sphHarm=0.00390625*SQRT(19.d0/pi)*(12155.d0*cos9-25740.d0*cos7+18018.d0*cos5-4620.d0*cos3+315.d0*COS(theta))
                RETURN
              CASE(1)
                sphHarm=-0.01171875*SQRT(47.5/pi)*EXP(im*phi)*SIN(theta)*(2431.d0*cos8-4004.d0*cos6+2002.d0*cos4-308.d0*cos2+7.d0)
                RETURN
              CASE(2)
                sphHarm=0.0234375*SQRT(1045.d0/pi)*EXP(2.d0*im*phi)*sin2*(221.d0*cos7-273.d0*cos5+91.d0*cos3-7.d0*COS(theta))
                RETURN
              CASE(3)
                sphHarm=-0.00390625*SQRT(21945.d0/pi)*EXP(3.d0*im*phi)*sin3*(221.d0*cos6-195.d0*cos4+39.d0*cos2-1.d0)
                RETURN
              CASE(4)
                sphHarm=0.0234375*SQRT(47547.5/pi)*EXP(4.d0*im*phi)*sin4*(17.d0*cos5-10.d0*cos3+COS(theta))
                RETURN
              CASE(5)
                sphHarm=-0.01171875*SQRT(2717.d0/pi)*EXP(5.d0*im*phi)*sin5*(85.d0*cos4-30.d0*cos2+1.d0)
                RETURN
              CASE(6)
                sphHarm=0.0078125*SQRT(40755.d0/pi)*EXP(6.d0*im*phi)*sin6*(17.d0*cos3-3.d0*COS(theta))
                RETURN
              CASE(7)
                sphHarm=-0.005859375*SQRT(13585.d0/pi)*EXP(7.d0*im*phi)*sin7*(17.d0*cos2-1.d0)
                RETURN
              CASE(8)
                sphHarm=0.01171875*SQRT(115472.5/pi)*EXP(8.d0*im*phi)*sin8*COS(theta)
                RETURN
              CASE(9)
                sphHarm=-0.001953125*SQRT(230945.d0/pi)*EXP(-9.d0*im*phi)*sin9
                RETURN
            END SELECT
          CASE(10)
            SELECT CASE(m)
              CASE(-10)
                sphHarm=0.0009765625*SQRT(969969.d0/pi)*EXP(-10.d0*im*phi)*sin10
                RETURN
              CASE(-9)
                sphHarm=0.001953125*SQRT(4849845.d0/pi)*EXP(-9.d0*im*phi)*sin9*COS(theta)
                RETURN
              CASE(-8)
                sphHarm=0.001953125*SQRT(127627.5/pi)*EXP(-8.d0*im*phi)*sin8*(19.d0*cos2-1.d0)
                RETURN
              CASE(-7)
                sphHarm=0.005859375*SQRT(85085.d0/pi)*EXP(-7.d0*im*phi)*sin7*(19.d0*cos3-3.d0*COS(theta))
                RETURN
              CASE(-6)
                sphHarm=0.0029296875*SQRT(5005.d0/pi)*EXP(-6.d0*im*phi)*sin6*(323.d0*cos4-102.d0*cos2+3.d0)
                RETURN
              CASE(-5)
                sphHarm=0.01171875*SQRT(1001.d0/pi)*EXP(-5.d0*im*phi)*sin5*(323.d0*cos5-170.d0*cos3+15.d0*COS(theta))
                RETURN
              CASE(-4)
                sphHarm=0.01171875*SQRT(2502.5/pi)*EXP(-4.d0*im*phi)*sin4*(323.d0*cos6-255.d0*cos4+45.d0*cos2-1.d0)
                RETURN
              CASE(-3)
                sphHarm=0.01171875*SQRT(5005.d0/pi)*EXP(-3.d0*im*phi)*sin3*(323.d0*cos7-357.d0*cos5+105.d0*cos3-7.d0*COS(theta))
                RETURN
              CASE(-2)
                sphHarm=0.005859375*SQRT(192.5/pi)*EXP(-2.d0*im*phi)*sin2*(4199.d0*cos8-6188.d0*cos6+2730.d0*cos4-364.d0*cos2+7.d0)
                RETURN
              CASE(-1)
                sphHarm=0.00390625*SQRT(577.5/pi)*EXP(-1.d0*im*phi)*SIN(theta)* &
                      & (4199.d0*cos9-7956.d0*cos7+4914.d0*cos5-1092.d0*cos3+63.d0*COS(theta))
                RETURN
              CASE(0)
                sphHarm=0.001953125*SQRT(21.d0/pi)*(46189.d0*cos10-109395.d0*cos8+90090.d0*cos6-30030.d0*cos4+3465.d0*cos2-63.d0)
                RETURN
              CASE(1)
                sphHarm=-0.00390625*SQRT(577.5/pi)*EXP(im*phi)*SIN(theta)* &
                      & (4199.d0*cos9-7956.d0*cos7+4914.d0*cos5-1092.d0*cos3+63.d0*COS(theta))
                RETURN
              CASE(2)
                sphHarm=0.005859375*SQRT(192.5/pi)*EXP(2.d0*im*phi)*sin2*(4199.d0*cos8-6188.d0*cos6+2730.d0*cos4-364.d0*cos2+7.d0)
                RETURN
              CASE(3)
                sphHarm=-0.01171875*SQRT(5005.d0/pi)*EXP(3.d0*im*phi)*sin3*(323.d0*cos7-357.d0*cos5+105.d0*cos3-7.d0*COS(theta))
                RETURN
              CASE(4)
                sphHarm=0.01171875*SQRT(2502.5/pi)*EXP(4.d0*im*phi)*sin4*(323.d0*cos6-255.d0*cos4+45.d0*cos2-1.d0)
                RETURN
              CASE(5)
                sphHarm=-0.01171875*SQRT(1001.d0/pi)*EXP(5.d0*im*phi)*sin5*(323.d0*cos5-170.d0*cos3+15.d0*COS(theta))
                RETURN
              CASE(6)
                sphHarm=0.0029296875*SQRT(5005.d0/pi)*EXP(6.d0*im*phi)*sin6*(323.d0*cos4-102.d0*cos2+3.d0)
                RETURN
              CASE(7)
                sphHarm=-0.005859375*SQRT(85085.d0/pi)*EXP(7.d0*im*phi)*sin7*(19.d0*cos3-3.d0*COS(theta))
                RETURN
              CASE(8)
                sphHarm=0.001953125*SQRT(127627.5/pi)*EXP(8.d0*im*phi)*sin8*(19.d0*cos2-1.d0)
                RETURN
              CASE(9)
                sphHarm=-0.001953125*SQRT(4849845.d0/pi)*EXP(9.d0*im*phi)*sin9*COS(theta)
                RETURN
              CASE(10)
                sphHarm=0.0009765625*SQRT(969969.d0/pi)*EXP(10.d0*im*phi)*sin10
                RETURN
            END SELECT
          CASE(11:)
            WRITE(*,*) "l>10 for spherical harmonics not supported yet"
            RETURN
        END SELECT
!
        END FUNCTION


elemental function __APPEND(erfinv,__PREC) (y) result(res)
    !*  ERFINV returns the value X for a given value Y such that
    !   Y = ERF(x), ie ERFINV implements the inverse function of the
    !   Error function ERF.
    !
    !   Note on implementation: ERFINV calls the inverse CDF of a
    !   standard-normal random variable implemented in NDTRI ported from
    !   the Cephes math library.
    !   We use the relationship
    !       NDTRI(z) = sqrt(2) * ERFINV(2*y - 1)
    !   to compute ERFINV as
    !       ERFINV(y) = NDTRI((y+1)/2) / sqrt(2)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in) :: y
    real (PREC) :: res

    res = ndtri((y + 1.0_PREC)/2.0_PREC) / sqrt(2.0_PREC)
end function


elemental function __APPEND(ndtri,__PREC) (y) result(res)
    !*  NDTRI returns the argument x for which the area under the Gaussian
    !   probability density function (integrated from -inf to x) is equal to
    !   y.
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in) :: y
    real (PREC) :: res

    real (PREC), parameter :: S2PI = sqrt(2.0_PREC * __APPEND(PI,__PREC))
        !   sqrt(2.0 * pi)
    real (PREC), parameter :: EXPM2 = exp(-2.0_PREC)

    real (PREC), parameter :: P0(*) = [ real (PREC) :: &
        -5.99633501014107895267e1_PREC, &
         9.80010754185999661536e1_PREC, &
        -5.66762857469070293439e1_PREC, &
         1.39312609387279679503e1_PREC, &
        -1.23916583867381258016e0_PREC &
    ]

    real (PREC), parameter :: Q0(*) = [ real (PREC) :: &
         1.95448858338141759834e0_PREC, &
         4.67627912898881538453e0_PREC, &
         8.63602421390890590575e1_PREC, &
        -2.25462687854119370527e2_PREC, &
         2.00260212380060660359e2_PREC, &
        -8.20372256168333339912e1_PREC, &
         1.59056225126211695515e1_PREC, &
        -1.18331621121330003142e0_PREC &
    ]

    real (PREC), parameter :: P1(*) = [ real (PREC) :: &
         4.05544892305962419923e0_PREC, &
         3.15251094599893866154e1_PREC, &
         5.71628192246421288162e1_PREC, &
         4.40805073893200834700e1_PREC, &
         1.46849561928858024014e1_PREC, &
         2.18663306850790267539e0_PREC, &
        -1.40256079171354495875e-1_PREC, &
        -3.50424626827848203418e-2_PREC, &
        -8.57456785154685413611e-4_PREC &
    ]

    real (PREC), parameter :: Q1(*) = [ real (PREC) :: &
         1.57799883256466749731e1_PREC, &
         4.53907635128879210584e1_PREC, &
         4.13172038254672030440e1_PREC, &
         1.50425385692907503408e1_PREC, &
         2.50464946208309415979e0_PREC, &
        -1.42182922854787788574e-1_PREC, &
        -3.80806407691578277194e-2_PREC, &
        -9.33259480895457427372e-4_PREC &
    ]

    real (PREC), parameter :: P2(*) = [ real (PREC) :: &
         3.23774891776946035970e0_PREC, &
         6.91522889068984211695e0_PREC, &
         3.93881025292474443415e0_PREC, &
         1.33303460815807542389e0_PREC, &
         2.01485389549179081538e-1_PREC, &
         1.23716634817820021358e-2_PREC, &
         3.01581553508235416007e-4_PREC, &
         2.65806974686737550832e-6_PREC, &
         6.23974539184983293730e-9_PREC &
    ]

    real (PREC), parameter :: Q2(*) = [ real (PREC) :: &
         6.02427039364742014255e0_PREC, &
         3.67983563856160859403e0_PREC, &
         1.37702099489081330271e0_PREC, &
         2.16236993594496635890e-1_PREC, &
         1.34204006088543189037e-2_PREC, &
         3.28014464682127739104e-4_PREC, &
         2.89247864745380683936e-6_PREC, &
         6.79019408009981274425e-9_PREC &
    ]

    real (PREC) :: x, y1, z, y2, x0, x1
    integer :: code

    if (y <= 0.0_PREC) then
        x = ieee_value (0.0_PREC, IEEE_NEGATIVE_INF)
        goto 100
    else if (y >= 1.0_PREC) then
        x = ieee_value (0.0_PREC, IEEE_POSITIVE_INF)
        goto 100
    end if

    code = 1
    y1 = y
    if (y1 > (1.0_PREC - EXPM2)) then
        y1 = 1.0_PREC - y1
        code = 0
    end if

    if (y1 > EXPM2) then
        y1 = y1 - 0.5_PREC
        y2 = y1 * y1
        x = y1 + y1 * (y2 * polevl(y2, P0(1:5)) / p1evl(y2, Q0(1:8)))
        x = x * S2PI
        goto 100
    end if

    x = sqrt(-2.0_PREC * log(y1))
    x0 = x - log(x) / x

    z = 1.0_PREC / x
    if (x < 8.0_PREC) then
        !* y1 > exp(-32) = 1.26e-14
        x1 = z * polevl(z, P1(1:9)) / p1evl(z, Q1(1:8))
    else
        x1 = z * polevl(z, P2(1:9)) / p1evl(z, Q2(1:8))
    end if

    x = x0 - x1
    if (code /= 0) then
        x = -x
    end if


100 continue
    res = x

end function



pure function __APPEND(polevl,__PREC) (x, coef) result(res)
    !*  POLEVL is a quick-and-dirty routine to evaluate a polynomial of
    !   degree n = size(coef) - 1 given by
    !       p(x) = a_0 + a_1 x + ... + a_n x^n
    !   with coefficient vector COEF at the point X.
    !
    !   This routine is a Fortran implementation of polevl() from the Cephes
    !   math library (as used in Scipy). This routine should not
    !   be used within this module.
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in) :: x
        !*  Point at which polynomial should be evaluated.
    real (PREC), intent(in), dimension(:) :: coef
        !*  Array of polynomial coefficients stored in REVERSE order, ie
        !   coef[1] = a_n
    real (PREC) :: res

    integer :: i

    res = coef(1)

    do i = 2, size(coef)
        res = res * x + coef(i)
    end do

end function



pure function __APPEND(p1evl,__PREC) (x, coef) result(res)
    !*  P1EVL is a quick-and-dirty routine to evaluate a polynomial of degree
    !   n = size(coef) given by
    !       p(x) = a_0 + a_1 x + ... + a_n x^n
    !   This routine assumes that a_n = 1 and is omitted from the COEF array,
    !   otherwise it is the same as POLEVL above.
    !
    !   This routine is a Fortran implementation of p1evl() from the Cephes
    !   math library (as used in Scipy). This routine should not
    !   be used within this module.
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in) :: x
    real (PREC), intent(in), dimension(:) :: coef
        !*  Array of polynomial coefficients stored in REVERSE order, ie
        !   coef[1] = a_{n-1}
    real (PREC) :: res

    integer :: i

    res = x + coef(1)

    do i = 2, size(coef)
        res = res * x + coef(i)
    end do

end function

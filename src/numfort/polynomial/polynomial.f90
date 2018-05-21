module numfort_polynomial

    use numfort_polynomial_create, only: polyder, polyshift, polyint

    use numfort_polynomial_polyfit
    use numfort_polynomial_polyroots
    use numfort_polynomial_polyval, only: polyval

    use numfort_polynomial_complete
    use numfort_polynomial_polyval_complete

    use numfort_polynomial_ppoly


end module

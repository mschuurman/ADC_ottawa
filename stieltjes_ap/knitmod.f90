  module knitmod

    use constants

    save

    integer                            :: npntb,npntc
    real(d)                            :: et
    real(d), dimension(:), allocatable :: fb,fc,eb,ec
    character(len=120)                 :: abound,acont

  end module knitmod

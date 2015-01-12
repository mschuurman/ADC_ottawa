  program h2p
    use multigrid

    call MultiGridInit(5,6,1)
    call SimpleGridNew('The Universe',(/  20,          15,          20 /), &
                             reshape( (/ -10d0, 10d0, -20d0, 20d0, -30d0, 30d0 /), (/ 2, 3 /)))
    call SimpleGridNew('Molecule 1'  ,(/   5,           5,           5 /), &
                             reshape( (/ -1d0,   1d0,  -2d0,   2d0, -3d0,   3d0 /), (/ 2, 3 /)))
    call SimpleGridNew('Molecule 2'  ,(/   5,           5,           5 /), &
                             reshape( (/  5d0,   6d0,   7d0,   8d0,  9d0,  10d0 /), (/ 2, 3 /)))
  end program h2p

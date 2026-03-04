```Fortran
!Pressure gradient
    work = aa * h0bf + (1.0 - aa) * h0bl  
    wkk(:,:,1) = (psa * OD0 + work * grav) * vit(:,:,1)
    do k = 1, npz
      wkk(:,:,k+1) = wkk(:,:,k) - gg(:,:,k) * dzp(k) * vit(:,:,k)
      wka(:,:,k  ) = 0.5D0 * (wkk(:,:,k) + wkk(:,:,k+1))
      pres(:,:,k) = wka(:,:,k)
      call op%grad%scalar(gradx, grady, wka(:,:,k), dg, k)
      pgf_x(:,:,k) = - gradx
      pgf_y(:,:,k) = - grady
    enddo
  !G'DH/DX
    do k = 1, npz
      wka(:,:,k) = (1.0 + ohbt * zkt(k)) * work * vit(:,:,k)
      call op%grad%scalar(gradx, grady, wka(:,:,k), dg, k)
      pres(:,:,k) = pres(:,:,k) + wka(:,:,k) * gg(:,:,k)
      pgf_x(:,:,k) = pgf_x(:,:,k) + gg(:,:,k) * gradx
      pgf_y(:,:,k) = pgf_y(:,:,k) + gg(:,:,k) * grady
    enddo
```
这部分对应的是

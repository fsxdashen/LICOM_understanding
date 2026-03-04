```Fortran
subroutine pgf(pgf_x,pgf_y,momentum,psa,aa)
```
该代码计算了两部分的压力梯度力，所参考的公式如下 (licom手册P31)

$$\frac{du}{dt}=-\frac{1}{\rho _0}\frac{\partial p}{\partial x}+g\prime \frac{\partial}{\partial x}\left[ \left( 1+\frac{H_m}{H_b}\eta \right) z_0 \right] +f*v+F_x$$

$$\frac{dv}{dt}=-\frac{1}{\rho _0}\frac{\partial p}{\partial y}+g\prime \frac{\partial}{\partial y}\left[ \left( 1+\frac{H_m}{H_b}\eta \right) z_0 \right] +f*u+F_y$$

其中 $H_m\eta\sim z_{k+0.5}(zkt)$

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
!G DH/DX
    do k = 1, npz
      wka(:,:,k) = (1.0 + ohbt * zkt(k)) * work * vit(:,:,k)
      call op%grad%scalar(gradx, grady, wka(:,:,k), dg, k)
      pres(:,:,k) = pres(:,:,k) + wka(:,:,k) * gg(:,:,k)
      pgf_x(:,:,k) = pgf_x(:,:,k) + gg(:,:,k) * gradx
      pgf_y(:,:,k) = pgf_y(:,:,k) + gg(:,:,k) * grady
    enddo
```
其中各变量的对应情况如下

$$wkk(:,:,1)=\frac{p_{as}}{\rho _0}+gz_0$$

之后每一层的wkk则通过以下公式进行计算

$$\frac{\partial p}{\partial z}=-\rho g\prime $$

其中 $gg\sim g\prime=-g\rho/\rho_0$

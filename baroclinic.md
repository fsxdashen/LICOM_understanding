```Fortran
subroutine pgf(pgf_x,pgf_y,momentum,psa,aa)
```
该代码计算了两部分的压力梯度力，所参考的公式如下 (licom手册P31)

$$\frac{\partial u}{\partial t}=-\vec{v}\cdot\nabla u-\frac{1}{\rho _0}\frac{\partial p}{\partial x}+g\prime \frac{\partial}{\partial x}\left[ \left( 1+\frac{H_m}{H_b}\eta \right) z_0 \right] + f*v+F_x$$

$$\frac{\partial v}{\partial t}=-\vec{v}\cdot\nabla v-\frac{1}{\rho _0}\frac{\partial p}{\partial y}+g\prime \frac{\partial}{\partial y}\left[ \left( 1+\frac{H_m}{H_b}\eta \right) z_0 \right] - f*u+F_y$$

其中 $H_m\eta\sim z_{k+0.5}(zkt)$, $gg\sim g\prime=-g\rho/\rho_0$

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

$$\frac{\partial p}{\partial z}=-\rho g\prime \Longrightarrow \frac{\partial}{\partial z}(\frac{p}{\rho})=-g\prime$$

```op%grad%scalar```算子可用于计算wka的梯度，其中x方向上的梯度存放在```gradx```, y方向上的梯度存放在```grady```

---

```Fortran
subroutine updatedluv(momentum,pgf_x,pgf_y)
```
在上一步已经在右端项上加入压力梯度的基础上，进一步加上科氏力
```Fortran
dlv(:,:,k) = dlv(:,:,k) + pgf_y(:,:,k) - a_f * uu_ct(:,:,k)
dlu(:,:,k) = dlu(:,:,k) + pgf_x(:,:,k) + a_f * vv_ct(:,:,k)
```

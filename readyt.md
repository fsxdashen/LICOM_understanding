```Fortran
DO K = 1, npz-1
  DO J = jsd, jed
    DO I = isd, ied
      TUP = ATB (I,J,K,1) - TO (K +1)
      SUP = ATB (I,J,K,2) - SO (K +1)
      TLO = ATB (I,J,K +1,1) - TO (k +1)
      SLO = ATB (I,J,K +1,2) - SO (K +1)
      RHOUP = op%others%dens(TUP, SUP, C(k+1,:))
      RHOLO = op%others%dens(TLO, SLO, C(k+1,:))
      rict(I,J,K) = VIT (I,J,K +1)* OD0* G * (RHOLO - RHOUP)*ODZT(K+1)
    END DO
  END DO
END DO
```
最终需要得到的条件如下

$$R_i=\frac{g}{\rho _0}\frac{\partial \rho _{pot}}{\partial z}/\left[ \left( \frac{\partial u}{\partial z} \right) ^2+\left( \frac{\partial v}{\partial z} \right) ^2 \right] $$

其中 $\rho_{pot}$为位密度，同样用 Bryan 和 Cox（1972）的九项式计算，垂直方向相邻两层的位密度计算参照两层中较深一层。

该代码先考虑除剪切平方项之前的部分，并得到中间变量**rict**，具体为

$$rict=\frac{g}{\rho _0}\frac{\delta \rho \left( T_{up},S_{up} \right) -\delta \rho \left( T_{low},S_{low} \right)}{\delta z}=\frac{g}{\rho _0}\frac{\delta \rho _{up}-\delta \rho _{low}}{\delta z}=\frac{g}{\rho _0}\frac{\partial \rho _{pot}}{\partial z}$$

后面的剪切平方项的计算出现在readyc.F90中，所使用的代码为
```Fortran
call get_s2trit(s2t,rit,tracer%rict,dg,u,v,epsln)
```
在这里实现最终对 $Ri$的计算

---

```Fortran
dlu(:,:,1) = 0.0D0
dlu(:,:,2) = 0.0D0
 
DO K = 1, npz
  DO J = jsd, jed
    DO I = isd, ied
      ABCD = GG (I,J,K)* dg%OHBT (I,J)* dg%DZP (K)
      DLU (I,J,1)= DLU (I,J,1) + ABCD
      DLU (I,J,2)= DLU (I,J,2) + ABCD * dg%ZKT (K)
    END DO
  END DO
END DO
 
 
DO J = jsd, jed
  DO I = isd, ied
    DLV (I,J,1)= (DLU (I,J,1) + DLU(I,J,2)*dg%OHBT (I,J))/ G
    DLV (I,J,2)= DLU (I,J,2)* dg%OHBT(I,J)*dg%OHBT(I,J)
  END DO
END DO
!
wgp(:,:) = dlv(:,:,1)
whx(:,:) = hbx(:,:)*dlv(:,:,2)*vit(:,:,1)
why(:,:) = hby(:,:)*dlv(:,:,2)*vit(:,:,1)
```
中间变量如下
| 名称 | 公式 | 说明 | 单位 |
|---|---|---|---|
| HBX(IMT,JMT) | $\frac{1}{a\sin\theta_{j+1/2}\Delta\lambda}\left(\overline{(H_b)^{\lambda}}\right)_{i-1/2,j+1/2}$ | 海底地形在 U 网格上的纬向梯度 |  |
| HBY(IMT,JMT) | $\frac{1}{a\Delta\theta_{j+1/2}}\left(\overline{(H_b)^{\theta}}\right)_{i-1/2,j+1/2}$ | 海底地形在 U 网格上的经向梯度 |  |
| OHBT(IMT,JMT) | $(H_b)^{-1}_{i,j}$ | T 网格水柱厚度的倒数 | $m^{-1}$ |
| DZP(KM) | $\Delta z_{k+1/2}=z_k-z_{k+1}$ | T 网格元的厚度 | $m$ |
| ZKP(KMP1) | $z_k$ | W 网格的深度，从 INDEX.DATA 读入。注意，此变量在 grids.F90 中定义 | $m$ |

最终需要得到的条件如下：

$$WGP(i,j) = \frac{\overline{g^{\prime}}}{g}+\frac{\overline{g^{\prime}\eta }}{g\eta _s} \qquad WHX(i,j) = \frac{\overline{g^{\prime}\eta }}{\eta _s^2}\frac{\partial \eta_s}{\partial y} \qquad WHY(i,j) = \frac{\overline{g^{\prime}\eta }}{\eta _s^2}\frac{\partial \eta _s}{\partial x}$$

模式垂直坐标采用了 $\eta$坐标系统，它和 $z$坐标的关系可表示为

$$\eta \equiv -\frac{z_0-z}{z_0+H_b}\times \eta _s;\eta _s=\frac{H_b}{H_m}$$

其中 $z_0$为海表起伏, $H_b$为准阶梯状的海底地形, $H_m$为最大地形深度。因此，在 $z=z_0$处, $\eta=0$;在 $z=-H_b$处, $\eta=-\eta_s$

这里需要用到一步小的近似，即

$$\eta \equiv -\frac{z_0-z}{z_0+H_b}\times \eta _s\approx \frac{z}{H_b}\eta _s=\frac{z}{H_b}\frac{H_b}{H_m}=\frac{z}{H_m}$$

$\eta$ 坐标系中深度积分定义为 

$$\overline{\left(  \right)}=\frac{1}{\eta _s}\int_{-\eta _s}^0{\left(  \right) d\eta}$$

下面开始根据代码进行推导

$$ABCD=g\prime \frac{1}{H_b}\Delta z$$

$$dlu1=\frac{1}{H_b}\int_{-H_b}^0{g\prime dz}\qquad\qquad   dlu2=\frac{1}{H_b}\int_{-H_b}^0{g\prime zdz}$$

$$dlv1=\frac{1}{g}\left( \frac{1}{H_b}\int_{-H_b}^0{g\prime dz}+\frac{1}{H_{b}^{2}}\int_{-H_b}^0{g\prime zdz} \right)$$

$$dlv2=\frac{1}{H_{b}^{3}}\int_{-H_b}^0{g\prime zdz}=\frac{1}{H_b}\frac{H_{m}^{2}}{H_{b}^{2}}\int_{-\eta _s}^0{g\prime \eta d\eta}=\frac{1}{H_b}\frac{\overline{g\prime \eta }}{\eta _s}$$

最后推导得到这一步需要的三个变量

$$wgp=\frac{1}{g}\left( \frac{1}{H_b}\int_{-H_b}^0{g\prime dz}+\frac{1}{H_{b}^{2}}\int_{-H_b}^0{g\prime zdz} \right) =\frac{1}{g}\left( \frac{H_m}{H_b}\int_{-\eta _s}^0{g\prime d\eta}+\frac{H_{m}^{2}}{H_{b}^{2}}\int_{-\eta _s}^0{g\prime \eta d\eta} \right) =\frac{1}{g}\left( \frac{1}{\eta _s}\int_{-\eta _s}^0{g\prime d\eta}+\frac{1}{\eta _s}\frac{1}{\eta _s}\int_{-\eta _s}^0{g\prime \eta d\eta} \right) =\frac{1}{g}\left( \overline{g\prime }+\frac{1}{\eta _s}\overline{g\prime \eta } \right) $$

$$whx=\frac{1}{H_b}\frac{\overline{g\prime \eta }}{\eta _s}\frac{\partial H_b}{\partial y}=\frac{H_m}{H_b}\frac{\overline{g\prime \eta }}{\eta _s}\frac{\partial}{\partial y}\left( \frac{H_b}{H_m} \right) =\frac{\overline{g\prime \eta }}{\eta _{s}^{2}}\frac{\partial \eta _s}{\partial y}$$

$$why=\frac{1}{H_b}\frac{\overline{g\prime \eta }}{\eta _s}\frac{\partial H_b}{\partial x}=\frac{H_m}{H_b}\frac{\overline{g\prime \eta }}{\eta _s}\frac{\partial}{\partial x}\left( \frac{H_b}{H_m} \right) =\frac{\overline{g\prime \eta }}{\eta _{s}^{2}}\frac{\partial \eta _s}{\partial x}$$

```
call get_s2trit(s2t,rit,tracer%rict,dg,u,v,epsln)
```
**输入**
- ```u, v```：三维速度场（各层速度），来自动量模块
- ```rict```：来自 readyt 的稳定度相关量（已先算好）
- ```dg```：网格/度量信息（包括 odzt、a_cosa、a_sina、vit 等）
- ```epsln```：极小值，避免除零

**输出**
- ```s2t(i,j,k)```：垂向剪切的“平方项”（含网格度量修正）
- ```rit(i,j,k)```：基于 rict / (s2t+eps) 得到的稳定度指标（类似局地 Richardson 数或其变体）

该子程序首先计算

$$\left( \frac{\partial u}{\partial z} \right) ^2+\left( \frac{\partial v}{\partial z} \right) ^2$$

对应源码：
```Fortran
riv1  = u(i,j,k) - u (I,J,K +1)
riv2  = v(i,j,k) - v (I,J,K +1)
s2t (i,j,k) = (riv1 **2 + riv2**2 - 2.0D0 * riv1 * riv2 * dg%a_cosa(i,j))/dg%a_sina(i,j)/dg%a_sina(i,j) * dg%odzt(K+1) * dg%odzt(K+1) * dg%vit(i,j,k+1)
```
由于使用的是立方球网格，因此中间包含很多修正项。

其次计算Richardson数( $Ri$，手册P26)，

$$Ri=\frac{g}{\rho _0}\frac{\partial \rho _{pot}}{\partial z}/\left[ \left( \frac{\partial u}{\partial z} \right) ^2+\left( \frac{\partial v}{\partial z} \right) ^2 \right] $$

对应源码：
```Fortran
rit(i,j,k)= dg%vit(I,J,K+1)*rict(i,j,k)/(s2t(i,j,k)+epsln)
```

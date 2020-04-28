/*lyapnov指数計算用の線形方程式をルンゲ＝クッタ法で解く*/
#include <stdio.h>
#include <math.h>

/*ベクトルの構造体*/
typedef struct Vector3D 		
{
  double x, y, z;
} V3D; 

/*関数プロトタイプ宣言*/
double f(double, double);
double g(double, double, double);
double h(double, double, double);
double df(double, double);
double dg(double, double, double, double, double);
double dh(double, double, double, double, double);
V3D init(double, double, double);

/*parameter*/	//外部変数
int    sgm = 10;	
int    r   = 28;
double b   = 8.0/3.0;

main()
{
  /*変数宣言*/
  int i,j, n;											
  double x0, y0, z0, xn, yn, zn;
  V3D w0[3], wn[3];
  double t=0.0, dt; 
  double kf[4], kg[4], kh[4];
  double kdf[4], kdg[4], kdh[4];

  /*初期値の設定*/
  x0 = 1.0;  y0 = 1.0;  z0 = 1.0;						
  w0[0] = init( 1.0, 0, 0 );
  w0[1] = init( 0, 1.0, 0 );
  w0[2] = init( 0, 0, 1.0 );
  
  /*刻み幅・ルンゲ＝クッタ法の計算回数*/
  dt = 0.020;  n  = 200000;         					

  /*Runge-Kutta*/
  for(i=0; i < n; i++)
  {
    t+=dt;      

    kf[0] = dt * f(x0, y0);
    kg[0] = dt * g(x0, y0, z0);
    kh[0] = dt * h(x0, y0, z0);
    
    kf[1] = dt * f(x0+kf[0]/2.0, y0+kg[0]/2.0);
    kg[1] = dt * g(x0+kf[0]/2.0, y0+kg[0]/2.0, z0+kh[0]/2.0);
    kh[1] = dt * h(x0+kf[0]/2.0, y0+kg[0]/2.0, z0+kh[0]/2.0);
    
    kf[2] = dt * f(x0+kf[1]/2.0, y0+kg[1]/2.0);
    kg[2] = dt * g(x0+kf[1]/2.0, y0+kg[1]/2.0, z0+kh[1]/2.0);
    kh[2] = dt * h(x0+kf[1]/2.0, y0+kg[1]/2.0, z0+kh[1]/2.0);
    
    kf[3] = dt * f(x0+kf[2], y0+kg[2]);
    kg[3] = dt * g(x0+kf[2], y0+kg[2], z0+kh[2]);
    kh[3] = dt * h(x0+kf[2], y0+kg[2], z0+kh[2]);
    
    xn = x0 + (kf[0] + 2.0*(kf[1]+kf[2]) + kf[3])/6.0;  
    yn = y0 + (kg[0] + 2.0*(kg[1]+kg[2]) + kg[3])/6.0;
    zn = z0 + (kh[0] + 2.0*(kh[1]+kh[2]) + kh[3])/6.0;
   
    for(j=0; j < 3; j++)
    { 
      kdf[0] = dt * df(w0[j].x, w0[j].y);
	  kdg[0] = dt * dg(x0, z0, w0[j].x, w0[j].y, w0[j].z);
	  kdh[0] = dt * dh(x0, y0, w0[j].x, w0[j].y, w0[j].z);

	  kdf[1] = dt * df( w0[j].x + kdf[0]/2.0, w0[j].y + kdg[0]/2.0);
	  kdg[1] = dt * dg( x0, z0, w0[j].x + kdf[0]/2.0, w0[j].y + kdg[0]/2.0, w0[j].z + kdh[0]/2.0);
	  kdh[1] = dt * dh( x0, y0, w0[j].x + kdf[0]/2.0, w0[j].y + kdg[0]/2.0, w0[j].z + kdh[0]/2.0);

	  kdf[2] = dt * df( w0[j].x + kdf[1]/2.0, w0[j].y +kdg[1]/2.0);
	  kdg[2] = dt * dg( x0, z0, w0[j].x + kdf[1]/2.0, w0[j].y +kdg[1]/2.0, w0[j].z + kdh[1]/2.0);
	  kdh[2] = dt * dh( x0, y0, w0[j].x + kdf[1]/2.0, w0[j].y +kdg[1]/2.0, w0[j].z + kdh[1]/2.0);

	  kdf[3] = dt * df( w0[j].x + kdf[2], w0[j].y +kdg[2]);
	  kdg[3] = dt * dg( x0, z0, w0[j].x + kdf[2], w0[j].y +kdg[2], w0[j].z + kdh[2]);
	  kdh[3] = dt * dh( x0, y0, w0[j].x + kdf[2], w0[j].y +kdg[2], w0[j].z + kdh[2]);

      wn[j].x = w0[j].x + (kdf[0] + 2.0*(kdf[1]+kdf[2]) + kdf[3])/6.0;  
      wn[j].y = w0[j].y + (kdg[0] + 2.0*(kdg[1]+kdg[2]) + kdg[3])/6.0;
      wn[j].z = w0[j].z + (kdh[0] + 2.0*(kdh[1]+kdh[2]) + kdh[3])/6.0;
    } 

    x0 = xn; y0 = yn; z0 = zn;  /*次の点のために基準点の入れ替え*/ 

    w0[0] = init( wn[0].x, wn[0].y, wn[0].z );
    w0[1] = init( wn[1].x, wn[1].y, wn[1].z );
    w0[2] = init( wn[2].x, wn[2].y, wn[2].z );
  }
}

/*lorenz-eqの定義*/
double f(double x, double y) 
{
  return ( sgm*(-x + y) );
}

double g(double x, double y, double z)
{
  return ( r*x - y - x*z );
}

double h(double x, double y, double z)
{
  return ( x*y - b*z );
}  

/*領域軸を求める線形微分方程式の定義*/
double df(double w_x, double w_y)
{
  return ( sgm*(-w_x + w_y) );
}

double dg(double x, double z, double w_x, double w_y, double w_z)
{
  return ( (r - z)*w_x - w_y - x*w_z );
}

double dh(double x, double y, double w_x, double w_y, double w_z)
{
  return ( y*w_x + x*w_y - b*w_z );
}  

//ｗの初期化関数
V3D init(double a,double b,double c)
{
  V3D v;
   
  v.x = a; 
  v.y = b; 
  v.z = c; 

  return v;
}

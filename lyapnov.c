/*lyapnov指数計算用プログラム*/
#include <stdio.h>
#include <math.h>
#define NMAX 3

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
double inn_pro(V3D, V3D);
V3D init(double, double, double);
void gram_schmidt(V3D []); 
void sort(double [], V3D []);

/*parameter（外部変数）*/	
double sgm = 10.0;	
double r   = 28.0;
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
  double l0[3], ln[3], sum[3]={0.0,0.0,0.0}, lyap[3];

  /*初期値の設定*/
  x0 = 1.0;  y0 = 1.0;  z0 = 1.0;						
  w0[0] = init( 1.0, 0, 0 );
  w0[1] = init( 0, 1.0, 0 );
  w0[2] = init( 0, 0, 1.0 );
  
  /*刻み幅・ルンゲ＝クッタ法の計算回数*/
//dt = 0.020;  n  = 200000;	元のプログラムの数値	
//dt = 0.01;  n  = 1000;	精度が高い	
  dt = 0.01;  n  = 1000;

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

      wn[j].x = w0[j].x +  (kdf[0] + 2.0*(kdf[1]+kdf[2]) + kdf[3])/6.0;  
      wn[j].y = w0[j].y +  (kdg[0] + 2.0*(kdg[1]+kdg[2]) + kdg[3])/6.0;
      wn[j].z = w0[j].z +  (kdh[0] + 2.0*(kdh[1]+kdh[2]) + kdh[3])/6.0;

    } 

    gram_schmidt(wn); /*グラム・シュミット法で軸を直交化*/ 

    /*軸の規格化*/
    for(j=0; j < 3; j++)   
    { 
      ln[j] = sqrt( inn_pro(wn[j],wn[j]) );	
      wn[j] = init( (wn[j].x)/ln[j], (wn[j].y)/ln[j], (wn[j].z)/ln[j]);  
    }

    sort(ln,wn);   /*伸び率の大きい順に軸の添字をソート*/

    /*基準の入れ替え*/ 
    x0 = xn; y0 = yn; z0 = zn;  
    for(j=0; j < 3; j++)
    {
      //printf("wn[%d]:\t%lf\t%lf\t%lf\n", j+1, wn[j].x, wn[j].y, wn[j].z);
      sum[j] += log(ln[j]);
	  w0[j] = wn[j]; 
    }
  }

  /*リアプノフ指数の計算・結果の表示*/  
  for(i=0; i < 3; i++)
  {
	lyap[i] = sum[i] / (n * dt);
    printf("λ[%d]＝%lf\n", i+1, lyap[i]);
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

/*ｗの初期化関数*/
V3D init(double a, double b, double c)
{
  V3D v;
   
  v.x = a; 
  v.y = b; 
  v.z = c; 

  return v;
}

/*内積計算用の関数*/
double inn_pro(V3D a, V3D b)
{
  double sum = 0.0;		/*sumの初期化*/

  sum += a.x * b.x;		/*内積を成分で計算*/
  sum += a.y * b.y;
  sum += a.z * b.z;

  return sum;
}

/*gram-schmidt法で直交化*/
void gram_schmidt(V3D w0[])
{
  /*変数宣言*/
  int i,k;
  V3D wn[3];
  
  /*gram-schmidt*法で各軸を直交化*/
  for(i=0; i<3; i++)
  {
    wn[i] = w0[i];

    for(k=0; k < i; k++)
	{
	  wn[i].x -= ( inn_pro(wn[k],w0[i])/inn_pro(wn[k],wn[k]) ) * wn[k].x;
	  wn[i].y -= ( inn_pro(wn[k],w0[i])/inn_pro(wn[k],wn[k]) ) * wn[k].y;
	  wn[i].z -= ( inn_pro(wn[k],w0[i])/inn_pro(wn[k],wn[k]) ) * wn[k].z;
	}

    w0[i] = wn[i]; /*main関数に直交化した軸を返す*/
  }
}

/*大きい順にソート*/
void sort(double a[], V3D b[])
{
  /*変数宣言*/
  int i,j;
  int tmp_index; 
  double tmp_r;
  V3D tmp_v;

  /*数値の大小関係のチェック*/
  for(i=0; i<NMAX; i++)
  {
    tmp_r = a[i]; 
    tmp_v = b[i];    
    tmp_index = i;
    for(j=i+1; j<NMAX; j++)
    {
      if(a[j] > tmp_r)
        {
          tmp_r = a[j];
          tmp_v = b[j];
          tmp_index = j;
        }
    } 

    /*入れ替え操作*/
    if(i != tmp_index) 
    {
      a[tmp_index] = a[i]; 
      a[i] = tmp_r;
      b[tmp_index] = b[i]; 
      b[i] = tmp_v;
    } 
  }
}

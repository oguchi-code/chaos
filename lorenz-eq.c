/*ローレンツ方程式の解を数値計算*/
/*データをlorenz.datというファイルに出力する*/
#include <stdio.h>

/*parameter*/
int    sgm = 10;
int    r   = 28;
double b   = 8.0/3.0;

/*Lorenz-eqの定義*/
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

int main(void)
{ 
  int i, n;						            /*変数宣言*/
  double x0, y0, z0, xn, yn, zn;
  double t=0.0, dt, kf[4], kg[4], kh[4];
  FILE *fp, *gp;					        /*ストリームポインタ*/
  
  x0 = 1.0;  y0 = 1.0;  z0 = 1.0;			/*初期値の設定*/
  dt = 0.020;  n  = 200000;         	    /*刻み幅・ルンゲ＝クッタ法の計算回数*/

  fp = fopen("lorenz.dat", "w");			/*ファイルのオープン*/
  gp = popen("gnuplot -persist", "w");      /*gnuplotの呼び出し(残留)*/
  
  /*Runge-Kutta*/
  for(i = 0 ; i <= n ; i++)
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
    
    fprintf(fp, "%lf\t%lf\t%lf\t%lf\t\n", t, xn, yn, zn); /*fpのファイルにデータを出力*/
    
    x0 = xn; y0 = yn; z0 = zn;                            /*次点のために基準の入替*/ 
  }
  
  fclose(fp);						                      /*ファイルのクローズ*/
  fprintf(gp, "splot \"lorenz.dat\" u 2:3:4 w l\n");      /*gnuplotに指示*/		
  pclose(gp);						                      /*gnuplotの終了*/

  return 0;
}

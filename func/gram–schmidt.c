#include <stdio.h>
#include <math.h>

/*ベクトルの構造体*/
typedef struct Vector3D 		
{
  double x, y, z;
} V3D; 

/*関数の実装*/
double inn_pro(V3D,V3D);
V3D init(double,double,double);

/*gram-schmidt法で正規直交化*/
main()
{
  /*変数宣言*/
  int i, k;
  double c;
  V3D w[3],wnew[3];
  
  /*初期化*/
  w[0] = init( 2.0, 0, 0);
  w[1] = init( 0, 1.0, 0);
  w[2] = init( 0, 0, 3.0);


  for(i=0; i<3; i++)
  {
    wnew[i] = w[i];
    for(k=0; k < i; k++)
	{
	  wnew[i].x -= ( inn_pro(wnew[k],w[i])/inn_pro(wnew[k],wnew[k]) ) * wnew[k].x;
	  wnew[i].y -= ( inn_pro(wnew[k],w[i])/inn_pro(wnew[k],wnew[k]) ) * wnew[k].y;
	  wnew[i].z -= ( inn_pro(wnew[k],w[i])/inn_pro(wnew[k],wnew[k]) ) * wnew[k].z;
	}

    c = sqrt( inn_pro(wnew[i],wnew[i]) );

	wnew[i].x /= c;  
	wnew[i].y /= c;
	wnew[i].z /= c;
  }

  for(i=0; i<3; i++)	
  {
    printf("wnew [%d] の x 成分:%f\n",i+1,wnew[i].x);
    printf("wnew [%d] の y 成分:%f\n",i+1,wnew[i].y);
    printf("wnew [%d] の z 成分:%f\n",i+1,wnew[i].z);
  }
}

/*内積計算用の関数*/
double inn_pro(V3D a,V3D b)
{
  double sum = 0.0;		/*sumの初期化*/

  sum += a.x * b.x;		/*内積を成分で計算*/
  sum += a.y * b.y;
  sum += a.z * b.z;

  return sum;
}

//構造体初期化関数
V3D init(double a,double b,double c)
{
  V3D v;
   
  v.x = a; 
  v.y = b; 
  v.z = c; 

  return v;
}

/*コールバック関数*/
/*関数の引数に使われる関数*/
/*ルンゲ・クッタ法を使用する関数を変更する方法*/
/*関数のポインタ*/
#include<stdio.h>

/*型宣言*/
typedef double (*Ftmp)(double x, double y);

/*関数プロトタイプ宣言*/
double add (double, double);
double sub (double, double);

/*main program*/
main()
{
  double x=10, y=20;
  int i;

  Ftmp func[2]={ add, sub };
/*  same results  */
/*  Ftmp func[2]={ &add, &sum }; */  
/*	func[0]=&add; func[1]=&sub;	 */


  for(i=0; i<2; i++)
    {
      printf("%d, answer, %f\n", i,func[i](x,y)); 
    }
}

/*	functions  */
double add(double x, double y)
{
	return x+y;
}

double sub(double x, double y)
{
	return x-y;
}

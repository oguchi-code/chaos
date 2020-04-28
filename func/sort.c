/*大きい順にソート*/
#include <stdio.h>
#define NMAX 5

main()
{
  /*変数宣言*/
  int i,j;
  double a[NMAX]; 
  int tmp, tmp_index; 

  /*数値の入力*/
  printf("sortしたい数値を入力してくだちい!\n");

  for(i=0; i<NMAX; i++)
  {
    scanf("%lf", a+i);
  } 

  /*数値の大小関係のチェック*/
  for(i=0; i<NMAX; i++)
  {
    tmp = a[i]; 
    tmp_index = i;
    for(j=i+1; j<NMAX; j++)
    {
      if(a[j] > tmp)
        {
          tmp = a[j];
          tmp_index = j;
        }
    } 

    /*入れ替え操作*/
    if(i != tmp_index) 
    {
      a[tmp_index] = a[i]; 
      a[i] = tmp;
    } 
  }

  /*結果の表示*/
  for(i=0; i<NMAX; i++)
  {
    printf("%lf ", a[i]);
  }
  printf("\n");
}

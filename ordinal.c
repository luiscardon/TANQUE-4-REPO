#include<stdio.h>
#include<stdlib.h>

int ordinal(int dd, int mm, int yy);

void main()
{
int dd,mm,yy,julian;
dd=16;
mm=3;
yy=17;
julian=ordinal(dd,mm,yy);

/*printf("\nenter a valid date (dd mm yy)");
scanf("%d%d%d",&dd,&mm,&yy);
julian=dd;
mm=mm-1;
switch(mm)
{
case 11 : julian=julian+30;
case 10 : julian=julian+31;
case 9 : julian=julian+30;
case 8 : julian=julian+31;
case 7 : julian=julian+31;
case 6: julian=julian+30;
case 5: julian=julian+31;
case 4 : julian=julian+30;
case 3 : julian=julian+31;
case 2 : if(yy%4==0)julian=julian+29;
             else julian=julian+28;
case 1 : julian=julian+31;
}
*/
printf("\nJulian date is %d\n",julian);

  
}

int ordinal(int dd, int mm, int yy)
{
int julian;

julian=dd;
mm=mm-1;
switch(mm)
{
case 11 : julian=julian+30;
case 10 : julian=julian+31;
case 9 : julian=julian+30;
case 8 : julian=julian+31;
case 7 : julian=julian+31;
case 6: julian=julian+30;
case 5: julian=julian+31;
case 4 : julian=julian+30;
case 3 : julian=julian+31;
case 2 : if(yy%4==0)julian=julian+29;
             else julian=julian+28;
case 1 : julian=julian+31;
}
return julian;

  
}
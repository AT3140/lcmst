#include<stdio.h>
#include<stdlib.h>

void fun(int *arr){
    //static int count =0;
    //count++;
    if(arr)
    printf("%d",arr[2]);
    else printf("hello");
}

int main(){
    //static int count =0;
    int arr[3]={3,2,5};
    fun(NULL);

    //printf("%d",a);
    return 0;
}
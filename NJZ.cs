using System;

class operator_i
{
    public Complex[,] bcinv(Complex[,] cpx){//原矩阵 求逆函数
        int flag,n=cpx.Rank;;//flag判断奇异性；n是矩阵维度
        Complex[,] _bcinv=new Complex[n,cpx.GetLength(0)];//逆矩阵
        double [,] ar=new double[n,n];//实部矩阵ar
        double [,] ai=new double[n,n];//虚部矩阵ai
        double d,p,t,q,s,b;//中间变量
        int[] _is = new int[n],js=new int[n];// 中间变量

        //n=cpx.Rank;
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                ar[i,j] = cpx[i,j].Real;
                ai[i,j] = cpx[i,j].Imaginary;
            }
        }
        flag=1;
        for(int k=0;k<n;k++){
            d=0.0;
            for(int i=k;i<n;i++){
                for(int j=k;j<n;j++){
                    p=ar[i,j]*ar[i,j]+ai[i,j]*ai[i,j];
                    if(p>d){
                        d=p;
                        _is[k]=i;
                        js[k]=j;
                    }
                }
            }
            if(d+1.0==1.0){
                flag=0;
                Console.WriteLine("flag=0,复矩阵奇异！" );
                return _bcinv;
            }
            for(int j=0;j<n;j++){
                t=ar[k,j];
                ar[k,j]=ar[_is[k],j];
                ar[_is[k],j]=t;
                t=ai[k,j];
                ai[k,j]=ai[_is[k],j];
                ai[_is[k],j]=t;
            }
            for(int i=0;i<n;i++){
                t=ar[i,k];
                ar[i,k]=ar[i,js[k]];
                ar[i,js[k]]=t;
                t=ai[i,k];
                ai[i,k]=ai[i,js[k]];
                ai[i,js[k]]=t;
            }
            ar[k,k]=ar[k,k]/d;
            ai[k,k]=-ai[k,k]/d;
            for(int j=0;j<n;j++){
                if(j!=k){
                    p=ar[k,j]*ar[k,k];
                    q=ai[k,j]*ai[k,k];
                    s=(ar[k,j]+ai[k,j])*(ar[k,k]+ai[k,k]);
                    ar[k,j]=p-q;
                    ai[k,j]=s-p-q;
                }
            }
            for(int i=0;i<n;i++){
                if(i!=k){
                    for(int j=0;j<n;j++){
                        if(j!=k){
                            p=ar[k,j]*ar[i,k];
                            q=ai[k,j]*ai[i,k];
                            s=(ar[k,j]+ai[k,j])*(ar[i,k]+ai[i,k]);
                            t=p-q;
                            b=s-p-q;
                            ar[i,j]=ar[i,j]-t;
                            ai[i,j]=ai[i,j]-b;
                        }
                    }
                }
            }
            for(int i=0;i<n;i++){
                if(i!=k){
                    p=ar[i,k]*ar[k,k];
                    q=ai[i,k]*ai[k,k];
                    s=(ar[i,k]+ai[i,k])*(ar[k,k]+ai[k,k]);
                    ar[i,k]=q-p;
                    ai[i,k]=p+q-s;
                }
            }
        }
        for(int k=n-1;k>=0;k--){
            for(int j=0;j<n;j++){
                t=ar[k,j];
                ar[k,j]=ar[js[k],j];
                ar[js[k],j]=t;
                t=ai[k,j];
                ai[k,j]=ai[js[k],j];
                ai[js[k],j]=t;
            }
            for(int i=0;i<n;i++){
                t=ar[i,k];
                ar[i,k]=ar[i,_is[k]];
                ar[i,_is[k]]=t;
                t=ai[i,k];
                ai[i,k]=ai[i,_is[k]];
                ai[i,_is[k]]=t;
            }
        }
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                _bcinv[i,j].Real=ar[i,j];
                _bcinv[i,j].Imaginary=ai[i,j];
            }
        }
        return _bcinv;
    }
    
}
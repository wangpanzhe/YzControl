/*
C    该程序取自 《FORTRAN 常用算法程序集》第二版 徐士良
C    参数说明：
C     N：整形，第一类贝塞尔函数的阶数，要求N>=0.当N<0时，按|N|
C     X：双精度实型变量，自变量值。
C     函数名MNSL1返回双精度实型函数值JN(X) // N,下角标
C
C    调用时 写为Y=MBSL1(N,X) 
*/
using System;
class FUNC_MBSL1
{
    public double MBSL1(int N,double X){
        double mbsl1;
        double T,Y,Z,P,Q,S,B0,B1;
        double[] A=new double[6] {57568490574.0,-13362590354.0,651619640.7, -11214424.18,77392.33017,-184.9052456};
        double[] B=new double[6] {57568490411.0,1029532985.0,9494680.718,59272.64853,267.8532712,1.0};
        double[] C=new double[6] {72362614232.0,-7895059235.0,242396853.1,-2972611.439,15704.4826,-30.16036606};
        double[] D=new double[6] {144725228443.0,2300535178.0,18583304.74,99447.43394,376.9991397,1.0};
        double[] E=new double[5] {1.0,-0.1098628627D-02,0.2734510407D-04,-0.2073370639D-05,0.2093887211D-06};
        double[] F=new double[5] {-0.1562499995D-01,0.1430488765D-03,-0.6911147651D-05,0.7621095161D-06,-0.934935152D-07};
        double[] G=new double[5] {1.0,0.183105D-02,-0.3516396496D-04,0.2457520174D-05,-0.240337019D-06};
        double[] H=new double[5] {0.4687499995D-01,-0.2002690873D-03,0.8449199096D-05,-0.88228987D-06,0.105787412D-06};

        T=Math.Abs(X);
        if(N<0) N=-N;
        if(N!=1){
            if(T<8.0){
                Y=T*T;
                P=A[5];
                Q=B[5];//下标修改
                for(int i=4;i>=0;i--){
                    P=P*Y+A[i];
	                Q=Q*Y+B[i];
                }
                P=P/Q;
            }
            else{
                Z=8.0/T;
	            Y=Z*Z;
                P=E[4];
	            Q=F[4];//下标修改
                for(int i=3;i>=0;i--){
                    P=P*Y+E[i];
	                Q=Q*Y+F[i];
                }
                S=T-0.785398164;
                P=P*Math.Cos(S)-Z*Q*Math.Sin(S); //----------角度弧度的对应。
	            P=P*Math.Sqrt(0.636619772/T);
            }
        }
        if(N==0){
            mbsl1=P;//-------------初始值未确定
            return mbsl1;
        }
        B0=P;
        if(T<8.0){
            Y=T*T;
	        P=C[5];
	        Q=D[5];//下标修改
            for(int i=4;i>=0;i--){
                P=P*Y+C[i];
	            Q=Q*Y+D[i];
            }
            P=X*P/Q;
        }
        else{
            Z=8.0/T;
	        Y=Z*Z;
	        P=G[4];
	        Q=H[4];
            for(int i=3;i>=0;i--){
                P=P*Y+G[i];
	            Q=Q*Y+H[i];
            }
            S=T-2.356194491;
	        P=P*Math.Cos(S)-Z*Q*Math.Sin(S);
	        P=P*X*Math.Sqrt(0.636619772/T)/T;
        }
        if(N==1){
            mbsl1=P;
            return mbsl1;
        }
        B1=P;
        if(X==0.0){
            mbsl1=0.0;
            return mbsl1;
        }
        S=2.0/T;
        if(T>1.0*N){
            if(X<0.0) B1=-B1;
            for(int i=1;i<=N-1;i++){
                P=S*i*B1-B0;
	            B0=B1;
	            B1=P;
            }
        }
        else{
            int M=Convert.ToInt32(Math.Sqrt(40.0*N));//修改类型
	        M=(N+M)/2;
	        M=M+M;
	        P=0.0;
	        Q=0.0;
	        B0=1.0;
	        B1=0.0;
            for(int i=M-1;i>=0;i--){
                T=S*(i+1)*B0-B1;
	            B1=B0;
	            B0=T;
                if(Math.Abs(B0)>1.0D+10){
                    B0=B0*1.0D-10;
	                B1=B1*1.0D-10;
	                P=P*1.0D-10;
	                Q=Q*1.0D-10;
                }
                if((i+2)%2==0) Q=Q+B0;
                if(i+1==N) P=B1;
            }
            Q=2.0*Q-B0;
	        P=P/Q;
        }
        if((X<0.0)&&(N%2==1)) P=-P;
        mbsl1=P;
        return mbsl1;
    }
}
/*
该程序取自 《FORTRAN 常用算法程序集》第二版 徐士良
C  参数说明：
C    PR、PI，双精度实型一位数组，长度为N,输入兼输出参数。
C        分别为输入序列的实部和虚部;同时也返回傅里叶变换后的模和幅角
C    N 整型，输入的点数，要求N=2^K(K>0)
C    K 整型，要求满足N=2^K(K>0)
C    FR,FI,双精度实型一维数组，长度为N,输出参数。
C         分别存放计算后的实部和虚部
C    L=0，表示傅里叶变换；L=1,表示逆傅里叶变换
C    IL=0,表示计算不Fourier或逆Fourier变换的模和幅角；
C        IL=1,表示计算
*/

using System;

class SUB_KKFFT
{
    public void KKFFT(double[] PR,double[] PI,int N,int K,double[] FR,double[] FI,int L,int IL){
        double  P,Q,S,VR,VI,PODDR,PODDI;

        int M,IS,J;
        for(int it=0;it<=N-1;it++){
            M=it;
            IS=0;
            for(int i=0;i<=K-1;i++){
                J=M/2;
	            IS=2*IS+(M-2*J);
	            M=J;
            }//10
            FR[it]=PR[IS];
	        FI[it]=PI[IS];
        }//20

        PR[0]=1.0;
	    PI[0]=0.0;
	    PR[1]=Math.Cos(6.283185306/N);
	    PI[1]=-Math.Sin(6.283185306/N);
	    if(L!=0) PI[1]=-PI[1];
        for(int i=2;i<N;i++){
            P=PR[i-1]*PR[1];
	        Q=PI[i-1]*PI[1];
	        S=(PR[i-1]+PI[i-1])*(PR[1]+PI[1]);
	        PR[i]=P-Q;
	        PI[i]=S-P-Q;
        }//30

        for(int it=-1;it<N-2;it+=2){
            VR=FR[it+1];
	        VI=FI[it+1];
	        FR[it+1]=VR+FR[it+2];
	        FI[it+1]=VI+FI[it+2];
	        FR[it+2]=VR-FR[it+2];
	        FI[it+2]=VI-FI[it+2];
        }//40
        M=N/2;
	    int NV=2;
        for(int l0=K-2;l0>=0;l0--){
            M=M/2;
	        NV=2*NV;
            for(int it=0;it<=(M-1)*NV;it=it+NV){
                for(int j=0;j<=(NV/2)-1;j++){
                    P=PR[M*j]*FR[it+j+NV/2];
	                Q=PI[M*j]*FI[it+j+NV/2];
	                S=PR[M*j]+PI[M*j];
	                S=S*(FR[it+j+NV/2]+FI[it+j+NV/2]);
	                PODDR=P-Q;
	                PODDI=S-P-Q;
	                FR[it+j+NV/2]=FR[it+j]-PODDR;
	                FI[it+j+NV/2]=FI[it+j]-PODDI;
	                FR[it+j]=FR[it+j]+PODDR;
	                FI[it+j]=FI[it+j]+PODDI;
                }
            }//60
        }//70

        if(L!=0){
            for(int i=0;i<N;i++){
                FR[i]=FR[i]/N;
	            FI[i]=FI[i]/N;
            }
        }
        if(IL!=0){
            for(int i=0;i<N;i++){
                PR[i]=Math.Sqrt(FR[i]*FR[i]+FI[i]*FI[i]);
	            PI[i]=Math.Atan(FI[i]/FR[i])*360.0/6.283185306;
            }//90
        }
    }
}
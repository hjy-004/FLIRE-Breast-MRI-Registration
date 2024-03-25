#include <math.h>
#include <mex.h>

/* ------------------------------------------------------------  */
/*                   FOR WINDOWS PORTABILITY                     */
/* ------------------------------------------------------------  */
#ifdef _WIN32
#include <float.h>
#define isnan _isnan
#define finite _finite
#define M_PI 3.141592
bool isinf(double x)
{
  if (_finite(x)==0)
    return false;
  else
    return true;	
}

int round(double val)
{
    double valf=floor(val);
    if (val>0)
    {
    	if ((val-valf)>=0.5)
    	 return (int)(valf+1);
      else
    	return (int)(valf);
    }
    else
    {
      if ((val-valf)>=0.5)
    	return (int)(valf);
      else
    	return (int)(valf+1);
    }
}
#endif
/* ------------------------------------------------------------  */

//Matlab is column indexing
int index(const int *ps, int i,int j, int k)
{return (k-1)*ps[1]*ps[0]+(j-1)*ps[0]+i-1;}

void ind2sub(const int *ps, const int ind, int &i, int &j, int &k)
{
    int ndx;
    ndx = ind;
    k = (int)floor(ndx / (ps[1]*ps[0])) +1;
    ndx = (int)fmod(ndx, ps[1]*ps[0]);
    
    j= (int)floor(ndx / ps[0] ) +1;
    ndx = (int)fmod(ndx, ps[0]);
     
    i=ndx+1;
}

bool inbound(const double *pvxl, const int *ps, int interpm, int pad)
{
    
    int iL=0, iH=0, jL=0, jH=0, kL=0, kH=0;
    
    switch(interpm)
    {
        case 0: // Nearest Neighbor
            iL=(int)round(pvxl[0]);
            iH=iL;
            jL=(int)round(pvxl[1]);
            jH=jL;
            kL=(int)round(pvxl[2]);
            kH=kL;
            break;
        case 1: //Linear Interpolation
            iL=(int)(pvxl[0]);
            iH=iL+1;
            jL=(int)(pvxl[1]);
            jH=jL+1;
            kL=(int)(pvxl[2]);
            kH=kL+1;
            break;
        case 2: // Cubic Interpolation
        case 3: // Key's function
            iL=(int)(pvxl[0])-1;
            iH=iL+3;
            jL=(int)(pvxl[1])-1;
            jH=jL+3;
            kL=(int)(pvxl[2])-1;
            kH=kL+3;
            break;
        case 4: // Cubic Spline interpolation
        case 5: // Sinc Interpolation
            iL=(int)(pvxl[0])+1-pad;
            iH=(int)(pvxl[0])+pad;
            jL=(int)(pvxl[1])+1-pad;
            jH=(int)(pvxl[1])+pad;
            kL=(int)(pvxl[2])+1-pad;
            kH=(int)(pvxl[2])+pad;            
            break;
    }
	// Using Matlab 1 based
	if (iH>ps[0])
		return false;

	if (jH>ps[1])
		return false;

	if (kH>ps[2])
		return false;

	 if (iL<1)
		return false;
 
	if (jL<1)
		return false;
 
	 if (kL<1)
		return false;	

	 return true;
}

double getVxlVal_Linear(const int *ps, const double *vxli, const double *pimg)
{
    //Linear interpolation
    int i;
    int im=(int)vxli[0];
    int ip=im+1;
    float tim=vxli[0]-im;
    float tip=1-tim;
    
    int jm=(int)vxli[1];
    int jp=jm+1;
    float tjm=vxli[1]-jm;
    float tjp=1-tjm;
    
    int km=(int)vxli[2];
    int kp=km+1;
    float tkm=vxli[2]-km;
    float tkp=1-tkm;
    
    double tf1[4],tf2[4];
    tf1[0]=pimg[index(ps,im,jm,km)];
    tf1[1]=pimg[index(ps,ip,jm,km)];
    tf1[2]=pimg[index(ps,ip,jp,km)];
    tf1[3]=pimg[index(ps,im,jp,km)];
    tf2[0]=pimg[index(ps,im,jm,kp)];
    tf2[1]=pimg[index(ps,ip,jm,kp)];
    tf2[2]=pimg[index(ps,ip,jp,kp)];
    tf2[3]=pimg[index(ps,im,jp,kp)];
    //fVec4 tf1(f0,f1,f2,f3), tf2(f4,f5,f6,f7);
    double tf3[4],t3[4],t4[4];
    for (i=0;i<4;i++)
        tf3[i]=tkp*tf1[i]+tkm*tf2[i];
    //fVec4 tf3=tf1*tkp+tf2*tkm;
    //fVec4 t3(tip, tim, tim, tip);
    t3[0]=tip;t3[1]=tim;t3[2]=tim;t3[3]=tip;
    //fVec4 t4(tjp, tjp, tjm, tjm);
    t4[0]=tjp;t4[1]=tjp;t4[2]=tjm;t4[3]=tjm;
    double val=0.;
    for (i=0;i<4;i++)
        val+=tf3[i]*t3[i]*t4[i];
    //fVec4 sumt=tf3*t3*t4;
    //val= (sumt.sum());
	return val;
}

// Cubic Interpolation
double cubicInterp(double x, const double *y)
{
    double w0=y[1];
    double w1=(-y[0]+y[2])/2.;
    double w2=(2*y[0]-5*y[1]+4*y[2]-y[3])/2.;
    double w3=(-y[0]+3*y[1]-3*y[2]+y[3])/2.;
//      double w0=y[1];
//      double w1=(-2*y[0]-3*y[1]+6*y[2]-y[3])/6.;
//      double w2=(3*y[0]-6*y[1]+3*y[2])/6.;
//      double w3=(-y[0]+3*y[1]-3*y[2]+y[3])/6.;
    return  w0+w1*x+w2*x*x+w3*x*x*x;
}

double getVxlVal_Cubic(const int *ps, const double *vxli, const double *pimg)
{
    int i,j,k;
    int im=(int)vxli[0];
    int ip=im+1;
    float tim=vxli[0]-im;
    float tip=1-tim;
    
    int jm=(int)vxli[1];
    int jp=jm+1;
    float tjm=vxli[1]-jm;
    float tjp=1-tjm;
    
    int km=(int)vxli[2];
    int kp=km+1;
    float tkm=vxli[2]-km;
    float tkp=1-tkm;
    
    double fj[4], fi[4], fk[4];
    for (k=km-1;k<=km+2;k++)
    {
        for (i=im-1;i<=im+2;i++)
        {
            for (j=jm-1;j<=jm+2;j++)
                fj[j-jm+1]=pimg[index(ps,i,j,k)];
            fi[i-im+1]=cubicInterp(tjm, fj);
        }
        fk[k-km+1]=cubicInterp(tim, fi);
    }
    return cubicInterp(tkm, fk);
}

// Key function interpolation
double keyu(double x)
{
    double ax=fabs(x);
    if ((ax>=0) && (ax<1))
        return 1.5*ax*ax*ax-2.5*ax*ax+1.;
    if ((ax>=1) && (ax<2))
        return -0.5*ax*ax*ax+2.5*ax*ax-4*ax+2.;
    return 0.;
}

double keyInterp(double x, const double *y)
{
    double w0=keyu(-1.-x);
    double w1=keyu(-x);
    double w2=keyu(1-x);
    double w3=keyu(2-x);
    return w0*y[0]+w1*y[1]+w2*y[2]+w3*y[3];
}

double getVxlVal_Key(const int *ps, const double *vxli, const double *pimg)
{
     int i,j,k;
    int im=(int)vxli[0];
    int ip=im+1;
    float tim=vxli[0]-im;
    float tip=1-tim;
    
    int jm=(int)vxli[1];
    int jp=jm+1;
    float tjm=vxli[1]-jm;
    float tjp=1-tjm;
    
    int km=(int)vxli[2];
    int kp=km+1;
    float tkm=vxli[2]-km;
    float tkp=1-tkm;
    
    double fj[4], fi[4], fk[4];
    for (k=km-1;k<=km+2;k++)
    {
        for (i=im-1;i<=im+2;i++)
        {
            for (j=jm-1;j<=jm+2;j++)
                fj[j-jm+1]=pimg[index(ps,i,j,k)];
            fi[i-im+1]=keyInterp(tjm, fj);
        }
        fk[k-km+1]=keyInterp(tim, fi);
    }
    return keyInterp(tkm, fk);
}


// Cubic Spline interpolation from numerical recipe
void spline(const double *x, const double *y, int n , double yp1, double ypn, double *y2)
{
    
    int i,k;
    double p, qn, sig, un;
    double *u=(double *)mxMalloc(sizeof(double)*(n-1));
    
    if (yp1>0.99e30)
        y2[0]=u[0]=0.0;
    else
    {
        y2[0]=-0.5;
        u[1]=(3./(x[2]-x[1]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
    }
    
    for (i=1;i<n-1;i++)
    {
        sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
        p=sig*y2[i-1]+2.0;
        y2[i]=(sig-1.0)/p;
        u[i]=(y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]);
        u[i]=(6.*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }
    
    if (ypn>0.99e30)
        qn=un=0.;
    else
    {
        qn=0.5;
        un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
    }
    y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
    
    for (k=n-2;k>=0;k--)
        y2[k]=y2[k]*y2[k+1]+u[k];
    
    mxFree(u);
}

void splint(const double *xa, const double *ya, const double *y2a, int n, int m, const double *x, double *y)
{
    int klo, khi, k;
    double h,b,a;
    for (int i=0;i<m;i++)
    {
        klo=0;
        khi=n-1;
        while((khi-klo)>1)
        {
            k=(khi+klo) >>1;
            if (xa[k]>x[i])
                khi=k;
            else
                klo=k;
        }
        
        h=xa[khi]-xa[klo];
        a=(xa[khi]-x[i])/h;
        b=(x[i]-xa[klo])/h;
        y[i]=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
    }
}

double getVxlVal_CSpline(const int *ps, const double *vxli, const double *pimg, int padding)
{
     int i,j,k;
    int im=(int)vxli[0];
    int ip=im+1;
    float tim=vxli[0]-im;
    float tip=1-tim;
    
    int jm=(int)vxli[1];
    int jp=jm+1;
    float tjm=vxli[1]-jm;
    float tjp=1-tjm;
    
    int km=(int)vxli[2];
    int kp=km+1;
    float tkm=vxli[2]-km;
    float tkp=1-tkm;
    
    
    double * xrow=(double *)mxMalloc(sizeof(double)*(2*padding));
    double * xcol=(double *)mxMalloc(sizeof(double)*(2*padding));
    double * xdep=(double *)mxMalloc(sizeof(double)*(2*padding));
    double *fi=(double *)mxMalloc(sizeof(double)*(2*padding));
    double *f2i=(double *)mxMalloc(sizeof(double)*(2*padding));
    double *fj=(double *)mxMalloc(sizeof(double)*(2*padding));
    double *f2j=(double *)mxMalloc(sizeof(double)*(2*padding));
    double *fk=(double *)mxMalloc(sizeof(double)*(2*padding));
    double *f2k=(double *)mxMalloc(sizeof(double)*(2*padding));
    
    double val=0.;
    for (k=kp-padding ;k<=kp+padding-1;k++)
    {
        xdep[k-(kp-padding)]=k;
        for (i=ip-padding;i<=ip+padding-1;i++)
        {
            xrow[i-(ip-padding)]=i;
            for (j=jp-padding;j<=jp+padding-1;j++)
            {
                xcol[j-(jp-padding)]=j;
                fj[j-(jp-padding)]=pimg[index(ps,i,j,k)];
            }
            spline(xcol, fj, 2*padding, 1e30, 1e30, f2j);
            splint(xcol, fj, f2j, 2*padding, 1, &vxli[1], &fi[i-(ip-padding)]);
        }
        spline(xrow, fi, 2*padding, 1e30, 1e30, f2i);
        splint(xrow, fi, f2i, 2*padding, 1, &vxli[0], &fk[k-(kp-padding)]);
    }
    
    spline(xdep, fk, 2*padding, 1e30, 1e30, f2k);
    splint(xdep, fk, f2k, 2*padding, 1, &vxli[2], &val);
    
    mxFree(fi);
    mxFree(fj);
    mxFree(fk);
    mxFree(f2i);
    mxFree(f2j);
    mxFree(f2k);
    mxFree(xrow);
    mxFree(xcol);
    mxFree(xdep);
    
    return val;
}

// SINC interpolation
double hamm_sinc(double x, double fullwidth)
{
    if (fabs(x)<1e-5)
        return 1.;
    else
    {
        double val=sin(M_PI*x)/(M_PI*x);
        val*= 0.54 + 0.46 * cos(2.0*M_PI*x/fullwidth);
        return val;
    }
}

double sincInterp(const double *x, const double *y, const double xi, const int padding)
{
    double val=0.;
    double w=0.;
    double sum_w=0.;
    for (int i=0;i<=2*padding-1;i++)
    {
        w=hamm_sinc(xi-x[i], 2*padding);
        sum_w=sum_w+w;
        val=val+w*y[i];
    }
    return val/sum_w;
}

double getVxlVal_Sinc(const int *ps, const double *vxli, const double *pimg, int padding)
{
     int i,j,k;
    int im=(int)vxli[0];
    int ip=im+1;
    float tim=vxli[0]-im;
    float tip=1-tim;
    
    int jm=(int)vxli[1];
    int jp=jm+1;
    float tjm=vxli[1]-jm;
    float tjp=1-tjm;
    
    int km=(int)vxli[2];
    int kp=km+1;
    float tkm=vxli[2]-km;
    float tkp=1-tkm;
    
    
    double * xrow=(double *)mxMalloc(sizeof(double)*(2*padding));
    double * xcol=(double *)mxMalloc(sizeof(double)*(2*padding));
    double * xdep=(double *)mxMalloc(sizeof(double)*(2*padding));
    double *fi=(double *)mxMalloc(sizeof(double)*(2*padding));
    double *fj=(double *)mxMalloc(sizeof(double)*(2*padding));
    double *fk=(double *)mxMalloc(sizeof(double)*(2*padding));
    
    double val=0.;
    for (k=kp-padding ;k<=kp+padding-1;k++)
    {
        xdep[k-(kp-padding)]=k;
        for (i=ip-padding;i<=ip+padding-1;i++)
        {
             xrow[i-(ip-padding)]=i;
             for (j=jp-padding;j<=jp+padding-1;j++)
             {
                 xcol[j-(jp-padding)]=j;
                 fj[j-(jp-padding)]=pimg[index(ps,i,j,k)];
             }
           fi[i-(ip-padding)] = sincInterp(xcol, fj, vxli[1], padding);
        }
        fk[k-(kp-padding)] = sincInterp(xrow, fi, vxli[0], padding);
    }
    val=sincInterp(xdep, fk, vxli[2], padding);
    
    mxFree(fi);
    mxFree(fj);
    mxFree(fk);
    mxFree(xrow);
    mxFree(xcol);
    mxFree(xdep);
    
    return val;
}

bool getvxlval(const double *pvxlo, const int *ps, const double *pimg, const double *pmv2v, int interpm, int padding, double &val)
{
	double vxli[4];
	int	i,j,k;
    val = 0.;
	//mapping to voli vxl coordinate
	// becareful about column indxing in matlab
	for (i=0;i<4;i++)
	{
		vxli[i]=0;
		for (j=0;j<4;j++)
			vxli[i]=vxli[i]+pmv2v[j*4+i]*pvxlo[j];
	}

    bool bInBound= inbound(vxli, ps, interpm, padding);
	if (bInBound)
	{
        switch(interpm)
        {
            case 0: // Nearest Neighbor
                val=pimg[index(ps, (int)round(vxli[0]), (int)round(vxli[1]), (int)round(vxli[2]))];
                break;
            case 1: //Linear Interpolation
                val = getVxlVal_Linear(ps, vxli, pimg);  
                break;
            case 2: // Cubic Interpolation
                val = getVxlVal_Cubic(ps, vxli, pimg); 
                break;
            case 3: // Key's function
                val = getVxlVal_Key(ps, vxli, pimg);
                break;       
            case 4: // Cubic Spline interpolation
                val= getVxlVal_CSpline(ps, vxli, pimg, padding);
                break;
            case 5: // Sinc Interpolation
                val=getVxlVal_Sinc(ps, vxli, pimg, padding);
                break;
        }
	}
	
    return bInBound; 
}

double getmin(int isize, const double *pval)
{
	double minv=1000000.;
	for (int i=0;i<isize;i++)
	{
		if (pval[i]<minv)
			minv=pval[i];
	}
	return minv;

}

double getmax(int isize, const double *pval)
{
	double maxv=-1000000.;
	for (int i=0;i<isize;i++)
	{
		if (pval[i]>maxv)
			maxv=pval[i];
	}
	return maxv;

}


int binind1d(double dbr, double dmin, int numbins, double val, double thmin, double thmax)
{
	if ((val>thmin) && (val<thmax))
	{
		int tmp=(int)((val-dmin)/dbr*0.99999);
		if (tmp<0)
			return 0;
		if (tmp>=numbins)
			return numbins-1;
		return tmp;
	}
	else
		return -1;

}

double gethist(int isize, const double *pvxl, double dmin, double dmax, int numbins, double *pval, double *pbins, double thmin, double thmax)
{
	double dbr= (dmax-dmin)/numbins;
	int i;
	for (i=0;i<numbins;i++)
		pbins[i]=i*dbr+dmin;

	//accumulating
	for (i=0;i<isize;i++)
	{
		int ind=binind1d(dbr, dmin, numbins,pvxl[i], thmin, thmax);
		if (ind>0)
			pval[ind]=pval[ind]+1;
	}

	//calculating entropy
	double *plogv=(double *)mxMalloc(sizeof(double)*(numbins));
	double asum=0.;
	double dt=0.;
	for (i=0;i<numbins;i++)
	{
		asum=asum+fabs(pval[i]);
		if (pval[i]>0)
			plogv[i]=log(pval[i]);
		else
			plogv[i]=0.;
		dt+=plogv[i]*pval[i];
	}
	mxFree(plogv);
	return -(dt/asum)+log(asum)+log(dbr);
}

bool binind2d(int numbins1, int numbins2, double dbr1, double dbr2, 
			double dmin1, double dmin2, float val1, float val2,  
			double thmin1, double thmax1, double thmin2, double thmax2, int* pind)
{
	
	if ((val1>thmin1) && (val1<thmax1) && (val2>thmin2) && (val2<thmax2))
	{
		int indi, indj;
		indi=(int)(0.99999*(val1-dmin1)/dbr1);
		indj=(int)(0.99999*(val2-dmin2)/dbr2);

		if (indi<0)
			indi=0;
		if (indi>=numbins1)
			indi=numbins1-1;
		if (indj<0)
			indj=0;
		if (indj>=numbins2)
			indj=numbins2-1;
		pind[0]=indi;
		pind[1]=indj;
		return true;
	}
	else
		return false;
}

double getjhist(int vsize, const double *pvxlval1, const double *pvxlval2, 
				double dmin1, double dmax1, int numbins1, double thmin1 ,double thmax1, 
				double dmin2, double dmax2, int numbins2, double thmin2 ,double thmax2, 
				double *pval, double *pbins1, double *pbins2, double *pslog10)
{
	double dbr1= (thmax1-thmin1)/numbins1;
	double dbr2= (thmax2-thmin2)/numbins2;
	int i;
	for (i=0;i<numbins1;i++)
		pbins1[i]=i*dbr1+thmin1;
	for (i=0;i<numbins2;i++)
		pbins2[i]=i*dbr2+thmin2;
	//accumulating
	int pind[2];
	pind[0]=0;pind[1]=0;
	for (i=0;i<vsize;i++)
	{
		if (binind2d(numbins1, numbins2, dbr1, dbr2, thmin1,thmin2,  pvxlval1[i],  pvxlval2[i],  
			 thmin1, thmax1, thmin2, thmax2, pind))
		{
			int ind=pind[1]*numbins2+pind[0];
			pval[ind]=pval[ind]+1;
		}	
	}
	//calculating entropy
	int tbins=numbins1*numbins2;
	double *plogv=(double *)mxMalloc(sizeof(double)*(tbins));
	double asum=0.;
	double dt=0.;
	*pslog10=0.;
	for (i=0;i<tbins;i++)
	{
		asum=asum+fabs(pval[i]);
		if (pval[i]>0)
			plogv[i]=log(pval[i]);
		else
			plogv[i]=0.;
		dt+=plogv[i]*pval[i];
		*pslog10+=log10(pval[i]+0.5);
	}
	mxFree(plogv);
	return -(dt/asum)+log(asum)+log(dbr1)+log(dbr2);
}
